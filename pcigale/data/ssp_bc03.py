# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

import numpy as np


class SspBC03(object):
    """Single Stellar Population as defined in Bruzual and Charlot (2003)

    This class holds the data associated with a single stellar population
    (SSP) as defined in Bruzual and Charlot (2003). Compare to the pristine
    SSP, the wavelength are given in nm (rather than Å), the time is given in

    """

    def __init__(self, imf, metallicity, time_grid, wavelength_grid,
                 color_table, lumin_table):
        """Create a new single stellar population as defined in Bruzual and
        Charlot (2003).

        Parameters
        ----------
        imf : string
            Initial mass function (IMF): either 'salp' for Salpeter (1955) or
            'chab' for Chabrier (2003).
        metallicity : float
            The metallicity. Possible values are 0.0001, 0.0004, 0.004, 0.008,
            0.02 (solar metallicity) and 0.05.
        time_grid : array of floats
            The time [Myr] grid used in the colors_table and the lumin_table.
        wavelength_grid : array of floats
            The wavelength [nm] grid used in the lumin_table.
        color_table: 2 axis array of floats
            Array containing information from some of the *.?color tables from
            Bruzual and Charlot (2003) at each time of the time_grid.
                * color_table[0]: Total mass in stars in solar mass
                * color_table[1]: Mass returned to the ISM by evolved stars
                                   in solar mass
                * color_table[2]: rate of H-ionizing photons (s-1)
                * color_table[3]: Amplitude of 4000 Å break (Bruzual 2003)
                * color_table[4]: Amplitude of 4000 Å narrow break (Balogh
                                   et al. 1999)
                * color_table[5]: Amplitude of 4000 Å break (Stoughton
                                   et al. 2002)
                * color_table[6]: Amplitude of Lyman discontinuity
        lumin_table : 2 axis array of floats
            Luminosity density in W/nm. The first axis is the wavelength and
            the second the time (index bases on the wavelength and time grids).
        """

        if imf in ['salp', 'chab']:
            self.imf = imf
        else:
            raise ValueError('IMF must be either sal for Salpeter or '
                             'cha for Chabrier.')
        self.metallicity = metallicity
        self.time_grid = time_grid
        self.wavelength_grid = wavelength_grid
        self.color_table = color_table
        self.lumin_table = lumin_table

    def convolve(self, sfh_time, sfh_sfr, norm=False):
        """Convolve the SSP with a Star Formation History

        Given a SFH (an time grid and the corresponding star formation rate
        SFR), this method convolves the color table and the SSP luminosity
        spectrum along the whole SFR.

        The time grid of the SFH is expected to be ordered and must not run
        beyong 20 Gyr (the maximum time for Bruzual and Charlot 2003 SSP).

        Parameters
        ----------
        sfh_time : array of floats
            Time grid [Myr] of the star formation history. It must be
            increasing and not run beyond 20 Gyr. The SFH will be regrided to
            the SSP time.
        sfh_sfr: array of floats
            Star Formation Rates in Msun/yr at each time of the SFH time grid.
        norm: boolean
            If true, the sfh will be normalised to 1 solar mass produced.

        Returns
        -------
        wavelength : array of floats
            Wavelength grid [nm] for the spectrum
        luminosity : array of floats
            Luminosity density [W/nm] at each wavelength.
        bc03_info: dictionary
            Dictionary containing various information from the *.?color tables:
            - "m_star": Total mass in stars in solar mass
            - "m_gas": Mass returned to the ISM by evolved stars
                       in solar mass
            - "n_ly": rate of H-ionizing photons (s-1)
            - "b_4000": Amplitude of 4000 Å break (Bruzual 2003)
            - "b4_vn": Amplitude of 4000 Å narrow break (Balogh et al. 1999)
            - "b4_sdss" : Amplitude of 4000 Å break (Stoughton et al. 2002)
            - "b_912": Amplitude of Lyman discontinuity

        """
        # We work on a copy of SFH (as we change it)
        sfh_time, sfh_sfr = np.copy((sfh_time, sfh_sfr))

        # Index, in the SSP time grid, of the time nearest to the age of
        # the SFH.
        idx = np.abs(self.time_grid - np.max(sfh_time)).argmin()

        # We regrid the SFH to the time grid of the SSP using a linear
        # interpolation. If the SFH does no start at 0, the first SFR values
        # will be set to 0.
        sfh_sfr = np.interp(self.time_grid[:idx + 1],
                            sfh_time, sfh_sfr,
                            left=0., right=0.)

        # Step between two item in the time grid in Myr
        step = self.time_grid[1] - self.time_grid[0]

        # If needed, we normalise the SFH to 1 solar mass produced.
        if norm:
            sfh_sfr = sfh_sfr / np.trapz(sfh_sfr * 1.e6,
                                         self.time_grid[:idx + 1])

        # As both the SFH and the SSP (limited to the age of the SFH) data now
        # share the same time grid, the convolution is just a matter of
        # reverting one and computing the sum of the one to one product; this
        # is done using the dot product.
        color_table = self.color_table[:, :idx + 1]
        lumin_table = self.lumin_table[:, :idx + 1]

        # The 1.e6 * step is because the SFH is in solar mass per year.
        color_info = 1.e6 * step * np.dot(color_table, sfh_sfr[::-1])
        luminosity = 1.e6 * step * np.dot(lumin_table, sfh_sfr[::-1])

        bc03_info = dict(zip(
            ["m_star", "m_gas", "n_ly", "b_4000", "b4_vn", "b4_sdss", "b_912"],
            color_info
        ))

        return self.wavelength_grid, luminosity, bc03_info
