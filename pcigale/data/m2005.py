## -*- coding: utf-8 -*-
# Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

import numpy as np


class M2005(object):
    """Single Stellar Population as defined in Maraston (2005)

    This class holds the data associated with a single stellar population
    (SSP) as defined in Maraston (2005). Compare to the pristine Maraston's
    SSP:

    - The time grid used ranges from 1 My to 13.7 Gyr with 1 My step. This
      excludes the metallicities -2.25 and 0.67 for which the pristine age
      grid starts at 1 Gyr.

    - The stellar masses and the luminosities ( luminosity vs age for a given
      wavelength) were all interpolated to this new grid.

    - The wavelength are given in nm rather than in Å.

    - The luminosities are given in W rather than erg/s.

     """

    def __init__(self, imf, metallicity, time_grid, wavelength_grid,
                 mass_table, spec_table):
        """Create a new single stellar population as defined in Maraston 2005.

        Parameters
        ----------
        imf: string
            Initial mass function (IMF): either 'salp' for single Salpeter
            (1955) or 'krou' for Kroupa (2001).
        metallicity: float
            The metallicity. Possible values are 0.001, 0.01, 0.02 (solar
            metallicity) and 0.04.
        time_grid: array of floats
            The time [Myr] grid used in the mass_table and the spec_table.
        wavelength_grid: array of floats
            The wavelength [nm] grid used in the spec_table.
        mass_table: (6, n) array of floats
            The 2D table giving the various stellar masses at a given age. The
            first axis is the king of mass, the second is the age based on the
            time_grid.
                * mass_table[0]: total star mass
                * mass_table[1]: alive star mass
                * mass_table[2]: white dwarf star mass
                * mass_table[3]: neutrino star mass
                * mass_table[4]: black hole star mass
                * mass_table[5]: mass in the turn off
        spec_table: (2, n) array of floats
            The 2D table giving the luminosity density [W/nm] at various time.
            The first axis is the age, base on the time_grid, the second is the
            wavelength, base on the wavelength_grid.

        """

        if imf in ['salp', 'krou']:
            self.imf = imf
        else:
            raise ValueError('IMF must be either salp for Salpeter or '
                             'krou for Krupa.')
        self.metallicity = metallicity
        self.time_grid = time_grid
        self.wavelength_grid = wavelength_grid
        self.mass_table = mass_table
        self.spec_table = spec_table

    def convolve(self, sfh_time, sfh_sfr, norm=False):
        """Convolve the SSP with a Star Formation History

        Given a SFH (an time grid and the corresponding star formation rate
        SFR), this method convolves the mass distribution and the SSP spectrum
        along the whole SFR.

        The time grid of the SFH is expected to be ordered and must not run
        beyond 13.7 Gyr (the maximum time for Maraston 2005 SSP).

        Parameters
        ----------
        sfh_time: array of floats
            Time grid [Myr)] of the star formation history. It must be
            increasing and not run beyond 13.7 Gyr. As the SFH will be
            regrided to the SSP time grid, it is better to have a SFH age grid
            compatible, i.e. with a precision limited to 1 Myr.
        sfh_sfr: array of floats
            Star Formation Rates in Msun/yr at each time of the SFH time grid.
        norm: boolean
            If true, the sfh will be normalised to 1 solar mass produced.

        Returns
        -------
        masses, spectra: array of floats, array of floats
            masses is an array of floats or and array of arrays of floats (the
            second axe is then the age):
                 - masses[0]: total stellar mass(es)
                 - masses[1]: alive star mass(es)
                 - masses[2]: white dwarf mass(es)
                 - masses[3]: neutrino star mass(es)
                 - masses[4]: black hole mass(es)
                 - masses[5]: mass(es) in the turn-off
            spectra holds:
                - spectra[0]: wavelengths in nm
                - spectra[1]: luminosity in W/nm

        """
        # We work on a copy of SFH (as we change it)
        sfh_time, sfh_sfr = np.copy((sfh_time, sfh_sfr))

        # Step between two item in the time grid in Myr
        step = self.time_grid[1] - self.time_grid[0]

        # Number of step to go to the age of the SFH on the SSP age grid.
        nb_steps = 1 + np.round((np.max(sfh_time) - self.time_grid[0]) / step)
        # We regrid the SFH to the time grid of the SSP using a linear
        # interpolation. If the SFH does no start at 0, the first SFR values
        # will be set to 0.
        sfh_sfr = np.interp(self.time_grid[:nb_steps],
                            sfh_time, sfh_sfr,
                            left=0., right=0.)

        # If needed, we normalise the SFH to 1 solar mass produced.
        if norm:
            sfh_sfr = sfh_sfr / np.trapz(sfh_sfr * 1.e6,
                                         self.time_grid[:nb_steps])

        # As both the SFH and the SSP (limited to the age of the SFH) data now
        # share the same time grid, the convolution is just a matter of
        # reverting one and computing the sum of the one to one product; this
        # is done using the dot product.
        mass_table = self.mass_table[:,:nb_steps]
        spec_table = self.spec_table[:,:nb_steps]

        # The 1.e6 * step is because the SFH is in solar mass per year.
        masses = 1.e6 * step * np.dot(mass_table, sfh_sfr[::-1])
        spectra = 1.e6 * step * np.dot(spec_table, sfh_sfr[::-1])

        return masses, spectra
