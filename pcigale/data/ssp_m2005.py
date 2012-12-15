# -*- coding: utf-8 -*-
"""
Copyright (C) 2012 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""

import numpy as np
from scipy import interpolate


class SspM2005(object):
    """Single Stellar Population as defined in Maraston (2005)

    This class holds the data associated with a single stellar population
    (SSP) as defined in Maraston (2005). Compare to the pristine Maraston's
    SSP:

    - The age grid used ranges from 1 My to 13.7 Gyr with 1 My step. This
      excludes the metallicities -2.25 and 0.67 for which the pristine age
      grid starts at 1 Gyr.

    - The stellar masses and the luminosities ( luminosity vs age for a given
      wavelength) were all interpolated to this new grid.

    - The wavelength are given in nm rather than in Å.

    - The luminosities are given in W rather than erg/s.

     """

    def __init__(self, imf, metallicity, time_grid, wavelength_grid,
                 mass_table, spec_table):
        """Create a new single stellar population as defined in Maraston 2005.

        Parametres
        ----------
        imf : string
            Initial mass function (IMF): either 'ss' for single Salpeter
            (1955) or 'kr' for Kroupa (2001).
        metallicity : float
            The metallicity [Z/H] defined as the abundance of heavy elements
            with respect to hydrogen, normalised to the solar values:
            [Z/H] = Log10(Z/Zsun) - Log10(H/Hsun). The possible values are:
                * +0.35 (corresponding to 2.0 Zsun)
                * +0.00 (corresponding to 1.0 Zsun)
                * -0.33 (corresponding to 0.5 Zsun)
                * -1.35 (corresponding to 1/50 Zsun)
        time_grid : array of floats
            The age [Gyr] grid used in the mass_table and the spec_table.
        wavelength_grid : array of floats
            The wavelength [nm] grid used in the spec_table.
        mass_table : (6, n) array of floats
            The 2D table giving the various stellar masses at a given age. The
            first axis is the king of mass, the second is the age based on the
            time_grid.
                * mass_table[0]: total star mass
                * mass_table[1]: alive star mass
                * mass_table[2]: white dwarf star mass
                * mass_table[3]: neutrino star mass
                * mass_table[4]: black hole star mass
                * mass_table[5]: mass in the turn off
        spec_table : (2, n) array of floats
            The 2D table giving the luminosity density [W/nm] at various age.
            The first axis is the age, base on the time_grid, the second is the
            wavelength, base on the wavelength_grid.

        """

        if imf in ['ss', 'kr']:
            self.imf = imf
        else:
            raise ValueError('IMF must be either ss for Salpeter or '
                             'kr for Krupa.')
        self.metallicity = metallicity
        self.time_grid = time_grid
        self.wavelength_grid = wavelength_grid
        self.mass_table = mass_table
        self.spec_table = spec_table

    def convolve(self, sfr, age, norm=False):
        """
        Given a Star Formation History (SFH) and an age, this method convolves
        the mass distribution and the spectrum at a given age.

        If 'age' is an array, the full convolution is computed and the method
        returns the linear interpolation of the convolution result at each
        age. This can take a few minutes.

        If 'age' is a float, the computation is a lot faster, based on the
        fact that the SFR is given at each age of the SSP age grid.

        Parametres
        ----------
        sfr : array of floats
            Star Formation Rates in Msun/y for each age of the SSP age grid.
        age : float or array of floats
            Age(s) in Gyr we want the mass distribution at. If an array is
            given, the various outputs will be arrays of the same length.
        norm: boolean
            If true, the sfr will be normalised to 1 solar mass produced at
            'age' (or max(age) if it's an array).

        Returns
        -------
        masses, spectra : array of floats, array of floats
            masses is an array of floats or and array of arrays of floats (the
            second axe is then the age):
                 - masses[0] : total stellar mass(es)
                 - masses[1] : alive star mass(es)
                 - masses[2] : white dwarf mass(es)
                 - masses[3] : neutrino star mass(es)
                 - masses[4] : black hole mass(es)
                 - masses[5] : mass(es) in the turn-off
            spectra holds:
                - spectra[0] : wavelengths in nm
                - spectra[1] : luminosity in W/nm

        """
        # We work on a copy of SFR (as we change it)
        sfr = np.copy(sfr)

        # The SFR must be on the same age grid as the SSP.
        if not len(sfr) == len(self.time_grid):
            raise ValueError("The star formation rate must be base"
                             "on the same grid than the SSP.")

        # Check if the age parametre is a unique float or an array.
        try:
            age = float(age)
            isAgeUnique = True
        except:
            age = np.array(age, dtype='float')
            isAgeUnique = False

        # Step between two item in the age grid in Gyr
        step = self.time_grid[1] - self.time_grid[0]

        if isAgeUnique:
            # This is the fast convolution technique
            # 1. We find the index of the nearest element to the given age in
            # the age grid.
            idx = np.abs(self.time_grid - age).argmin()

            # 2. If needed, we normalise the SFR to 1 solar mass formed at
            # 'age'.
            if norm:
                sfr = sfr / (1.e6 * sfr[:idx + 1].sum())

            # 3. We must convolve the mass evolution array with the SFR at a
            # given time (i.e. a given index). As both tables share the same
            # age grid, it's just a matter of slicing the arrays to the given
            # index, reverting one and computing the sum of the one to one
            # product. This is done using the dot product.
            # The 1.e9 * step is because the SFR is in solar mass per year.
            masses = 1.e9 * step * np.dot(self.mass_table[:, :idx + 1],
                                          sfr[:idx + 1][::-1])

            #4.  We do the same thing for the spectre.
            spectra = 1.e9 * step * np.dot(self.spec_table[:, :idx + 1],
                                           sfr[:idx + 1][::-1])

        else:
            # If the age parametre is an array, we do the full convolution.
            # TODO: We can stop the convolution at the youngest age.

            # If needed, we normalise the SFR to 1 solar mass formed at the
            # end of the 'age' lapse.
            if norm:
                idx = np.abs(self.time_grid - max(age)).argmin()
                sfr = sfr / (1.e6 * sfr[:idx + 1].sum())

            # We convolve the mass table with the SFR for each kind of mass.
            # Because of the way numpy full convolution work, we must take in
            # the result the first slice with the length of the age grid. The
            # SFR is multiplicated by step*1e9 to convert it in Msun/Gyr
            conv_masses = [
                np.convolve(table, sfr * step * 1.e9)[0:len(self.time_grid)]
                for table in self.mass_table]
            conv_masses = np.array(conv_masses)

            # We then interpolate the convolved mass table to the given age
            # and return it.
            masses = interpolate.interp1d(self.time_grid, conv_masses)(age)

            # We convolve the spectrum table with the SFR.
            conv_spectra = [
                np.convolve(table, sfr * step * 1e9)[0:len(self.time_grid)]
                for table in self.spec_table]
            conv_spectra = np.array(conv_spectra)
            spectra = interpolate.interp1d(self.time_grid, conv_spectra)(age)

        return masses, spectra
