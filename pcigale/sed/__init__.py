# -*- coding: utf-8 -*-
"""Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

This class represents a Spectral Energy Distribution (SED) as used by pcigale.
Such SED is characterised by:

- sfh: a tuple (time [Myr], Star Formation Rate [Msun/yr]) representing the
  Star Formation History of the galaxy.

- modules: a list of tuples (module name, parameter dictionary) containing all
  the pcigale modules the SED 'went through'.

- wavelength_grid: the grid of wavelengths [nm] used for the luminosities.

- contribution_names: the list of the names of the luminosity contributions
  making part of the SED.

- luminosities: a two axis numpy array containing all the luminosity density
  [W/nm] contributions to the SED. The index in the first axis corresponds to
  the contribution (in the contribution_names list) and the index of the
  second axis corresponds to the wavelength in the wavelength grid.

- lines: a dictionnary containing the emission lines associated with the SED.
  A dictionnary is used to allow the storage of various sets of lines. The
  lines are stored in lists of tuples (wavelength [nm], luminosity [W], width
  [km.s-1]).

- info: a dictionnary containing various information about the SED.

- mass_proportional_info: the list of keys in the info dictionnary whose value
  is proportional to the galaxy mass.

"""


import numpy as np
from . import utils
from itertools import chain
from scipy.constants import c
from scipy.interpolate import interp1d


class SED(object):
    """Spectral Energy Distribution with associated information
    """

    def __init__(self, sfh=None):
        """Create a new SED

        Parameters
        ----------
        sfh : (numpy.array, numpy.array)
            Star Formation History: tuple of two numpy array, the first is the
            time in Myr and the second is the Star Formation Rate in Msun/yr.
            If no SFH is given, it's set to None.

        """
        self.sfh = sfh
        self.modules = []
        self.wavelength_grid = None
        self.contribution_names = []
        self.luminosities = None
        self.lines = {}
        self.info = {}
        self.mass_proportional_info = []

    @property
    def sfh(self):
        """Return a copy of the star formation history
        """
        if self._sfh is None:
            return None
        else:
            return np.copy(self._sfh)

    @sfh.setter
    def sfh(self, value):
        self._sfh = value

    @property
    def wavelength_grid(self):
        """ Return a copy of the wavelength grid
        """
        if self._wavelength_grid is None:
            return None
        else:
            return np.copy(self._wavelength_grid)

    @wavelength_grid.setter
    def wavelength_grid(self, value):
        self._wavelength_grid = value

    @property
    def luminosities(self):
        """ Return a copy of the luminosity contributions
        """
        if self._luminosities is None:
            return None
        else:
            return np.copy(self._luminosities)

    @luminosities.setter
    def luminosities(self, value):
        self._luminosities = value

    @property
    def luminosity(self):
        """Total luminosity of the SED

        Return the total luminosity density vector, i.e. the sum of all the
        contributions in W/nm.
        """
        if self._luminosities is None:
            return None
        else:
            return self._luminosities.sum(0)

    def lambda_fnu(self, redshift=0, redshift_spectrum=False):
        """
        Return the (redshifted if asked) total Fν flux density vs wavelength
        spectrum of the SED.

        Parameters
        ----------
        redshift : float, default = 0
            If 0 (the default), the flux at 10 pc is computed.
        redshift_spectrum : boolean, default = None
            If true, the spectrum will be redshifted before computing the
            flux. The default is False because we generally use a specific
            module to apply the redshift.

        Returns
        -------
        wavelength, f_nu : tuple of array of floats
            The wavelength is in nm and the Fν is in mJy.

        """
        # Fλ flux density in W/m²/nm
        f_lambda = utils.luminosity_to_flux(self.luminosity, redshift)

        # The 1.e29 factor is to convert from W/m²/Hz to mJy
        f_nu = (f_lambda * self.wavelength_grid
                * 1.e-9 * self.wavelength_grid / c
                * 1.e29)

        if redshift_spectrum:
            wavelength = utils.redshift_wavelength(self.wavelength_grid,
                                                   redshift)
        else:
            wavelength = np.copy(self.wavelength_grid)

        return wavelength, f_nu

    def add_info(self, key, value, mass_proportional=False, force=False):
        """
        Add a key / value to the information dictionary

        If the key is present in the dictionary, it will raise an exception.
        Use this method (instead of direct value assignment ) to avoid
        overriding a yet present information.

        Parameters
        ----------
        key : any immutable
           The key used to retrieve the information.
        value : anything
           The information.
        mass_proportional : boolean
           If True, the added variable is set as proportional to the
           mass.
        force : boolean
           If false (default), adding a key that yet exists in the info
           dictionary will raise an error. If true, doing this will update
           the associated value.

        """
        if (key not in self.info) or force:
            self.info[key] = value
            if mass_proportional:
                self.mass_proportional_info.append(key)
        else:
            raise KeyError("The information %s is yet present "
                           "in the SED. " % key)

    def add_module(self, module_name, module_conf):
        """Add a new module information to the SED.

        Parameters
        ----------
        module_name : string
            Name of the module. This name can be suffixed with anything
            using a dot.
        module_conf : dictionary
            Dictionary containing the module parameters.

        TODO: Complete the parameter dictionary with the default values from
              the module if they are not present.

        """
        self.modules.append((module_name, module_conf))

    def add_contribution(self, contribution_name, results_wavelengths,
                         results_lumin):
        """
        Add a new luminosity contribution to the SED.

        The luminosity contribution of the module is added to the contribution
        table doing an interpolation between the current wavelength grid and
        the grid of the module contribution. During the interpolation,
        everything that is outside of the concerned wavelength domain has its
        luminosity set to 0. Also, the name of the contribution is added to
        the contribution names array.

        Parameters
        ----------
        contribution_name : string
            Name of the contribution added. This name is used to retrieve the
            luminosity contribution and allows one module to add more than
            one contribution.

        results_wavelengths : array of floats
            The vector of the wavelengths of the module results (in nm).

        results_lumin : array of floats
            The vector of the Lλ luminosities (in W/nm) of the module results.

        """
        self.contribution_names.append(contribution_name)

        # If the SED luminosity table is empty, then there is nothing to
        # compute.
        if self.luminosities is None:
            self.wavelength_grid = np.copy(results_wavelengths)
            self.luminosities = np.copy(results_lumin)
        else:
            # Compute the new wavelength grid for the spectrum.
            new_wavelength_grid = utils.best_grid(results_wavelengths,
                                                  self.wavelength_grid)

            # Interpolate each luminosity component to the new wavelength grid
            # setting everything outside the wavelength domain to 0.
            new_luminosities = interp1d(self.wavelength_grid,
                                        self.luminosities,
                                        bounds_error=False,
                                        fill_value=0.)(new_wavelength_grid)

            # Interpolate the added luminosity array to the new wavelength
            # grid
            interp_lumin = np.interp(
                new_wavelength_grid,
                results_wavelengths,
                results_lumin,
                right=0,
                left=0)

            self.wavelength_grid = new_wavelength_grid
            self.luminosities = np.vstack((new_luminosities, interp_lumin))

    def get_lumin_contribution(self, name):
        """Get the luminosity vector of a given contribution

        If the name of the contribution is not unique in the SED, the flux of
        the last one is returned.

        Parameters
        ----------
        name : string
            Name of the contribution

        Returns
        -------
        luminosities : array of floats
            Vector of the luminosity density contribution based on the SED
            wavelength grid.

        """
        # Find the index of the _last_ name element
        idx = (len(self.contribution_names) - 1
               - self.contribution_names[::-1].index(name))
        return self.luminosities[idx]

    def compute_fnu(self, transmission, lambda_eff,
                    redshift=0, redshift_spectrum=False):
        """
        Compute the Fν flux density corresponding the filter which
        transmission is given.

        As the SED stores the Lλ luminosity density, we first compute the Fλ
        flux density. Fλ is the integration of the Lλ luminosity multiplied by
        the filter transmission, normalised to this transmission and corrected
        by the luminosity distance of the source. This is done by the
        pcigale.sed.utils.luminosity_to_flux function.

        Fλ = luminosity_to_flux( integ( LλT(λ)dλ ) / integ( T(λ)dλ ) )

        Fλ is in W/m²/nm. At redshift 0, the flux is computed at 10 pc. Then,
        to compute Fν, we make the approximation:

        Fν = λeff / c . λeff . Fλ

        Fν is computed in W/m²/Hz and then converted to mJy.

        If the SED spectrum does not cover all the filter response table,
        -99 is returned.

        Parameters
        ----------
        transmission : 2D array of floats
            A numpy 2D array containing the filter response profile
            wavelength[nm] vs transmission).

        lambda_eff : float
            Effective wavelength of the filter in nm.

        redshift : float
            The redshift of the galaxy. If 0, the flux is computed at 10 pc.

        redshift_spectrum : boolean
            If true, the spectrum will be redshifted before computing the
            flux. The default is False because we generally use a specific
            module to apply the redshift.

        Return
        ------
        fnu : float
            The integrated Fν density in mJy.
        """

        lambda_min = min(transmission[0])
        lambda_max = max(transmission[0])

        # FIXME Shouldn't it be the reverse
        if ((min(self.wavelength_grid) > lambda_min) or
                (max(self.wavelength_grid) < lambda_max)):
            f_nu = -99.

        else:
            if redshift_spectrum:
                wavelength = utils.redshift_wavelength(self.wavelength_grid,
                                                       redshift)
            else:
                wavelength = np.copy(self.wavelength_grid)

            l_lambda = self.luminosity

            # We regrid both spectrum and filter to the best wavelength grid
            # to avoid interpolating a high wavelength density curve to a low
            # density one. Also, we limit the work wavelength domain to the
            # filter one, taking care the presence of λmin and λman in the
            # used wavelength grid.
            wavelength_r = utils.best_grid(wavelength, transmission[0])
            if lambda_min not in wavelength_r:
                wavelength_r.append(lambda_min)
            if lambda_max not in wavelength_r:
                wavelength_r.append(lambda_max)
            wavelength_r.sort()
            wavelength_r = wavelength_r[wavelength_r <= lambda_max]
            wavelength_r = wavelength_r[wavelength_r >= lambda_min]

            l_lambda_r = np.interp(wavelength_r, wavelength, l_lambda)
            transmission_r = np.interp(wavelength_r, transmission[0],
                                       transmission[1])

            # TODO: Can we avoid to normalise as the filter transmission is
            # yet normalised?
            f_lambda = utils.luminosity_to_flux(
                (np.trapz(transmission_r * l_lambda_r, wavelength_r) /
                 np.trapz(transmission_r, wavelength_r)),
                redshift
            )

            # Add the Fλ fluxes from the spectral lines.
            for line in chain(*self.lines.values()):
                if (line[0] >= lambda_min) and (line[0] <= lambda_max):
                    pass  # TODO write the code

            # Fν in W/m²/Hz. The 1.e-9 factor is because λ is in nm.
            f_nu = lambda_eff * f_lambda * lambda_eff * 1.e-9 / c

            # Conversion from W/m²/Hz to mJy
            f_nu *= 1.e+29

        return f_nu
