# -*- coding: utf-8 -*-
"""
Copyright (C) 2012 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""


import numpy as np
from . import utils
from scipy.constants import c


class SED(object):
    """
    Spectral Energy Distribution with associated information

    This class represents a Spectral Energy Distribution (SED) as constructed
    by pCigale. Such a SED is characterised by:

    - A list of tuples (module name, parametre dictionnary) describing all the
      pCigale modules the SED 'went through'.

    - The wavelengths grid (in nm).

    - An array of luminosity densities (in W/nm) containing the luminosity
      contribution (in positive or negative) of each module the SED 'went
      through'. The first axis corresponds to the index in the contribution
      list (see below) and the second axis corresponds to the wavelength grid.

    - The list of the contributions that are in the above array. This list is
      separated from the list of the modules so that one module can result in
      various contributions to the SED.

    - A dictionnary of arbitrary information associated with the SED.

    """

    def __init__(self):
        self.modules = []
        self.wavelength_grid = None
        self.lumin_contributions = None
        self.contribution_names = []
        self.info = {}

    @property
    def wavelength_grid(self):
        """ Return a copy of the wavelength grid to avoid involuntary
        modifications.
        """
        if self._wavelength_grid is None:
            return None
        else:
            return np.copy(self._wavelength_grid)

    @wavelength_grid.setter
    def wavelength_grid(self, value):
        self._wavelength_grid = value

    @property
    def lumin_contributions(self):
        """ Return a copy of the luminosity contribution table to avoid
        involuntary modifications.
        """
        if self._lumin_contributions is None:
            return None
        else:
            return np.copy(self._lumin_contributions)

    @lumin_contributions.setter
    def lumin_contributions(self, value):
        self._lumin_contributions = value

    @property
    def luminosity(self):
        """
        Return the total luminosity density vector, i.e. the sum of all the
        contributions in W/nm.
        """
        return self.lumin_contributions.sum(0)

    def lambda_fnu(self, redshift=0):
        """
        Return the (redshifted if asked) total Fν flux density vs wavelength
        spectrum of the SED.

        Parametres
        ----------
        redshift : float, default = 0
            If 0 (the default), the flux at 10 pc is computed.

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

        wavelength = utils.redshift_wavelength(self.wavelength_grid, redshift)

        return wavelength, f_nu

    def add_info(self, key, value):
        """
        Add a key / value to the information dictionary

        If the key is present in the dictionnary, it will raise an exception.
        Use this method (instead of direct value assignment ) to avoid
        overriding a yet present information.

        Parametres
        ----------
        key : any immutable
           The key used to retrieve the information.
        value : anything
           The information.

        """
        if key not in self.info:
            self.info[key] = value
        else:
            raise KeyError("The information %s is yet present "
                           "in the SED. " % key)

    def add_component(self, module_name, module_conf, contribution_name,
                      results_wavelengths, results_lumin, infos):
        """
        Add a new module contribution to the SED

        The module name and parametres are added to the module list of the
        SED. If the module adds some information to the SED, it is added to
        the info dictionnary.

        The luminosity contribution of the module is added to the contribution
        table doing an interpolation between the current wavelength grid and
        the grid of the module contribution. During the interpolation,
        everything that is outside of the concerned wavelength domain has its
        luminosity set to 0. Also, the name of the contribution is added to
        the contribution names array.

        Parametres
        ----------
        module_name : string
            Name of the SED creation module used.

        module_conf : dictionary
            The dictionnary containing the module parametres.

        contribution_name : string
            Name of the contribution added. This name is used to retrieve the
            luminosity contribution and allows one module to add more than
            one contribution.

        results_wavelengths : array of floats
            The vector of the wavelengths of the module results (in nm).

        results_lumin : array of floats
            The vector of the Lλ luminosities (in W/nm) of the module results.

        infos : dictionary
            The dictionary of the informations added by the module to the SED.
            If a key of the dictionary is yet present in __info, its value
            will be overwritten.

        """
        self.modules.append((module_name, module_conf))
        self.contribution_names.append(contribution_name)
        for key, value in infos.items():
            self.add_info(key, value)

        # If the SED luminosity table is empty, then there is nothing to
        # compute.
        if self.lumin_contributions is None:
            self.wavelength_grid = np.copy(results_wavelengths)
            self.lumin_contributions = np.array([np.copy(results_lumin)])
        else:
            # Compute the new wavelength grid for the spectrum.
            new_wavelength_grid = utils.best_grid(results_wavelengths,
                                                  self.wavelength_grid)

            # Interpolate each luminosity component to the new wavelength grid
            # setting everything outside the wavelength domain to 0.
            new_lumin_table = None
            for old_lumin in self.lumin_contributions:
                interp_lumin = np.interp(
                    new_wavelength_grid,
                    self.wavelength_grid,
                    old_lumin,
                    right=0,
                    left=0)

                if new_lumin_table is None:
                    new_lumin_table = np.array([interp_lumin])
                else:
                    new_lumin_table = np.vstack((new_lumin_table,
                                                 interp_lumin))

            # Add the new module luminosities
            interp_lumin = np.interp(
                new_wavelength_grid,
                results_wavelengths,
                results_lumin,
                right=0,
                left=0)
            new_lumin_table = np.vstack((new_lumin_table, interp_lumin))

            self.wavelength_grid = new_wavelength_grid
            self.lumin_contributions = new_lumin_table

    def get_lumin_contribution(self, name):
        """Get the luminosity vector of a given contribution

        If the name of the contribution is not unique in the SED, the flux of
        the last one is returned.

        Parametres
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
        return self.lumin_contributions[idx]

    def compute_fnu(self, transmission, lambda_eff, redshift=0):
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

        Parametres
        ----------
        transmission : 2D array of floats
            A numpy 2D array containing the filter response profile
            wavelength[nm] vs transmission).

        lambda_eff : float
            Effective wavelength of the filter in nm.

        redshift : float
            The redshift of the galaxy. If 0, the flux is computed at 10 pc.

        Return
        ------
        fnu : float
            The integrated Fν density in mJy.
        """

        lambda_min = min(transmission[0])
        lambda_max = max(transmission[0])

        if ((min(self.wavelength_grid) > lambda_min) or
                (max(self.wavelength_grid) < lambda_max)):
            f_nu = -99.

        else:
            wavelength = utils.redshift_wavelength(self.wavelength_grid,
                                                   redshift)
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

            # Fν in W/m²/Hz. The 1.e-9 factor is because λ is in nm.
            f_nu = lambda_eff * f_lambda * lambda_eff * 1.e-9 / c

            # Conversion from W/m²/Hz to mJy
            f_nu *= 1.e+29

        return f_nu