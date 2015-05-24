# -*- coding: utf-8 -*-
# Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

"""
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

- info: a dictionary containing various information about the SED.

- mass_proportional_info: the list of keys in the info dictionary whose value
  is proportional to the galaxy mass.

"""

import numpy as np
from collections import OrderedDict
from . import utils
from .io.vo import save_sed_to_vo
from scipy.constants import c, parsec
from scipy.interpolate import interp1d
from ..data import Database


# Time lapse used to compute the average star formation rate. We use a
# constant to keep it easily changeable for advanced user while limiting the
# number of parameters. The value is in Myr.
AV_LAPSE = 100


class SED(object):
    """Spectral Energy Distribution with associated information
    """

    def __init__(self, sfh=None):
        """Create a new SED

        Parameters
        ----------
        sfh: (numpy.array, numpy.array)
            Star Formation History: tuple of two numpy array, the first is the
            time in Myr and the second is the Star Formation Rate in Msun/yr.
            If no SFH is given, it's set to None.

        """
        self.sfh = sfh
        self.modules = []
        self.wavelength_grid = None
        self.contribution_names = []
        self.luminosity = None
        self.luminosities = None
        self.info = OrderedDict()
        self.mass_proportional_info = []
        self.cache_filters = {}

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

        if value:
            sfh_time, sfh_sfr = value
            sfh_age = sfh_time[-1] - sfh_time
            self._sfh = value
            self.add_info("sfh.sfr", sfh_sfr[-1], True)
            self.add_info("sfh.sfr10Myrs", np.mean(sfh_sfr[-10:]), True)
            self.add_info("sfh.sfr100Myrs", np.mean(sfh_sfr[-100:]), True)
            self.add_info("sfh.age", sfh_time[-1], False)

    @property
    def fnu(self):
        """Total Fν flux density of the SED

        Return the total Fν density vector, i.e the total luminosity converted
        to Fν flux in mJy.
        """

        # Fλ flux density in W/m²/nm
        f_lambda = utils.luminosity_to_flux(self.luminosity,
                                            self.info
                                            ['universe.luminosity_distance'])

        # Fν flux density in mJy
        f_nu = utils.lambda_flambda_to_fnu(self.wavelength_grid, f_lambda)

        return f_nu

    def add_info(self, key, value, mass_proportional=False, force=False):
        """
        Add a key / value to the information dictionary

        If the key is present in the dictionary, it will raise an exception.
        Use this method (instead of direct value assignment ) to avoid
        overriding an already present information.

        Parameters
        ----------
        key: any immutable
           The key used to retrieve the information.
        value: anything
           The information.
        mass_proportional: boolean
           If True, the added variable is set as proportional to the
           mass.
        force: boolean
           If false (default), adding a key that already exists in the info
           dictionary will raise an error. If true, doing this will update
           the associated value.

        """
        if (key not in self.info) or force:
            self.info[key] = value
            if mass_proportional:
                self.mass_proportional_info.append(key)
        else:
            raise KeyError("The information %s is already present "
                           "in the SED. " % key)

    def add_module(self, module_name, module_conf):
        """Add a new module information to the SED.

        Parameters
        ----------
        module_name: string
            Name of the module. This name can be suffixed with anything
            using a dot.
        module_conf: dictionary
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
        contribution_name: string
            Name of the contribution added. This name is used to retrieve the
            luminosity contribution and allows one module to add more than
            one contribution.

        results_wavelengths: array of floats
            The vector of the wavelengths of the module results (in nm).

        results_lumin: array of floats
            The vector of the Lλ luminosities (in W/nm) of the module results.

        """
        self.contribution_names.append(contribution_name)

        # If the SED luminosity table is empty, then there is nothing to
        # compute.
        if self.luminosity is None:
            self.wavelength_grid = results_wavelengths.copy()
            self.luminosity = results_lumin.copy()
            self.luminosities = results_lumin.copy()
        else:
            # If the added luminosity contribution changes the SED wavelength
            # grid, we interpolate everything on a common wavelength grid.
            if (results_wavelengths.size != self.wavelength_grid.size or
                    not np.all(results_wavelengths == self.wavelength_grid)):
                # Compute the new wavelength grid for the spectrum.
                new_wavelength_grid = utils.best_grid(results_wavelengths,
                                                      self.wavelength_grid)

                # Interpolate each luminosity component to the new wavelength
                # grid setting everything outside the wavelength domain to 0.
                new_luminosities = interp1d(self.wavelength_grid,
                                            self.luminosities,
                                            bounds_error=False,
                                            assume_sorted=True,
                                            fill_value=0.)(new_wavelength_grid)

                # Interpolate the added luminosity array to the new wavelength
                # grid
                interp_lumin = np.interp(new_wavelength_grid,
                                         results_wavelengths, results_lumin,
                                         left=0., right=0.)

                self.wavelength_grid = new_wavelength_grid
                self.luminosities = np.vstack((new_luminosities, interp_lumin))
                self.luminosity = self.luminosities.sum(0)
            else:
                self.luminosities = np.vstack((self.luminosities,
                                               results_lumin))
                self.luminosity += results_lumin

    def get_lumin_contribution(self, name):
        """Get the luminosity vector of a given contribution

        If the name of the contribution is not unique in the SED, the flux of
        the last one is returned.

        Parameters
        ----------
        name: string
            Name of the contribution

        Returns
        -------
        luminosities: array of floats
            Vector of the luminosity density contribution based on the SED
            wavelength grid.

        """
        # Find the index of the _last_ name element
        idx = (len(self.contribution_names) - 1 -
               self.contribution_names[::-1].index(name))
        return self.luminosities[idx]

    def compute_fnu(self, filter_name):
        """
        Compute the Fν flux density in a given filter

        As the SED stores the Lλ luminosity density, we first compute the Fλ
        flux density. Fλ is the integration of the Lλ luminosity multiplied by
        the filter transmission, normalised to this transmission and corrected
        by the luminosity distance of the source.

        Fλ is in W/m²/nm. At redshift 0, the flux is computed at 10 pc. Then,
        to compute Fν, we make the approximation:

        Fν = λeff / c . λeff . Fλ

        Fν is computed in W/m²/Hz and then converted to mJy.

        If the SED spectrum does not cover all the filter response table,
        -99 is returned.

        Parameters
        ----------
        filter_name: string
            Name of the filter to integrate into. It must be presnt in the
            database.

        Return
        ------
        fnu: float
            The integrated Fν density in mJy.
        """

        wavelength = self.wavelength_grid

        # First we try to fetch the filter's wavelength, transmission and
        # effective wavelength from the cache. The two keys are the size of the
        # spectrum wavelength grid and the name of the filter. The first key is
        # necessary because different spectra may have different sampling. To
        # avoid having the resample the filter every time on the optimal grid
        # (spectrum+filter), we store the resampled filter. That way we only
        # have to resample to spectrum.
        key = (wavelength.size, filter_name)
        if key in list(self.cache_filters.keys()):
            wavelength_r, transmission_r, lambda_eff = self.cache_filters[key]
        else:
            with Database() as db:
                filter_ = db.get_filter(filter_name)
            trans_table = filter_.trans_table
            lambda_eff = filter_.effective_wavelength
            lambda_min = filter_.trans_table[0][0]
            lambda_max = filter_.trans_table[0][-1]

            # Test if the filter covers all the spectrum extent. If not then
            # the flux is not defined
            if ((wavelength[0] > lambda_min) or (wavelength[-1] < lambda_max)):
                return -99.

            # We regrid both spectrum and filter to the best wavelength grid
            # to avoid interpolating a high wavelength density curve to a low
            # density one. Also, we limit the work wavelength domain to the
            # filter one.
            w = np.where((wavelength >= lambda_min) &
                         (wavelength <= lambda_max))
            wavelength_r = utils.best_grid(wavelength[w], trans_table[0])
            transmission_r = np.interp(wavelength_r, trans_table[0],
                                       trans_table[1])

            self.cache_filters[key] = (wavelength_r, transmission_r,
                                       lambda_eff)

        l_lambda_r = np.interp(wavelength_r, wavelength, self.luminosity)

        if 'universe.luminosity_distance' in self.info.keys():
            dist = self.info['universe.luminosity_distance']
        else:
            dist = 10. * parsec

        f_lambda = utils.luminosity_to_flux(
            np.trapz(transmission_r * l_lambda_r, wavelength_r),
            dist)

        # Return Fν in mJy. The 1e-9 factor is because λ is in nm and 1e29 for
        # convert from W/m²/Hz to mJy.
        return lambda_eff * lambda_eff * f_lambda * 1e-9 / c * 1e29

    def to_votable(self, filename, mass=1.):
        """
        Save the SED to a VO-table file

        Parameters
        ----------
        filename: string
            Name of the VO-table file
        mass: float
            Galaxy mass in solar mass. When need, the saved data will be
            multiplied by this mass.

        """
        save_sed_to_vo(self, filename, mass)

    def copy(self):
        sed = SED()
        sed.sfh = (self._sfh[0].copy(), self._sfh[1].copy())
        sed.modules = self.modules[:]
        if self.wavelength_grid is not None:
            sed.wavelength_grid = self.wavelength_grid.copy()
            sed.luminosity = self.luminosity.copy()
            sed.luminosities = self.luminosities.copy()
        else:
            sed.wavelength_grid = None
            sed.luminosity = None
            sed.luminosities = None
        sed.contribution_names = self.contribution_names[:]
        sed.info = self.info.copy()
        sed.mass_proportional_info = self.mass_proportional_info[:]
        sed.cache_filters = self.cache_filters

        return sed
