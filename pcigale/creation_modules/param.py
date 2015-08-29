# -*- coding: utf-8 -*-
# Copyright (C) 2015 Laboratoire d'Astrophysique de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Denis Burgarella

"""
Module that estimates other parameters, e.g., UV slope, Lick indices, etc.
==========================================================================

This module estimates additional parameters of interest
close to the observation, e.g., the ultraviolet slope (beta), the rest-frame
far-ultraviolet luminosity, any type of indices (Lick), etc.
This module will be the last one (after redshifting) to account
for all the physical processes at play to build the received total emission.

"""

from collections import OrderedDict
import numpy as np
from . import CreationModule


class Param(CreationModule):
    """Other parameters

    This module does not need any input.
    It computes some parameters from the models:
    beta_UV, Lum_UV, etc.

    """

    parameter_list = OrderedDict([
        ("OK", (
            "boolean",
            "Are you here?.",
            "True"
        ))
    ])

    def process(self, sed):
        """Computes the parameters for each model.

        Parameters
        ----------
        sed: pcigale.sed.SED object

        """
        # Retrieve the final computed SED using all the previous modules
        # including the IGM and the redshifting. In other words,
        # this module must be the last one. Note that it does require
        # an SFH and an SSP module but nothing else (except redshifting)

        redshift = sed.info['redshift']
        # Wavelengths are in nanometers.
        wavelength = sed.wavelength_grid
        # Luminosity is is W/nm.
        luminosity = sed.luminosity

        # Attenuated (observed) UV slopes beta
        # as defined in Calzetti et al. (1994, ApJ 429, 582, Tab. 2)
        # that excludes the 217.5 nm bump wavelength range and other spectral features

        w_calz94 = np.where((wavelength >= 126.8 * (1. + redshift)) &
                            (wavelength <= 128.4 * (1. + redshift)) |
                            (wavelength >= 130.9 * (1. + redshift)) &
                            (wavelength <= 131.6 * (1. + redshift)) |
                            (wavelength >= 134.2 * (1. + redshift)) &
                            (wavelength <= 137.1 * (1. + redshift)) |
                            (wavelength >= 140.7 * (1. + redshift)) &
                            (wavelength <= 151.5 * (1. + redshift)) |
                            (wavelength >= 156.2 * (1. + redshift)) &
                            (wavelength <= 158.3 * (1. + redshift)) |
                            (wavelength >= 167.7 * (1. + redshift)) &
                            (wavelength <= 174.0 * (1. + redshift)) |
                            (wavelength >= 176.0 * (1. + redshift)) &
                            (wavelength <= 183.3 * (1. + redshift)) |
                            (wavelength >= 186.6 * (1. + redshift)) &
                            (wavelength <= 189.0 * (1. + redshift)) |
                            (wavelength >= 193.0 * (1. + redshift)) &
                            (wavelength <= 195.0 * (1. + redshift)) |
                            (wavelength >= 240.0 * (1. + redshift)) &
                            (wavelength <= 258.0 * (1. + redshift)))

        # Attenuated (observed) FUV luminosity L_FUV, in the GALEX FUV band,
        # i.e., lambda_eff = 152.8 nm and effective bandwidth = 11.4 nm

        w_FUV = np.where((wavelength >= (152.8 - 11.4) * (1. + redshift)) &
                         (wavelength <= (152.8 + 11.4) * (1. + redshift)) )

        # Strength of the D_4000 break using Balogh et al. (1999, ApJ 527, 54), i.e.,
        # ratio of the flux in the red continuum to that in the blue continuum:
        # Blue continuum: 385.0-395.0 nm& red continuum: 410.0-410.0 nm.

        w_D4000blue = np.where((wavelength >= 385.0 * (1. + redshift)) &
                               (wavelength <= 395.0 * (1. + redshift)))
        w_D4000red  = np.where((wavelength >= 400.0 * (1. + redshift)) &
                               (wavelength <= 410.0 * (1. + redshift)))

        regression_calz94 = np.polyfit(np.log10(10.*wavelength[w_calz94]),
                                       np.log10(1e7/10.*luminosity[w_calz94]), 1)
        beta_calz94 = regression_calz94[0]

        L_FUV = 152.8*(1. + redshift)*np.mean(luminosity[w_FUV])

        D_4000 = np.mean(luminosity[w_D4000red]) / np.mean(luminosity[w_D4000blue])

        sed.add_info("param.beta_calz94", beta_calz94)
        sed.add_info("param.FUV_luminosity", L_FUV, True)
        sed.add_info("param.D_4000", D_4000)

# CreationModule to be returned by get_module
Module = Param
