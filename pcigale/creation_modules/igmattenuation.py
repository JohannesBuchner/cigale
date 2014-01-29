# -*- coding: utf-8 -*-
# Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

"""
Redshift and IGM attenuation module
===================================

This module implements the effect of redshift and the attenuation caused by the
inter-galactic medium. It uses code developed for the Large Synoptic Survey
Telescope. http://dev.lsstcorp.org/trac/

"""

import numpy as np
from collections import OrderedDict
from . import CreationModule
from ..sed import utils
from ..extern.lsst import Sed as lsst


class IGMAtt(CreationModule):
    """Redshift a SED and add IGM attenuation

    This module adds both the effect of redshift and inter-galactic medium
    (IGM) attenuation to a SED object.

    """

    parameter_list = OrderedDict([
        ("redshift", (
            'float',
            "Redshift to apply to the galaxy.",
            0.
        )),
        ("dimming", (
            'boolean',
            "If set to true, the cosmological dimming is applied "
            "to the fluxes.",
            'True'
        )),
        ("rtau", (
            'float',
            "Parameter which scales the tau value at each wavelength.",
            1.
        ))
    ])

    def process(self, sed):
        """Add the redshift + IGM attenuation effect to the SED

        Parameters
        ----------
        sed  : pcigale.sed.SED object
        parameters : dictionary
           Dictionary with the module parameters (redshift and rtau)

        """
        redshift = float(self.parameters["redshift"])
        dimming = (self.parameters['dimming'].lower() == "true")
        rtau = float(self.parameters["rtau"])

        if redshift == 0:
            # If redshift is 0, we do nothing
            pass
        else:
            # We need a lsst.Sed object to use its methods
            lsstSed = lsst.Sed()

            # First, we get the redshifted spectrum of the galaxy
            wavelen, flambda = lsstSed.redshiftSED(redshift,
                                                   dimming,
                                                   sed.wavelength_grid,
                                                   sed.luminosity)

            wavelen, red_l_lambda = lsstSed.addIGMattenuation(
                redshift,
                rtau,
                wavelen=wavelen,
                flambda=flambda
            )
            # We only want to add the redshift + IGM 'effect' to the SED
            # object (even if this has no physical meaning). As the
            # redshifting change the wavelength grid, we must regrid both
            # before subtracting.
            new_wavelen = utils.best_grid(sed.wavelength_grid, wavelen)

            init_l_lambda = np.interp(new_wavelen, sed.wavelength_grid,
                                      sed.luminosity)
            red_l_lambda = np.interp(new_wavelen, wavelen, red_l_lambda)

            igm_effect = red_l_lambda - init_l_lambda

            sed.add_module(self.name, self.parameters)

            sed.add_info("redshift" + self.postfix, redshift)
            sed.add_info('rtau' + self.postfix, rtau)

            sed.add_contribution(
                self.name,
                new_wavelen,
                igm_effect
            )

            sed.redshift = redshift

# CreationModule to be returned by get_module
Module = IGMAtt
