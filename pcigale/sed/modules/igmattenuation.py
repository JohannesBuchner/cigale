#-*- coding: utf-8 -*-
"""
Copyright (C) 2012 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""

import numpy as np
from . import common
from .. import utils
from ...extern.lsst import Sed as lsst


class Module(common.SEDCreationModule):
    """Redshift a SED and add IGM attenuation

    This module adds both the effect of redshift and inter-galactic medium
    (IGM) attenuation to a SED object. It is based on the code from the Large
    Synoptic Survey Telescope. http://dev.lsstcorp.org/trac/

    """

    parametre_list = {
        "redshift": (
            'float',
            None,
            "Redshift to apply to the galaxy.",
            0.
        ),
        "dimming": (
            'boolean',
            None,
            "If set to true, the cosmological dimming is applied "
            "to the fluxes.",
            True
        ),
        "rtau": (
            'float',
            None,
            "Parametre which scale the tau value at each wavelength.",
            1.
        )
    }

    def _process(self, sed, parametres):
        """Add the redshift + IGM attenuation effect to the SED

        Parametres
        ----------
        sed  : pcigale.sed.SED object
        parametres : dictionnary
           Dictionnary with the module parametres (redshift and rtau)

        """

        if parametres['redshift'] == 0:
            # If redshift is 0, we do nothing
            pass
        else:
            # We need a lsst.Sed object to use its methods
            lsstSed = lsst.Sed()

            # First, we get the redshifted spectrum of the galaxy
            wavelen, flambda = lsstSed.redshiftSED(parametres['redshift'],
                                                   parametres['dimming'],
                                                   sed.wavelength_grid,
                                                   sed.luminosity)

            wavelen, red_l_lambda = lsstSed.addIGMattenuation(
                parametres['redshift'],
                parametres['rtau'],
                wavelen=wavelen,
                flambda=flambda
            )
            # We only want to add the redshift + IGM 'effect' to the SED
            # object (even if this has no physical meaning). As the
            # redshifting change the wavelength grid, we must regrid both
            # before substracting.
            new_wavelen = utils.best_grid(sed.wavelength_grid, wavelen)

            init_l_lambda = np.interp(new_wavelen, sed.wavelength_grid,
                                      sed.luminosity)
            red_l_lambda = np.interp(new_wavelen, wavelen, red_l_lambda)

            igm_effect = red_l_lambda - init_l_lambda

            # Base name for adding information to the SED.
            name = self.name or 'igmattenuation'

            sed.add_module(name, parametres)

            sed.add_info(name + '_redshift', parametres['redshift'])
            sed.add_info(name + '_rtau', parametres['rtau'])

            sed.add_contribution(
                name,
                new_wavelen,
                igm_effect
            )
