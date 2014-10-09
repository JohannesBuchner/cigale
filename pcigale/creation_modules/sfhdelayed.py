# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2014 Laboratoire d'Astrophysique de Marseille
# Copyright (C) 2014 University of Cambridge
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly, Véronique Buat & Médéric Boquien

"""
Delayed tau model for star formation history
============================================

This module implements a star formation history (SFH) described as a delayed
rise of the SFR up to a maximum, followed by an exponential decrease.

"""

import numpy as np
from collections import OrderedDict
from . import CreationModule

# Time lapse used in the age grid in Myr. If should be consistent with the
# time lapse in the SSP modules.
AGE_LAPSE = 1


class SFHDelayed(CreationModule):
    """Delayed tau model for Star Formation History

    This module sets the SED star formation history (SFH) proportional to time,
    with a declining exponential parametrised with a time-scale (tau_main).

    """

    parameter_list = OrderedDict([
        ("tau_main", (
            "float",
            "e-folding time of the main stellar population model in Myr.",
            None
        )),
        ("age", (
            "integer",
            "Age of the oldest stars in the galaxy in Myr. The precision "
            "is 1 Myr.",
            None
        )),
         ("sfr_A", (
            "float",
            "Multiplicative factor controlling the amplitude of SFR.",
            1.
        )),
        ("normalise", (
            "boolean",
            "Normalise the SFH to produce one solar mass.",
            "False"
        ))
   ])

    out_parameter_list = OrderedDict([
        ("tau_main", "e-folding time of the main stellar population model "
                     "in Myr."),
        ("age", "Age of the oldest stars in the galaxy in Myr."),
        ("sfr_A", "Multiplicative factor controlling the amplitude of SFR.")
    ])

    def process(self, sed):
        """
        Parameters
        ----------
        sed : pcigale.sed.SED object

        """
        tau_main = float(self.parameters["tau_main"])
        age = int(self.parameters["age"])
        sfr_A = int(self.parameters["sfr_A"])
        normalise = (self.parameters["normalise"].lower() == "true")

        # Time grid and age. If needed, the age is rounded to the inferior Myr
        time_grid = np.arange(AGE_LAPSE, age + AGE_LAPSE, AGE_LAPSE)

        # Main SFR
        sfr = sfr_A * time_grid / tau_main**2 * np.exp(-time_grid / tau_main)

        # Normalise the SFH to 1 solar mass produced if asked to.
        if normalise:
            sfr = sfr / np.trapz(sfr * 1e6, time_grid)

        sed.add_module(self.name, self.parameters)

        # Add the sfh and the output parameters to the SED.
        sed.sfh = (time_grid, sfr)
        sed.add_info("sfh.tau_main", tau_main)

# CreationModule to be returned by get_module
Module = SFHDelayed
