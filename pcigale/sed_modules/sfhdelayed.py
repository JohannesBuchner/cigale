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

from collections import OrderedDict

import numpy as np

from . import SedModule


class SFHDelayed(SedModule):
    """Delayed tau model for Star Formation History

    This module sets the SED star formation history (SFH) proportional to time,
    with a declining exponential parametrised with a time-scale (tau_main).

    """

    parameter_list = OrderedDict([
        ("tau_main", (
            "cigale_list()",
            "e-folding time of the main stellar population model in Myr.",
            2000.
        )),
        ("age", (
            "cigale_list(dtype=int, minvalue=0.)",
            "Age of the oldest stars in the galaxy in Myr. The precision "
            "is 1 Myr.",
            5000
        )),
        ("sfr_A", (
            "float(min=0.)",
            "Multiplicative factor controlling the amplitude of SFR.",
            1.
        )),
        ("normalise", (
            "boolean()",
            "Normalise the SFH to produce one solar mass.",
            True
        ))
    ])

    def _init_code(self):
        self.tau_main = float(self.parameters["tau_main"])
        age = int(self.parameters["age"])
        sfr_A = float(self.parameters["sfr_A"])
        normalise = bool(self.parameters["normalise"])

        # Time grid and SFR
        time_grid = np.arange(age)
        self.sfr = time_grid * np.exp(-time_grid / self.tau_main) / \
                   self.tau_main**2

        # Compute the galaxy mass and normalise the SFH to 1 solar mass
        # produced if asked to.
        self.sfr_integrated = np.sum(self.sfr) * 1e6
        if normalise:
            self.sfr /= self.sfr_integrated
            self.sfr_integrated = 1.
        else:
            self.sfr *= sfr_A
            self.sfr_integrated *= sfr_A

    def process(self, sed):
        """
        Parameters
        ----------
        sed : pcigale.sed.SED object

        """
        sed.add_module(self.name, self.parameters)

        # Add the sfh and the output parameters to the SED.
        sed.sfh = self.sfr
        sed.add_info("sfh.integrated", self.sfr_integrated, True)
        sed.add_info("sfh.tau_main", self.tau_main)

# SedModule to be returned by get_module
Module = SFHDelayed
