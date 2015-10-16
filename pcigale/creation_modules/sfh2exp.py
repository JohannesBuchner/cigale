# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

"""
Double decreasing exponential star formation history module
===========================================================

This module implements a star formation history (SFH) composed of two
decreasing exponentials.

"""

from collections import OrderedDict

import numpy as np

from . import CreationModule


class Sfh2Exp(CreationModule):
    """Double decreasing exponential Star Formation History

    This module sets the SED star formation history (SFH) as a combination of
    two exp(-t/τ) exponentials.

    """

    parameter_list = OrderedDict([
        ("tau_main", (
            "float",
            "e-folding time of the main stellar population model in Myr.",
            6000.
        )),
        ("tau_burst", (
            "float",
            "e-folding time of the late starburst population model in Myr.",
            50.
        )),
        ("f_burst", (
            "float",
            "Mass fraction of the late burst population.",
            0.01
        )),
        ("age", (
            "integer",
            "Age of the main stellar population in the galaxy in Myr."
            "The precision is 1 Myr.",
            5000.
        )),
        ("burst_age", (
            "integer",
            "Age of the late burst in Myr. Precision is 1 Myr.",
            20.
        )),
        ("sfr_0", (
            "float",
            "Value of SFR at t = 0 in M_sun/yr.",
            1.
        )),
        ("normalise", (
            "boolean",
            "Normalise the SFH to produce one solar mass.",
            "True"
        )),
    ])

    def process(self, sed):
        """Add a double decreasing exponential Star Formation History.

        Parameters
        ----------
        sed: pcigale.sed.SED object

        """
        tau_main = float(self.parameters["tau_main"])
        tau_burst = float(self.parameters["tau_burst"])
        f_burst = float(self.parameters["f_burst"])
        age = int(self.parameters["age"])
        burst_age = int(self.parameters["burst_age"])
        sfr_0 = int(self.parameters["sfr_0"])
        normalise = (self.parameters["normalise"].lower() == "true")

        # Time grid and age. If needed, the age is rounded to the inferior Myr
        time_grid = np.arange(1, age + 1)
        time_grid_burst = np.arange(1, burst_age + 1)

        # SFR for each component
        sfr = np.exp(-time_grid / tau_main)
        sfr_burst = np.exp(-time_grid_burst / tau_burst)

        # Height of the late burst to have the desired produced mass fraction
        sfr_burst *= f_burst / (1.-f_burst) * np.sum(sfr) / np.sum(sfr_burst)

        # We add the age burst exponential for ages superior to age -
        # burst_age
        sfr[-time_grid_burst[-1]:] += sfr_burst

        # Compute the galaxy mass and normalise the SFH to 1 solar mass
        # produced if asked to.
        galaxy_mass = np.sum(sfr) * 1e6
        if normalise:
            sfr /= galaxy_mass
            galaxy_mass = 1.
        else:
            sfr *= sfr_0
            galaxy_mass *= sfr_0

        sed.add_module(self.name, self.parameters)

        # Add the sfh and the output parameters to the SED.
        sed.sfh = (time_grid, sfr)
        sed.add_info("galaxy_mass", galaxy_mass, True)
        sed.add_info("sfh.tau_main", tau_main)
        sed.add_info("sfh.tau_burst", tau_burst)
        sed.add_info("sfh.f_burst", f_burst)
        sed.add_info("sfh.burst_age", burst_age)

# CreationModule to be returned by get_module
Module = Sfh2Exp
