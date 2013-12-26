# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly <yannick.roehlly@oamp.fr>

import numpy as np
from collections import OrderedDict
from . import common

# Time lapse used in the age grid in Myr. If should be consistent with the
# time lapse in the SSP modules.
AGE_LAPSE = 1


class Module(common.SEDCreationModule):
    """Create a double decreasing exponential Star Formation History

    This module create a star formation history (star formation rate vs galaxy
    age) composed of two exp(-t/τ). The SFH is added to the info dictionary of
    the SED as a tuple (age, SFR).

    """

    parameter_list = OrderedDict([
        ("tau_main", (
            "float",
            "e-folding time of the main stellar population model in Myr.",
            None
        )),
        ("tau_burst", (
            "float",
            "e-folding time of the late starburst population model in Myr.",
            None
        )),
        ("f_burst", (
            "float",
            "Mass fraction of the late burst population.",
            None
        )),
        ("age", (
            "integer",
            "Age of the oldest stars in the galaxy in Myr. The precision "
            "is 1 Myr.",
            None
        )),
        ("burst_age", (
            "integer",
            "Age of the late burst in Myr. Precision is 1 Myr.",
            None
        ))
    ])

    out_parameter_list = OrderedDict([
        ("tau_main", "e-folding time of the main stellar population model "
                     "in Myr."),
        ("tau_burst", "e-folding time of the late starburst population model "
                      "in Myr."),
        ("f_burst", "Produced mass fraction of the late burst population."),
        ("age", "Age of the oldest stars in the galaxy in Myr."),
        ("burst_age", "Age of the late burst in Myr.")
    ])

    def process(self, sed):
        """Add a double decreasing exponential Star Formation History.

        Parameters
        ----------
        sed : pcigale.sed.SED object

        """

        tau_main = float(self.parameters["tau_main"])
        tau_burst = float(self.parameters["tau_burst"])
        f_burst = float(self.parameters["f_burst"])
        age = int(self.parameters["age"])
        burst_age = int(self.parameters["burst_age"])

        # Time grid and age. If needed, the age is rounded to the inferior Myr
        time_grid = np.arange(AGE_LAPSE, age + AGE_LAPSE, AGE_LAPSE)
        age = np.max(time_grid)

        # Main exponential
        sfr = np.exp(-time_grid / tau_main)

        # Height of the late burst to have the desired produced mass fraction
        # (assuming that the main burst as a height of 1).
        burst_height = (
            f_burst * tau_main * (1 - np.exp(-age / tau_main))
            / (
                (1 - f_burst) * tau_burst * np.exp(-burst_age / tau_burst)
            ))

        # We add the age burst exponential for ages superior to age -
        # burst_age
        mask = (time_grid >= (age - burst_age))
        sfr[mask] = sfr[mask] + burst_height * np.exp(
            (-time_grid[mask] + age - burst_age) / tau_burst)

        # We normalise the SFH to have one solar mass produced.
        sfr = sfr / np.trapz(sfr * 1.e6, time_grid)

        sed.add_module(self.name, self.parameters)

        # Add the sfh and the output parameters to the SED.
        sed.sfh = (time_grid, sfr)
        sed.add_info("sfh_tau_main" + self.postfix, tau_main)
        sed.add_info("sfh_tau_burst" + self.postfix, tau_burst)
        sed.add_info("sfh_f_burst" + self.postfix, f_burst)
        sed.add_info("sfh_burst_age" + self.postfix, burst_age)
