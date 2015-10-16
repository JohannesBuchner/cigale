# Copyright (C) 2015 Centre de données Astrophysiques de Marseille
# Copyright (C) 2015 Institute of Astronomy
# Copyright (C) 2015 Universidad de Antofagasta
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Denis Burgarella & Médéric Boquien

"""
Periodic SFH in the form of rectangles, and decaying or delayed exponentials
============================================================================

# This module implements a periodic star formation history (SFH) formed by
regularly-spaced star formation events. Each even can either be rectangular, a
decaying exponential, or "delayed".

"""

from collections import OrderedDict

import numpy as np

from . import CreationModule


class SfhPeriodic(CreationModule):
    """Several regularly-spaced short delayed-SFH SF events

    This module sets the SED star formation history (SFH) as a combination of
    several regularly-spaced short SF events.

    """

    parameter_list = OrderedDict([
        ("type_bursts", (
            "integer",
            "Type of the individual star formation episodes. 0: exponential, "
            "1: delayed, 2: rectangle.",
            0
        )),
        ("delta_bursts", (
            "integer",
            "Elapsed time between the beginning of each burst in Myr. The "
            "precision is 1 Myr.",
            50
        )),
        ("tau_bursts", (
            "integer",
            "Duration (rectangle) or e-folding time of all short events in "
            "Myr. The precision is 1 Myr.",
            20
        )),
        ("age", (
            "integer",
            "Age of the main stellar population in the galaxy in Myr. The "
            "precision is 1 Myr.",
            1000
        )),
        ("sfr_A", (
            "float",
            "Multiplicative factor controlling the amplitude of SFR (valid "
            "for each event).",
            1.
        )),
        ("normalise", (
            "boolean",
            "Normalise the SFH to produce one solar mass.",
            "True"
        )),
    ])

    def _init_code(self):
        self.type_bursts = int(self.parameters["type_bursts"])
        self.delta_bursts = int(self.parameters["delta_bursts"])
        self.tau_bursts = int(self.parameters["tau_bursts"])
        age = int(self.parameters["age"])
        sfr_A = float(self.parameters["sfr_A"])
        normalise = (self.parameters["normalise"].lower() == "true")

        self.time_grid = np.arange(0, age)
        self.sfr = np.zeros_like(self.time_grid, dtype=np.float)

        if self.type_bursts == 0:
            burst = np.exp(-self.time_grid/self.tau_bursts)
        elif self.type_bursts == 1:
            burst = np.exp(-self.time_grid/self.tau_bursts) * \
                    self.time_grid/self.tau_bursts**2
        elif self.type_bursts == 2:
            burst = np.zeros_like(self.time_grid)
            burst[:self.tau_bursts+1] = 1.
        else:
            raise Exception("Burst type {} unknown.".format(self.type_bursts))

        for t_burst in np.arange(0, age, self.delta_bursts):
            self.sfr += burst
            burst = np.roll(burst, self.delta_bursts)
            burst[:self.delta_bursts] = 0.

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
        """Add a star formation history formed by several regularly-spaced SF
        events.

        ** Parameters **

        sed: pcigale.sed.SED object

        """

        sed.add_module(self.name, self.parameters)

        sed.sfh = (self.time_grid, self.sfr)
        sed.add_info("sfh.integrated", self.sfr_integrated, True)
        sed.add_info("sfh.type_bursts", self.type_bursts)
        sed.add_info("sfh.delta_bursts", self.delta_bursts)
        sed.add_info("sfh.tau_bursts", self.tau_bursts)

# CreationModule to be returned by get_module
Module = SfhPeriodic
