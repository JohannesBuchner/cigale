# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Denis Burgarella

"""
SFH resembling a Dirac Comb with regularly-spaced delayed-SFH  events
=====================================================================

This module implements a star formation history (SFH) formed by several 
regularly-spaced short regularly-spaced delayed-SFH SF events.

"""

import numpy as np
from collections import OrderedDict
from . import CreationModule

# Time lapse used in the age grid in Myr. If should be consistent with the
# time lapse in the SSP modules.
AGE_LAPSE = 1


class SfhComb(CreationModule):
    """Several regularly-spaced short delayed-SFH SF events

    This module sets the SED star formation history (SFH) as a combination of
    several regularly-spaced short SF events.

    """

    parameter_list = OrderedDict([
        ("N_events", (
            "integer",
            "Number of individual star formation events. ",
            5
        )),
        ("tau_events", (
            "float",
            "e-folding time of all short events in Myr.",
            20.
        )),
        ("age", (
            "integer",
            "Age of the main stellar population in the galaxy in Myr."
            "The precision is 1 Myr.",
            1000
        )),
        ("age_last", (
            "integer",
            "Time since the end of the last SF event [Myr]. "
            "Warning: Depending on the parameters, the last burst might not be finished. "
            "The precision is 1 Myr.",
            100
        )),
        ("sfr_A", (
            "float",
            "Multiplicative factor controlling the amplitude of SFR (valid for each event).",
            1.
        )),
        ("normalise", (
            "boolean",
            "Normalise the SFH to produce one solar mass.",
            "True"
        )),
    ])

    out_parameter_list = OrderedDict([
        ("N_events", "Number of individual star formation events"),
        ("tau_events", "e-folding time of all short events in Myr.."),
        ("age", "Age of the stellar population in the galaxy in Myr."),
        ("age_last", "Time since the end of the last SF event in Myr."),
        ("sfr_A", "Height of each SF event in M_sun/yr.")
    ])

    def process(self, sed):
        """Add a star formation history formed by several 
        regularly-spaced delayed-SFH short SF events.

        ** Parameters **
        
        sed: pcigale.sed.SED object

        """
        N_events = int(self.parameters["N_events"])
        tau_events = float(self.parameters["tau_events"])
        age = int(self.parameters["age"])
        age_last = int(self.parameters["age_last"])
        sfr_A = int(self.parameters["sfr_A"])
        normalise = (self.parameters["normalise"].lower() == "true")

        # Time grid and age. If needed, the age is rounded to the inferior Myr
        time_grid = np.arange(AGE_LAPSE, age + AGE_LAPSE, AGE_LAPSE)
        age = np.max(time_grid)
        sfr = np.zeros(len(time_grid))

        # Build the Dirac comb of events.
        # Ages of the sequence of SF events in Myr, 
        # N_events over "age" with the last one starting at "age_last"
        time_event = np.linspace(0., age-age_last, num = N_events)

        # Add each delayed SFH
        for i_time in range(N_events):
           sfr[time_event[i_time]:] = \
              sfr_A * (time_grid[time_event[i_time]:]-time_event[i_time]) / tau_events**2 * \
              np.exp(-(time_grid[time_event[i_time]:]-time_event[i_time]) / tau_events)

        # Compute the galaxy mass and normalise the SFH to 1 solar mass
        # produced if asked to.
        galaxy_mass = np.trapz(sfr * 1e6, time_grid)
        if normalise:
            sfr = sfr / galaxy_mass
            galaxy_mass = 1.

        sed.add_module(self.name, self.parameters)

        # Add the sfh and the output parameters to the SED.
        sed.sfh = (time_grid, sfr)
        sed.add_info("galaxy_mass", galaxy_mass, True)
        sed.add_info("sfh.N_events", N_events)
        sed.add_info("sfh.tau_events", tau_events)
        sed.add_info("sfh.age_last", age_last)
        sed.add_info("sfh.age", age)

# CreationModule to be returned by get_module
Module = SfhComb
