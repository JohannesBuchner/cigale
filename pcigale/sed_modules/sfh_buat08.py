# -*- coding: utf-8 -*-

# Copyright (C) 2015 Laboratoire d'Astrophysique de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Alessandro Boselli, Yannick Roehlly

"""
Physically motivated star formation history
===========================================

This module implements the star formation history (SFH) described in Buat et
al. (2008). It's an analytical star formation history resulting from the non
linear combinations of a star formation law, an infall history and its mass
dependence.

In Buat et al. (2005), the authors have fitted the polynomial formulae
log10(SFR(t)) = a + b log10(t) + c t^0.5 to their chemical evolution generated
SFH and give the values of the a, b and c parameters for different values of
the rotational velocity of the galaxy. We use this velocity as input parameter
and interpolate the values of a, b and c.

"""

from collections import OrderedDict

import numpy as np

from . import SedModule


class SfhBuat08(SedModule):
    """Chemical evolution motivated Star Formation History

    This module implements a chemical evolution motivated star formation
    history. The rotational velocity, meaningful for nearby galaxy, is used as
    input parameters.

    """

    parameter_list = OrderedDict([
        ("velocity", (
            "cigale_list(minvalue=80., maxvalue=360.)",
            "Rotational velocity of the galaxy in km/s. Must be between 80 "
            "and 360 (included).",
            200.
        )),
        ("age", (
            "cigale_list(dtype=int, minvalue=0.)",
            "Age of the oldest stars in the galaxy. The precision "
            "is 1 Myr.",
            5000
        )),
        ("normalise", (
            "boolean()",
            "Normalise the SFH to produce one solar mass.",
            True
        ))
    ])

    def _init_code(self):
        self.velocity = float(self.parameters["velocity"])
        self.age = int(self.parameters["age"])
        normalise = bool(self.parameters["normalise"])

        # Time grid and age. If needed, the age is rounded to the inferior Myr
        time_grid = np.arange(self.age)

        # Values from Buat et al. (2008) table 2
        paper_velocities = np.array([80., 150., 220., 290., 360.])
        paper_as = np.array([6.62, 8.74, 10.01, 10.81, 11.35])
        paper_bs = np.array([0.41, 0.98, 1.25, 1.35, 1.37])
        paper_cs = np.array([0.36, -0.20, -0.55, -0.74, -0.85])

        # Interpolation of a, b, c corresponding to the velocity.
        a = np.interp(self.velocity, paper_velocities, paper_as)
        b = np.interp(self.velocity, paper_velocities, paper_bs)
        c = np.interp(self.velocity, paper_velocities, paper_cs)

        # Main SFR
        t = (time_grid+1) / 1000  # The time is in Gyr in the formulae
        self.sfr = 10.**(a + b * np.log10(t) + c * t**.5) / 1.e9

        # Compute the galaxy mass and normalise the SFH to 1 solar mass
        # produced if asked to.
        self.sfr_integrated = np.sum(self.sfr) * 1e6
        if normalise:
            self.sfr /= self.sfr_integrated
            self.sfr_integrated = 1.

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
        sed.add_info("sfh.velocity", self.velocity)

# SedModule to be returned by get_module
Module = SfhBuat08
