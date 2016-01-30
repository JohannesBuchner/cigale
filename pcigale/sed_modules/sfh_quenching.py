# -*- coding: utf-8 -*-

# Copyright (C) 2015 Laboratoire d'Astrophysique de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

"""
Star Formation History quenching module
=======================================

This module performs a quenching on the Star Formation History. Below a given
age, the Star Formation Rate in multiplied by 1 - quenching_factor and is set
constant.

"""

from collections import OrderedDict

import numpy as np

from . import SedModule


class SfhQuench(SedModule):
    """Star Formation History Quenching

    This module implements a quenching of the Star Formation History.

    """

    parameter_list = OrderedDict([
        ("quenching_age", (
            "cigale_list(dtype=int, minvalue=0.)",
            "Look-back time when the quenching happens in Myr.",
            0
        )),
        ("quenching_factor", (
            "cigale_list(minvalue=0., maxvalue=1.)",
            "Quenching factor applied to the SFH. After the quenching time, "
            "the SFR is multiplied by 1 - quenching factor and made constant. "
            "The factor must be between 0 (no quenching) and 1 (no more star "
            "formation).",
            0.
        )),
        ("normalise", (
            "boolean()",
            "Normalise the SFH to produce one solar mass.",
            True
        ))
    ])

    def _init_code(self):
        self.quenching_age = int(self.parameters["quenching_age"])
        self.quenching_factor = float(self.parameters["quenching_factor"])
        self.normalise = bool(self.parameters["normalise"])

    def process(self, sed):
        """
        Parameters
        ----------
        sed : pcigale.sed.SED object

        """
        # Read the star formation history of the SED
        time, sfr = sed.sfh

        if self.quenching_age > time[-1]:
            raise Exception("[sfh_quenching] The quenching age is greater "
                            "than the galaxy age. Please fix your parameters.")

        # We assume the time in the star formation history is evenly spaced to
        # compute the reverse index (i.e. from the end of the array) of the SFH
        # step corresponding to the quenching age.
        # We make the computation only if the quenching age and the quenching
        # factor are not 0.
        if self.quenching_age > 0 and self.quenching_factor > 0.:
            age_lapse = time[1] - time[0]
            quenching_idx = np.int(self.quenching_age/age_lapse)
            sfr[-quenching_idx:] = sfr[-quenching_idx] * (
                1 - self.quenching_factor)

            # Compute the galaxy mass and normalise the SFH to 1 solar mass
            # produced if asked to.
            sfr_integrated = np.sum(sfr) * 1e6
            if self.normalise:
                sfr /= sfr_integrated
                sfr_integrated = 1.

            sed.sfh = (time, sfr)
            sed.add_info("sfh.integrated", sfr_integrated, True, force=True)

        sed.add_module(self.name, self.parameters)

        sed.add_info("sfh.quenching_age", self.quenching_age)
        sed.add_info("sfh.quenching_factor", self.quenching_factor)

# SedModule to be returned by get_module
Module = SfhQuench
