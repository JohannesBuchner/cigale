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

from . import CreationModule


# Time lapse used in the age grid in Myr. If should be consistent with the
# time lapse in the SSP modules.
AGE_LAPSE = 1


class SfhQuench(CreationModule):
    """Star Formation History Quenching

    This module implements a quenching of the Star Formation History.

    """

    parameter_list = OrderedDict([
        ("quenching_age", (
            "integer",
            "Age of the galaxy at which the quenching happens in Myr.",
            0
        )),
        ("quenching_factor", (
            "float",
            "Quenching factor applied to the SFH. After the quenching time, "
            "the SFR is multiplied by 1 - quenching factor and made constant. "
            "The factor must be between 0 (no quenching) and 1 (no more star "
            "formation).",
            0
        )),
        ("normalise", (
            "boolean",
            "Normalise the SFH to produce one solar mass.",
            "True"
        ))
    ])

    def process(self, sed):
        """
        Parameters
        ----------
        sed : pcigale.sed.SED object

        """
        quenching_age = int(self.parameters["quenching_age"])
        quenching_factor = float(self.parameters["quenching_factor"])
        normalise = (self.parameters["normalise"].lower() == "true")

        # Read the star formation history of the SED
        time, sfr = sed.sfh

        # We assume the time in the star formation history is evenly spaced to
        # compute the reverse index (i.e. from the end of the array) of the SFH
        # step corresponding to the quenching age.
        # We make the computation only if the quenching age and the quenching
        # factor are not 0.
        if quenching_age and quenching_factor:
            age_lapse = time[1] - time[0]
            quenching_idx = np.int(quenching_age/age_lapse)
            sfr[-quenching_idx:] = sfr[-quenching_idx] * (
                1 - quenching_factor)

            # Compute the galaxy mass and normalise the SFH to 1 solar mass
            # produced if asked to.
            galaxy_mass = np.trapz(sfr, time) * 1e6
            if normalise:
                sfr /= galaxy_mass
                galaxy_mass = 1.

            sed.sfh = (time, sfr)
            sed.add_info("galaxy_mass", galaxy_mass, True, force=True)

        sed.add_module(self.name, self.parameters)

        sed.add_info("sfh.quenching_age", quenching_age)
        sed.add_info("sfh.quenching_factor", quenching_factor)

# CreationModule to be returned by get_module
Module = SfhQuench
