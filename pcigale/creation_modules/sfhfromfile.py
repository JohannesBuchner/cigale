# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

"""
Read star formation history from file module
============================================

This module reads the star formation history in a file.

"""

from collections import OrderedDict

import numpy as np

from ..utils import read_table
from . import CreationModule


class SfhFromFile(CreationModule):
    """Module reading the SFH from a file

    This module is used to read the Star Formation Histories from a FITS or
    VO-Table file. The first column must contain the time values (in Myr) and
    each other column may contain the Star Formation Rates (in solar mass per
    year) corresponding. Each SFR may be cut and normalised to 1 solar mass
    produced at the desired age.

    """

    parameter_list = OrderedDict([
        ("filename", (
            "str",
            "Name of the file containing the SFH. The first column must be "
            "the time [Myr] and the other column must contain the SFR "
            "[Msun/yr].",
            None
        )),
        ("sfr_column", (
            "integer",
            "List of column numbers where the star formation rates will "
            "be read.",
            None
        )),
        ("age", (
            "integer",
            "Age [Myr] where each SFH will be looked at.",
            None
        )),
        ("normalise", (
            "boolean",
            "Normalise the SFH to one solar mass produced at the given age.",
            "True"
        ))
    ])

    def process(self, sed):
        """Add the SFH read from the file.

        Parameters
        ----------
        sed: pcigale.sed.SED object
        parameters: dictionary containing the parameters

        """
        filename = self.parameters['filename']
        table = read_table(filename)

        time_grid = table.columns[0].data

        # -1 because Python indexes start to 0.
        sfr_column_number = int(self.parameters['sfr_column']) - 1
        sfr = table.columns[sfr_column_number].data

        age = int(self.parameters['age'])
        normalise = (self.parameters["normalise"].lower() == "true")

        # We cut the SFH to the desired age.
        sfr = sfr[time_grid <= age]
        time_grid = time_grid[time_grid <= age]

        # Compute the galaxy mass and normalise the SFH to 1 solar mass
        # produced if asked to.
        sfr_integrated = np.sum(sfr) * 1e6
        if normalise:
            sfr /= sfr_integrated
            sfr_integrated = 1.

        sed.add_module(self.name, self.parameters)
        sed.sfh = (time_grid, sfr)
        sed.add_info("sfh.integrated", sfr_integrated, True)
        sed.add_info("sfh.id", sfr_column_number+1)

# CreationModule to be returned by get_module
Module = SfhFromFile
