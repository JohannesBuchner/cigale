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
            "the time in Myr, starting from 0 with a step of 1 Myr. The other "
            "columns must contain the SFR in Msun/yr."
            "[Msun/yr].",
            None
        )),
        ("sfr_column", (
            "integer",
            "List of column indices of the SFR. The first SFR column has the "
            "index 1.",
            None
        )),
        ("age", (
            "integer",
            "Age in Myr at which the SFH will be looked at.",
            None
        )),
        ("normalise", (
            "boolean",
            "Normalise the SFH to one solar mass produced at the given age.",
            "True"
        ))
    ])

    def _init_code(self):
        filename = self.parameters['filename']
        normalise = (self.parameters["normalise"].lower() == "true")
        age = int(self.parameters['age'])
        self.sfr_column_number = int(self.parameters['sfr_column'])

        table = read_table(filename)
        self.sfr = table.columns[self.sfr_column_number].data.astype(np.float)
        self.time_grid = table.columns[0].data.astype(np.int)
        if self.time_grid[0] != 0:
            raise Exception("The time grid must start from 0.")
        if np.all(self.time_grid[1:]-self.time_grid[:-1] == 1) == False:
            raise Exception("The time step must be 1 Myr. Computed models will "
                            "be wrong.")

        # We cut the SFH to the desired age.
        self.sfr = self.sfr[self.time_grid <= age]
        self.time_grid = self.time_grid[self.time_grid <= age]

        # Compute the galaxy mass and normalise the SFH to 1 solar mass
        # produced if asked to.
        self.sfr_integrated = np.sum(self.sfr) * 1e6
        if normalise:
            self.sfr /= self.sfr_integrated
            self.sfr_integrated = 1.

    def process(self, sed):
        """Add the SFH read from the file.

        Parameters
        ----------
        sed: pcigale.sed.SED object
        parameters: dictionary containing the parameters

        """

        sed.add_module(self.name, self.parameters)
        sed.sfh = (self.time_grid, self.sfr)
        sed.add_info("sfh.integrated", self.sfr_integrated, True)
        sed.add_info("sfh.index", self.sfr_column_number)

# CreationModule to be returned by get_module
Module = SfhFromFile
