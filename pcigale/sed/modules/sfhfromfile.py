# -*- coding: utf-8 -*-
"""
Copyright (C) 2013 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""

import atpy
import numpy as np
from . import common


class Module(common.SEDCreationModule):
    """Module reading the SFH from a file

    This module is used to read the Star Formation Histories from a FITS or
    VO-Table file. The first column must contain the time values (in Myr) and
    each other column may contain the Star Formation Rates (in solar mass per
    year) corresponding. Each SFR may be cut and normalised to 1 solar mass
    produced at the desired age.

    """

    parameter_list = {
        "filename": (
            "str",
            "Name of the file containing the SFH. The first column must be "
            "the time [Myr] and the other column must contain the SFR "
            "[Msun/yr].",
            None
        ),
        "sfr_column": (
            "int",
            "List of column numbers where the star formation rates will "
            "be read..",
            None
        ),
        "age": (
            "float",
            "List of ages [Myr] where each SFH will be looked at.",
            None
        )
    }

    def _process(self, sed, parameters):
        """Add the SFH read from the file.

        Parameters
        ----------
        sed  : pcigale.sed.SED object
        parameters : dictionary containing the parameters

        """
        filename = parameters['filename']
        table = atpy.Table(filename, verbose=False)

        time_column_name = table.columns.keys[0]
        time_grid = table[time_column_name]

        # -1 because Python indexes start to 0.
        sfr_column_number = int(parameters['sfr_column']) - 1
        sfr_column_name = table.columns.keys[sfr_column_number]
        sfr = table[sfr_column_name]

        age = parameters['age']

        # We cut the SFH to the desired age.
        sfr = sfr[time_grid <= age]
        time_grid = time_grid[time_grid <= age]

        # The we normalise it to 1 solar mass produced.
        sfr = sfr / np.trapz(sfr * 1.e6, time_grid)

        # Base name for adding information to the SED.
        name = self.name or 'loadfile'

        sed.add_module(name, parameters)
        sed.add_info(name + "_sfh", (time_grid, sfr))
        sed.add_info(name + "_age", age)
        sed.add_info(name + "_sfh_id", sfr_column_name)
