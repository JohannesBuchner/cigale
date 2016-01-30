# -*- coding: utf-8 -*-
# Copyright (C) 2014 Médéric Boquien
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Médéric Boquien

"""
Various utility functions for pcigale analysis modules
"""

from datetime import datetime
import collections
import itertools
import os
import shutil

import numpy as np
from astropy import log
from astropy.table import Table, Column

log.setLevel('ERROR')

# Directory where the output files are stored
OUT_DIR = "out/"


def backup_dir(directory=OUT_DIR):
    if os.path.exists(directory):
        new_name = datetime.now().strftime("%Y%m%d%H%M") + "_" + directory
        os.rename(directory, new_name)
        print("The existing {} directory was renamed to {}".format(
            directory,
            new_name
        ))
    os.mkdir(directory)
    shutil.copyfile('pcigale.ini', directory + 'pcigale.ini')
    shutil.copyfile('pcigale.ini.spec', directory + 'pcigale.ini.spec')


def save_fluxes(model_fluxes, model_parameters, filters, names, filename,
                directory=OUT_DIR, out_format='ascii.commented_header'):
    """Save fluxes and associated parameters into a table.

    Parameters
    ----------
    model_fluxes: RawArray
        Contains the fluxes of each model.
    model_parameters: RawArray
        Contains the parameters associated to each model.
    filters: list
        Contains the filter names.
    names: List
        Contains the parameters names.
    filename: str
        Name under which the file should be saved.
    directory: str
        Directory under which the file should be saved.
    out_format: str
        Format of the output file

    """
    out_fluxes = np.ctypeslib.as_array(model_fluxes[0])
    out_fluxes = out_fluxes.reshape(model_fluxes[1])

    out_params = np.ctypeslib.as_array(model_parameters[0])
    out_params = out_params.reshape(model_parameters[1])

    out_table = Table(np.hstack((out_fluxes, out_params)),
                      names=filters + list(names))

    out_table.add_column(Column(np.arange(model_fluxes[1][0]), name='id'),
                         index=0)

    out_table.write("{}/{}".format(directory, filename), format=out_format)
