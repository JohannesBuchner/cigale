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
from astropy.table import Table, Column

# Directory where the output files are stored
OUT_DIR = "out/"


class ParametersHandler(object):
    """Class to handle the parameters to generate a parameter list on-the-fly.
    """

    def __init__(self, modules, params):
        """Instantiate the class.

        Parameters
        ----------
        modules: list
            Contains the modules in the order they are called
        params: OrderedDict
            Contains a list of parameters for each module

        """
        self.modules = modules
        self.parameters = [self._param_dict_combine(dictionary)
                           for dictionary in params]
        self.shape = tuple(len(parameter) for parameter in self.parameters)
        self.size = int(np.product(self.shape))

    def _param_dict_combine(self, dictionary):
        """Given a dictionary associating to each key an array, returns all the
        possible dictionaries associating a single element to each key.

        Parameters
        ----------
        dictionary: dict
            Dictionary associating an array to its (or some of its) keys.

        Returns
        -------
        combination_list: list of dictionaries
            List of dictionaries with the same keys but associating one element
            to each.

        """
        # We make a copy of the dictionary as we are modifying it.
        dictionary = collections.OrderedDict(dictionary)

        # First, we must ensure that all values are lists; when a value is a
        # single element, we put it in a list.
        # We must take a special care of strings, because they are iterable.

        for key, value in dictionary.items():
            if ((not isinstance(value, collections.Iterable)) or
                    isinstance(value, str)):
                dictionary[key] = [value]

        # We use itertools.product to make all the possible combinations from the
        # value lists.
        key_list = dictionary.keys()
        value_array_list = [dictionary[key] for key in key_list]
        combination_list = [collections.OrderedDict(zip(key_list, combination))
                            for combination in
                            itertools.product(*value_array_list)]

        return combination_list

    def from_index(self, index):
        """Provides the parameters of a model given a 1D index.

        Parameters
        ----------
        index: int
            1D index of the model for which we want the parameters

        Returns
        -------
        params: list
            Parameters of the model corresponding to the index

        """
        # Our problem is isomorph to the conversion between a linear and an nD
        # index of an nD array. Thankfully numpy's unravel_index does the
        # conversion from a 1D index to nD indices.
        indices = np.unravel_index(index, self.shape)
        params = [self.parameters[module][param_idx]
                  for module, param_idx in enumerate(indices)]

        return params

    def index_module_changed(self, index1, index2):
        """Find the index of the first module affected by a change of parameters.

        Parameters
        ----------
        index1: int
            First index
        index2: int
            Second index

        Returns
        -------
        module_idx: int
            Index of the first module that has a different parameter

        """
        indices = np.unravel_index((index1, index2), self.shape)
        for module_idx, (i, j) in enumerate(indices):
            if i != j:
                return module_idx

        return len(self.shape)


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


def save_fluxes(model_fluxes, model_parameters, filters, names, filename,
                directory=OUT_DIR, out_format='ascii'):
    """Save fluxes and associated parameters into a table.

    Parameters
    ----------
    model_fluxes: RawArray
        Contains the fluxes of each model.
    model_parameters: RawArray
        Contains the parameters associated to each model.
    filters: OrderedDict
        Contains the filters.
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
                      names=list(filters.keys()) + list(names))

    out_table.write("{}/{}".format(directory, filename), format=out_format)
