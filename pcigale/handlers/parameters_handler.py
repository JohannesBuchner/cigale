# -*- coding: utf-8 -*-
# Copyright (C) 2014 University of Cambridge
# Copyright (C) 2016 Universidad de Antofagasta
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Médéric Boquien

import collections
import itertools
import numpy as np


class ParametersHandler(object):
    """Class to abstract the call to the relevant parameters handler depending
    how the physical parameters of the models are provided (directly in the
    pcigale.ini file ).

    A ParametersHandler allows to generate a list containing the parameters for
    a modules whose index is passed as an argument. It also allows to know what
    modules have changed their parameters given their indices. Because the
    order of the modules is optimised to minimise the computations, this allows
    to keep the cache of partially computed models in SedWarehouse as small as
    possible by weeding out partial models that will not be used anymore."""

    def __new__(object, configuration):
            return ParametersHandlerGrid(configuration)


class ParametersHandlerGrid(object):
    """Class to generate a parameters handler for a systematic grid using the
    parameters given in the pcigale.ini file."""

    def __init__(self, configuration):
        """Instantiate the class.

        Parameters
        ----------
        configuration: dictionary
            Contains the modules in the order they are called

        """
        self.modules = configuration['creation_modules']
        self.parameters = [self._param_dict_combine(dictionary) for dictionary
                           in configuration['creation_modules_params']]
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
        dictionary = dict(dictionary)

        # First, we must ensure that all values are lists; when a value is a
        # single element, we put it in a list.
        # We must take a special care of strings, because they are iterable.

        for key, value in dictionary.items():
            if ((not isinstance(value, collections.Iterable)) or
                    isinstance(value, str)):
                dictionary[key] = [value]

        # We use itertools.product to make all the possible combinations from
        # the value lists.
        key_list = dictionary.keys()
        value_array_list = [dictionary[key] for key in key_list]
        combination_list = [dict(zip(key_list, combination))
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
