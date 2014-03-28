# -*- coding: utf-8 -*-
# Copyright (C) 2012 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

import itertools
import collections


def param_dict_combine(dictionary):
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
    # We must take a special care of strings, because they are iterable. We
    # define our own string_type to work both with Python 2 and Python 3.
    try:
        string_type = basestring
    except NameError:
        string_type = str

    for key, value in dictionary.items():
        if ((not isinstance(value, collections.Iterable)) or
                isinstance(value, string_type)):
            dictionary[key] = [value]

    # We use itertools.product to make all the possible combinations from the
    # value lists.
    key_list = dictionary.keys()
    value_array_list = [dictionary[key] for key in key_list]
    combination_list = [collections.OrderedDict(zip(key_list, combination))
                        for combination in
                        itertools.product(*value_array_list)]

    return combination_list
