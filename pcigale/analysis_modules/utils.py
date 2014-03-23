# -*- coding: utf-8 -*-
# Copyright (C) 2014 Médéric Boquien
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Médéric Boquien

"""
Various utility functions for pcigale analysis modules
"""

import numpy as np


def find_changed_parameters(list_parameters):
    """
    Given a list of parameters dictionaries corresponding each to a given SED,
    find which parameters have changed between two adjacent items in the list.

    When used for to find which SEDs to discard for partial cache clearing,
    this relies on the assumption that the list is properly ordered. In this
    case, when a parameter changes, the corresponding SEDs can be discarded.
    This should work when the list of dictionaries is generating using
    itertools.product().

    Parameters
    ----------
    list_parameters: list of list of dictionaries
        Each item is a list of dictionaries containing the parameters for each
        module.

    Return
    ------
    An array of integers containing the index of the module that has to be
    discarded. When several parameters have changed we return the lowest index.
    The cache cleaning routine can then just discard all SED with at least as
    many modules. This will keep the cache small if used consistently.

    """
    changed = [None] * len(list_parameters)
    for i in range(len(list_parameters)-1):
        for idx, (par, par_next) in enumerate(zip(list_parameters[i],
                                                  list_parameters[i+1])):
            for k in par.keys():
                if par[k] != par_next[k]:
                    if changed[i] is not None:
                        print('Warning! It went wrong in the cache cleaning')
                    changed[i] = idx
                    break
            if changed[i] is not None:
                break
    changed[-1] = 0

    return np.array(changed, dtype=np.int)
