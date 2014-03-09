# -*- coding: utf-8 -*-
# Copyright (C) 2014 Médéric Boquien
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Médéric Boquien

"""
Various utility functions for pcigale analysis modules
"""


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
    A list a tuples with the same size as the input list. Each tuple contains
    the parameter that has changed and its value. When several parameters have
    changed, it selects only the one that would discard the most models.
    """
    changed = [None] * len(list_parameters)
    for i in range(len(list_parameters)-1):
        for par, par_next in zip(list_parameters[i], list_parameters[i+1]):
            for k in par.keys():
                if par[k] != par_next[k]:
                    if changed[i] is not None:
                        print('Warning! It went wrong in the cache cleaning')
                    changed[i] = (k, par[k])
                    break
            if changed[i] is not None:
                break
    # TODO: handle the special case of the last element
    return changed
