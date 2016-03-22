# -*- coding: utf-8 -*-
# Copyright (C) 2016 Universidad de Antofagasta
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Médéric Boquien

from collections import Iterable
import numpy as np
import validate as vdt


class VdtValueRepeatedError(vdt.VdtValueError):
    """The value supplied was repeated."""

    def __init__(self, value):
        """
        >>> raise VdtValueRepeatedError('jedie')
        Traceback (most recent call last):
        VdtValueRepeatedError: the value "jedie" is too long.
        """
        vdt.ValidateError.__init__(self, 'the value {} is repeated.'.format(value))


def is_cigale_list(inobject, dtype='float', minvalue=None, maxvalue=None,
                   options=None):
    """Function that returns a list made of float or int numpy. It also checks
    if any element is lower/high then a set minimum/maximum value or if any
    does not belong to a predefined list. If any of the tests fails, an
    exception is thrown. All the parameters besides "value" are strings because
    configobj/validate pass arguments as such.

    Parameters
    ----------
    inobject: list or string
        Object that we need to convert and verify
    dtype: string
        Desired type of the object
    minvalue: string
        Minimum value of all the elements
    maxvalue: string
        Maximum value of all the elements
    options: string
        All the values that the array can take

    Returns
    -------
    outobject: list
        List of the requested type and respecting the constraints
    """

    if dtype == 'float':
        dtype = float
        conv_type = lambda x: float(x)
    elif dtype == 'int':
        dtype = int
        conv_type = lambda x: int(float(x))
    else:
        raise Exception("Unsupported data type. Only float and int are "
                        "supported.")

    if isinstance(inobject, str):
        if inobject.startswith('eval '):
            outobject = eval(inobject[4:])
            # If the evaluation lead to a single value, we put it in a list.
            if not isinstance(outobject, Iterable):
                outobject = [outobject]
        elif inobject.startswith('range '):
            start, stop, step = [conv_type(item) for item in inobject[5:].split()]
            outobject = np.arange(start, stop+step, step)
        else:
            # We need to return a list to combine the list of possible values
            # for each parameter.
            outobject = [inobject]
    else:
        outobject = inobject

    try:  # We convert to an array as it is more convenient for tests
        outobject = np.asarray(outobject, dtype)
    except ValueError:
        raise vdt.VdtValueError(inobject)

    # Test if there are two identical elements
    unique_outobject, counts = np.unique(outobject, return_counts=True)
    if outobject.size > unique_outobject.size:
        raise VdtValueRepeatedError(unique_outobject[counts > 1])

    # Test if any value is smaller than the minimum value
    if minvalue is not None:
        minvalue = conv_type(minvalue)
        if np.any(outobject < minvalue):
            raise vdt.VdtValueTooSmallError(outobject[outobject < minvalue])

    # Test if any value is larger than the maximum value
    if maxvalue is not None:
        maxvalue = conv_type(maxvalue)
        if np.any(outobject > maxvalue):
            raise vdt.VdtValueTooBigError(outobject[outobject > maxvalue])

    # Test if any value is not in the list of allowed options
    if options is not None:
        options = {conv_type(option) for option in options.split('&')}
        if not set(outobject).issubset(options):
            raise vdt.VdtValueError(set(outobject) - options)

    return list(outobject)


def is_cigale_string_list(inobject):
    """Function that returns a list of strings. If inobject is already a list,
    it is returned as is. Otherwise the string is put into a list. If the
    string is empty, it returns an empty list.

    Parameters
    ----------
    inobject: list or string
        Object that we need to convert and verify

    Returns
    -------
    outobject: list
        Object transformed into a list

    """
    if isinstance(inobject, str):
        if len(inobject) == 0:
            return []
        else:
            return [inobject]

    return inobject


# Dictionary of functions that are to be passed to the validator so it knows
# about extra tests that are not covered with the base functions
functions = {'cigale_list': is_cigale_list,
             'cigale_string_list': is_cigale_string_list}
