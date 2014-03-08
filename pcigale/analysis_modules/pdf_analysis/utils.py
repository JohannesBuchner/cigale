# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2014 Yannick Roehlly <yannick@iaora.eu>
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

import numpy as np
from scipy.stats import gaussian_kde
from scipy.linalg import LinAlgError


def gen_pdf(values, probabilities, grid):
    """Generate a probability density function

    For a list of values and associated probabilities, this function
    generates a probability density using a weighted gaussian kernel
    density estimation.

    This part should definitely be improved as it is done in the simplest
    way: each value is repeated (probability * 100) times and the standard
    scipy gaussian KDE is used.

    Parameters
    ----------
    values : array like of floats
        The values of the variable.
    probabilities : array like of floats
        The probability associated with each value
    grid : array like of float
        The list of values to which the probability will be evaluated.

    Returns
    -------
    The list of probabilities evaluated at each value of the grid. If there
    is only one input value with a probability superior to 0, the scipy KDE
    algorithm will fail and None is returned.

    """

    # We use masked arrays because in the analysis module the arrays are
    # already masked and the mask is important.
    values = np.ma.array(values, dtype=float)
    probabilities = np.ma.array(probabilities, dtype=float)

    probabilities = np.ma.array(np.around(100. * probabilities),
                                dtype=int, copy=True)

    combined_values = []

    # We must convert the values masked array to list and test each value
    # against None because 0. is a valid value. For the probabilities,
    # we don't have this problem because if the probability is 0 we won't add
    # the value.
    for val_idx, val in enumerate(values.tolist()):
        if val is not None and probabilities[val_idx]:
            combined_values += [val] * probabilities[val_idx]

    try:
        result = gaussian_kde(combined_values)(grid)
    except LinAlgError:
        result = None

    return result
