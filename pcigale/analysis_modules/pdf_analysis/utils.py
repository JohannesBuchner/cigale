# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2014 Yannick Roehlly <yannick@iaora.eu>
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

import numpy as np
from scipy.stats import gaussian_kde
from scipy.linalg import LinAlgError
from matplotlib import pyplot as plt
from copy import deepcopy
from ...sed.cosmology import cosmology
from ...creation_modules import get_module as get_creation_module


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


def gen_best_sed_fig(wave, fnu, filters_wave, filters_model, filters_obs):
    """Generate a figure for plotting the best models

    Parameters
    ----------
    wave : array-like of floats
        The wavelength grid of the model spectrum.
    fnu : array-like of floats
        The Fnu spectrum of the model at each wavelength.
    filters_wave : array-like of floats
        The effective wavelengths of the various filters.
    filters_model : array-like of floats
        The model fluxes in each filter.
    filters_obs : array-like of floats
        The observed fluxes in each filter.

    Returns
    -------
    A matplotlib.plt.figure.

    """

    try:
        figure = plt.figure()
        ax = figure.add_subplot(111)
        ax.loglog(wave, fnu, "-b", label="Model spectrum")
        ax.loglog(filters_wave, filters_model, "ob", label="Model fluxes")
        ax.loglog(filters_wave, filters_obs, "or", label="Observation fluxes")
        ax.set_xlabel("Wavelength [nm]")
        ax.set_ylabel("Flux [mJy]")
        ax.legend(loc=0)

        return figure

    except ValueError:
        # If the SED can't be plot in x and y logarithm scaled, that means
        # that we have either negative wavelength or flux and that something
        # has gone wrong.
        return None
