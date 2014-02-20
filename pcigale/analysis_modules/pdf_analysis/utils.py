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


def gen_compute_fluxes_at_redshift(sed, filters, redshifting_module,
                                   igm_module):
    """"Generate function to compute the fluxes of a SED at a given redshift

    Given a SED, a list of filters and a redshift module name, this generator
    returns a function computing the fluxes of the SED in all the filters at
    a given redshift.  If the SED is older than the Universe at the given
    redshift, the fluxes returned by this function will be all -99.

    If the redshift module name is None, the SED is not redshifted by the
    returned function. This means that it will always return the same fluxes,
    whatever its redshift parameter. This can be used for computing
    photometric redshifts.

    Parameters
    ----------
    sed : pcigale.sed
        The pcigale SED object.
    filters : list of pcigale.data.filters
        List of pcigale filters objects.
    redshifting_module : picgale.creation_modules.Module
        A pcigale SED creation module to redshift the SED. It must accept a
        redshift parameter (or be None).
    igm_module : picgale.creation_modules.Module
        A pcigale SED creation module to add the IGM attenuation to the SED.
        Is not used if the redshifting module is None.

    Return
    ------
    gen_fluxes : function
        Function computing the fluxes of the SED in all filters at a given
        redshift.

    """
    # The returned function is memoized.
    cache = {}

    def gen_fluxes(redshift):
        """Compute the flux of the SED in various filters.

        Parameters
        ----------
        redshift : float

        Returns
        -------
        array fo floats

        """

        # If the function is generated without a redshift module, it always
        # computes the fluxes at the SED redshift (the SED may be already
        # redshifted).
        if not redshifting_module:
            redshift = sed.redshift

        if redshift not in cache:
            # Age of the Universe at redshift.
            # astropy 0.3 cosmology functions return quantities
            try:
                age_at_redshift = cosmology.age(redshift).value * 1000
            except AttributeError:
                age_at_redshift = cosmology.age(redshift) * 1000

            if sed.info["age"] > age_at_redshift:
                cache[redshift] = -99 * np.ones(len(filters))
            else:

                if redshifting_module:
                    red_sed = deepcopy(sed)
                    redshifting_module.parameters["redshift"] = redshift
                    redshifting_module.process(red_sed)
                    if igm_module:
                        igm_module.process(red_sed)
                else:
                    red_sed = sed

                cache[redshift] = np.array(
                    [red_sed.compute_fnu(filt_.trans_table,
                                         filt_.effective_wavelength)
                     for filt_ in filters])

        return cache[redshift]

    return gen_fluxes


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
