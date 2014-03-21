# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013-2014 Institute of Astronomy
# Copyright (C) 2014 Yannick Roehlly <yannick@iaora.eu>
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly & Médéric Boquien

import numpy as np
from .utils import save_best_sed, save_pdf, save_chi2
import pcigale.analysis_modules.myglobals as gbl
import time


# Probability threshold: models with a lower probability are excluded from
# the moments computation.
MIN_PROBABILITY = 1e-20


def sed(model_params, changed):
    """Worker process to retrieve a SED and return the relevant data

    Parameters
    ----------
    model_params: list
        Parameters of the creation modules

    Returns
    -------
    The model fluxes in each band, the values of the analysed variables, the
    redshift and all the related information. Note, we could retrieve the
    values of the analysed variables afterwards but this would be done in the
    main process. We should benchmark this to see whether it makes any
    significant difference as returning onlt the sed.info would make things
    cleaner here and no more dirty on the caller's side (as we have to unpack
    the list returned by starmap anyway.

    """
    sed = gbl.warehouse.get_sed(gbl.creation_modules, model_params)
    gbl.warehouse.partial_clear_cache(changed)

    if 'age' in sed.info and sed.info['age'] > sed.info['universe.age']:
        model_fluxes = -99. * np.ones(len(gbl.filters))
        model_variables = -99. * np.ones(len(gbl.analysed_variables))
    else:
        model_fluxes = np.array([sed.compute_fnu(filter_.trans_table,
                                                 filter_.effective_wavelength)
                                 for filter_ in gbl.filters.values()])
        model_variables = np.array([sed.info[name]
                                    for name in gbl.analysed_variables])
    redshift = sed.info['redshift']
    info = list(sed.info.values()) # Prevents from returning a view

    with gbl.n_computed.get_lock():
        gbl.n_computed.value += 1
        n_computed = gbl.n_computed.value
    if n_computed % 100 == 0 or n_computed == gbl.n_models:
        t_elapsed = time.time() - gbl.t_begin
        print("{}/{} models computed in {} seconds ({} models/s)".
              format(gbl.n_computed.value, gbl.n_models,
                     np.around(t_elapsed, decimals=1),
                     np.around(n_computed/t_elapsed, decimals=1)),
              end="\r")

    return model_fluxes, model_variables, redshift, info


def analysis(obs):
    """Worker process to analyse the PDF and estimate parameters values

    Parameters
    ----------
    obs: row
        Input data for an individual object

    Returns
    -------
    The analysed parameters (values+errors), best raw and reduced χ², best
    normalisation factor, info of the best SED, fluxes of the best SED
    """

    w = np.where(gbl.w_redshifts[gbl.redshifts[np.abs(obs['redshift'] -
                                               gbl.redshifts).argmin()]])

    model_fluxes = gbl.model_fluxes[w[0], :]
    model_variables = gbl.model_variables[w[0], :]

    obs_fluxes = np.array([obs[name] for name in gbl.filters])
    obs_errors = np.array([obs[name + "_err"] for name in gbl.filters])

    # Some observations may not have flux value in some filters, in
    # that case the user is asked to put -9999 as value. We mask these
    # values. Note, we must mask obs_fluxes after obs_errors.
    obs_errors = np.ma.masked_where(obs_fluxes < -9990., obs_errors)
    obs_fluxes = np.ma.masked_less(obs_fluxes, -9990.)

    # Normalisation factor to be applied to a model fluxes to best fit
    # an observation fluxes. Normalised flux of the models. χ² and
    # likelihood of the fitting. Reduced χ² (divided by the number of
    # filters to do the fit).
    norm_facts = (
        np.sum(model_fluxes * obs_fluxes / (obs_errors * obs_errors), axis=1) /
        np.sum(model_fluxes * model_fluxes / (obs_errors * obs_errors), axis=1)
    )
    norm_model_fluxes = model_fluxes * norm_facts[:, np.newaxis]

    # χ² of the comparison of each model to each observation.
    chi2_ = np.sum(
        np.square((obs_fluxes - norm_model_fluxes) / obs_errors),
        axis=1)

    # We define the reduced χ² as the χ² divided by the number of
    # fluxes used for the fitting.
    chi2_red = chi2_ / obs_fluxes.count()

    # We use the exponential probability associated with the χ² as
    # likelihood function.
    # WARNING: shouldn't this be chi2_red rather?
    likelihood = np.exp(-chi2_/2)
    # For the analysis, we consider that the computed models explain
    # each observation. We normalise the likelihood function to have a
    # total likelihood of 1 for each observation.
    likelihood /= np.sum(likelihood)
    # We don't want to take into account the models with a probability
    # less that the threshold.
    likelihood = np.ma.masked_less(likelihood, MIN_PROBABILITY)
    # We re-normalise the likelihood.
    likelihood /= np.sum(likelihood)

    # We take the mass-dependent variable list from the last computed
    # sed.
    for index, variable in enumerate(gbl.analysed_variables):
        if variable in gbl.mass_proportional_info:
            model_variables[:, index] *= norm_facts

    ##################################################################
    # Variable analysis                                              #
    ##################################################################

    # We compute the weighted average and standard deviation using the
    # likelihood as weight. We first build the weight array by
    # expanding the likelihood along a new axis corresponding to the
    # analysed variable.
    weights = likelihood[:, np.newaxis].repeat(len(gbl.analysed_variables),
                                               axis=1)

    # Analysed variables average and standard deviation arrays.
    analysed_averages = np.ma.average(model_variables, axis=0,
                                      weights=weights)

    analysed_std = np.ma.sqrt(np.ma.average(
        (model_variables - analysed_averages[np.newaxis, :])**2, axis=0,
        weights=weights))

    # We define the best fitting model for each observation as the one
    # with the least χ².
    best_index = chi2_.argmin()

    if gbl.save_best_sed:
        save_best_sed(obs['id'], gbl.creation_modules,
                      gbl.creation_modules_params[w[0][best_index]],
                      norm_facts[best_index])
    if gbl.save_chi2:
        save_chi2(obs['id'], gbl.analysed_variables, model_variables, chi2_red)
    if gbl.save_pdf:
        save_pdf(obs['id'], gbl.analysed_variables, model_variables,
                 likelihood)

    with gbl.n_computed.get_lock():
        gbl.n_computed.value += 1
        n_computed = gbl.n_computed.value
    if n_computed % 100 == 0 or n_computed == gbl.n_obs:
        t_elapsed = time.time() - gbl.t_begin
        print("{}/{} objects analysed in {} seconds ({} objects/s)".
              format(gbl.n_computed.value, gbl.n_obs,
                     np.around(t_elapsed, decimals=1),
                     np.around(n_computed/t_elapsed, decimals=1)),
              end="\r")

    return (analysed_averages,
            analysed_std,
            np.array(model_fluxes[best_index, :]),  # do NOT remove np.array()
            list(gbl.model_info[w[0][best_index]]),
            norm_facts[best_index],
            chi2_[best_index],
            chi2_red[best_index])
