# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013-2014 Institute of Astronomy
# Copyright (C) 2014 Yannick Roehlly <yannick@iaora.eu>
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly & Médéric Boquien

import time

import numpy as np

from .utils import save_best_sed, save_pdf, save_chi2
from ...warehouse import SedWarehouse

# Probability threshold: models with a lower probability are excluded from
# the moments computation.
MIN_PROBABILITY = 1e-20


def init_sed(params, filters, analysed, redshifts, fluxes, variables,
             t_begin, n_computed):
    """Initializer of the pool of processes. It is mostly used to convert
    RawArrays into numpy arrays. The latter are defined as global variables to
    be accessible from the workers.

    Parameters
    ----------
    params: ParametersHandler
        Handles the parameters from a 1D index.
    filters: OrderedDict
        Contains filters to compute the fluxes.
    analysed: list
        Variable names to be analysed.
    redshifts: RawArray and tuple containing the shape
        Redshifts of individual models. Shared among workers.
    fluxes: RawArray and tuple containing the shape
        Fluxes of individual models. Shared among workers.
    variables: RawArray and tuple containing the shape
        Values of the analysed variables. Shared among workers.
    n_computed: Value
        Number of computed models. Shared among workers.
    t_begin: float
        Time of the beginning of the computation.

    """
    global gbl_model_redshifts, gbl_model_fluxes, gbl_model_variables
    global gbl_n_computed, gbl_t_begin, gbl_params, gbl_previous_idx
    global gbl_filters, gbl_analysed_variables, gbl_warehouse

    gbl_model_redshifts = np.ctypeslib.as_array(redshifts[0])

    gbl_model_fluxes = np.ctypeslib.as_array(fluxes[0])
    gbl_model_fluxes = gbl_model_fluxes.reshape(fluxes[1])

    gbl_model_variables = np.ctypeslib.as_array(variables[0])
    gbl_model_variables = gbl_model_variables.reshape(variables[1])

    gbl_n_computed = n_computed
    gbl_t_begin = t_begin

    gbl_params = params

    gbl_previous_idx = -1

    gbl_filters = filters
    gbl_analysed_variables = analysed

    gbl_warehouse = SedWarehouse(cache_type="memory")

def init_analysis(params, filters, analysed, redshifts, fluxes, variables,
                  t_begin, n_computed, analysed_averages, analysed_std,
                  best_fluxes, best_parameters, best_chi2, best_chi2_red, save,
                  n_obs):
    """Initializer of the pool of processes. It is mostly used to convert
    RawArrays into numpy arrays. The latter are defined as global variables to
    be accessible from the workers.

    Parameters
    ----------
    params: ParametersHandler
        Handles the parameters from a 1D index.
    filters: OrderedDict
        Contains filters to compute the fluxes.
    analysed: list
        Variable names to be analysed
    redshifts: RawArray and tuple containing the shape.
        Redshifts of individual models. Shared among workers.
    fluxes: RawArray
        Fluxes of individual models. Shared among workers.
    variables: RawArray and tuple containing the shape
        Values of the analysed variables. Shared among workers.
    t_begin: float
        Time of the beginning of the computation.
    n_computed: Value
        Number of computed models. Shared among workers.
    analysed_averages: RawArray
        Analysed values for each observation.
    analysed_std: RawArray
        Standard devriation values for each observation.
    best_fluxes: RawArray
        Best fluxes for each observation.
    best_parameters: RawArray
        Best parameters for each observation.
    best_chi2: RawArray
        Best χ² for each observation.
    best_chi2_red: RawArray
        Best reduced χ² for each observation.
    save: dictionary
        Contains booleans indicating whether we need to save some data related
        to given models.
    n_obs: int
        Number of observations.

    """
    init_sed(params, filters, analysed, redshifts, fluxes, variables,
             t_begin, n_computed)
    global gbl_redshifts, gbl_w_redshifts, gbl_analysed_averages
    global gbl_analysed_std, gbl_best_fluxes, gbl_best_parameters
    global gbl_best_chi2, gbl_best_chi2_red, gbl_save, gbl_n_obs

    gbl_analysed_averages = np.ctypeslib.as_array(analysed_averages[0])
    gbl_analysed_averages = gbl_analysed_averages.reshape(analysed_averages[1])

    gbl_analysed_std = np.ctypeslib.as_array(analysed_std[0])
    gbl_analysed_std = gbl_analysed_std.reshape(analysed_std[1])

    gbl_best_fluxes = np.ctypeslib.as_array(best_fluxes[0])
    gbl_best_fluxes = gbl_best_fluxes.reshape(best_fluxes[1])

    gbl_best_parameters = np.ctypeslib.as_array(best_parameters[0])
    gbl_best_parameters = gbl_best_parameters.reshape(best_parameters[1])

    gbl_best_chi2 = np.ctypeslib.as_array(best_chi2[0])

    gbl_best_chi2_red = np.ctypeslib.as_array(best_chi2_red[0])

    gbl_redshifts = np.unique(gbl_model_redshifts)
    gbl_w_redshifts = {redshift: gbl_model_redshifts == redshift
                       for redshift in gbl_redshifts}

    gbl_save = save
    gbl_n_obs = n_obs

def sed(idx):
    """Worker process to retrieve a SED and affect the relevant data to shared
    RawArrays.

    Parameters
    ----------
    idx: int
        Index of the model to retrieve its parameters from the global variable

    """
    global gbl_previous_idx
    if gbl_previous_idx > -1:
        gbl_warehouse.partial_clear_cache(
            gbl_params.index_module_changed(gbl_previous_idx, idx))
    gbl_previous_idx = idx

    sed = gbl_warehouse.get_sed(gbl_params.modules,
                                gbl_params.from_index(idx))

    if 'age' in sed.info and sed.info['age'] > sed.info['universe.age']:
        model_fluxes = -99. * np.ones(len(gbl_filters))
        model_variables = -99. * np.ones(len(gbl_analysed_variables))
    else:
        model_fluxes = np.array([sed.compute_fnu(filter_.trans_table,
                                                 filter_.effective_wavelength)
                                 for filter_ in gbl_filters.values()])
        model_variables = np.array([sed.info[name]
                                    for name in gbl_analysed_variables])

    gbl_model_redshifts[idx] = sed.info['redshift']
    gbl_model_fluxes[idx, :] = model_fluxes
    gbl_model_variables[idx, :] = model_variables

    with gbl_n_computed.get_lock():
        gbl_n_computed.value += 1
        n_computed = gbl_n_computed.value
    if n_computed % 100 == 0 or n_computed == gbl_params.size:
        t_elapsed = time.time() - gbl_t_begin
        print("{}/{} models computed in {} seconds ({} models/s)".
              format(n_computed, gbl_params.size,
                     np.around(t_elapsed, decimals=1),
                     np.around(n_computed/t_elapsed, decimals=1)),
              end="\r")


def analysis(idx, obs):
    """Worker process to analyse the PDF and estimate parameters values

    Parameters
    ----------
    idx: int
        Index of the observation. This is necessary to put the computed values
        at the right location in RawArrays
    obs: row
        Input data for an individual object

    """

    w = np.where(gbl_w_redshifts[gbl_redshifts[np.abs(obs['redshift'] -
                                               gbl_redshifts).argmin()]])

    model_fluxes = np.ma.masked_less(gbl_model_fluxes[w[0], :], -90.)
    model_variables = gbl_model_variables[w[0], :]

    obs_fluxes = np.array([obs[name] for name in gbl_filters])
    obs_errors = np.array([obs[name + "_err"] for name in gbl_filters])

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

    # We use the exponential probability associated with the χ² as
    # likelihood function.
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

    ##################################################################
    # Variable analysis                                              #
    ##################################################################

    # We define the best fitting model for each observation as the one
    # with the least χ².
    best_index = chi2_.argmin()

    # We compute once again the best sed to obtain its info
    global gbl_previous_idx
    if gbl_previous_idx > -1:
        gbl_warehouse.partial_clear_cache(
            gbl_params.index_module_changed(gbl_previous_idx,
                                            w[0][best_index]))
    gbl_previous_idx = w[0][best_index]

    sed = gbl_warehouse.get_sed(gbl_params.modules,
                                gbl_params.from_index([w[0][best_index]]))

    # We correct the mass-dependent parameters
    for key in sed.mass_proportional_info:
        sed.info[key] *= norm_facts[best_index]
    for index, variable in enumerate(gbl_analysed_variables):
        if variable in sed.mass_proportional_info:
            model_variables[:, index] *= norm_facts

    # We compute the weighted average and standard deviation using the
    # likelihood as weight. We first build the weight array by
    # expanding the likelihood along a new axis corresponding to the
    # analysed variable.
    weights = likelihood[:, np.newaxis].repeat(len(gbl_analysed_variables),
                                               axis=1)

    # Analysed variables average and standard deviation arrays.
    analysed_averages = np.ma.average(model_variables, axis=0,
                                      weights=weights)
    analysed_std = np.ma.sqrt(np.ma.average(
        (model_variables - analysed_averages[np.newaxis, :])**2, axis=0,
        weights=weights))

    # TODO Merge with above computation after checking it is fine with a MA.
    gbl_analysed_averages[idx, :] = analysed_averages
    gbl_analysed_std[idx, :] = analysed_std

    gbl_best_fluxes[idx, :] = norm_model_fluxes[best_index, :]
    gbl_best_parameters[idx, :] = list(sed.info.values())
    gbl_best_chi2[idx] = chi2_[best_index]
    gbl_best_chi2_red[idx] = chi2_[best_index] / obs_fluxes.count()

    if gbl_save['best_sed']:
        save_best_sed(obs['id'], sed, norm_facts[best_index])
    if gbl_save['chi2']:
        save_chi2(obs['id'], gbl_analysed_variables, model_variables, chi2_ /
                  obs_fluxes.count())
    if gbl_save['pdf']:
        save_pdf(obs['id'], gbl_analysed_variables, model_variables,
                 likelihood)

    with gbl_n_computed.get_lock():
        gbl_n_computed.value += 1
        n_computed = gbl_n_computed.value
    if n_computed % 100 == 0 or n_computed == gbl_n_obs:
        t_elapsed = time.time() - gbl_t_begin
        print("{}/{} objects analysed in {} seconds ({} objects/s)".
              format(n_computed, gbl_n_obs, np.around(t_elapsed, decimals=1),
                     np.around(n_computed/t_elapsed, decimals=1)),
              end="\r")
