# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013-2014 Institute of Astronomy
# Copyright (C) 2014 Laboratoire d'Astrophysique de Marseille, AMU
# Copyright (C) 2014 Yannick Roehlly <yannick@iaora.eu>
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly, Médéric Boquien & Denis Burgarella

import time

import numpy as np
from scipy.special import erf
import scipy.stats as st

from .utils import (save_best_sed, save_pdf, save_chi2, compute_chi2,
                    weighted_param)
from ...warehouse import SedWarehouse

# Probability threshold: models with a lower probability  are excluded from the
# moments computation.
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
    filters: List
        Contains the names of the filters to compute the fluxes.
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

    gbl_warehouse = SedWarehouse()


def init_analysis(params, filters, analysed, redshifts, fluxes, variables,
                  t_begin, n_computed, analysed_averages, analysed_std,
                  best_fluxes, best_parameters, best_chi2, best_chi2_red, save,
                  lim_flag, n_obs):
    """Initializer of the pool of processes. It is mostly used to convert
    RawArrays into numpy arrays. The latter are defined as global variables to
    be accessible from the workers.

    Parameters
    ----------
    params: ParametersHandler
        Handles the parameters from a 1D index.
    filters: list
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
    global gbl_redshifts, gbl_analysed_averages, gbl_analysed_std
    global gbl_best_fluxes, gbl_best_parameters, gbl_best_chi2
    global gbl_best_chi2_red, gbl_save, gbl_n_obs, gbl_lim_flag, gbl_keys

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

    gbl_redshifts = gbl_model_redshifts[np.unique(gbl_model_redshifts,
                                                  return_index=True)[1]]

    gbl_save = save
    gbl_lim_flag = lim_flag

    gbl_n_obs = n_obs
    gbl_keys = None


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

    if 'sfh.age' in sed.info and sed.info['sfh.age'] > sed.info['universe.age']:
        model_fluxes = np.full(len(gbl_filters), np.nan)
        model_variables = np.full(len(gbl_analysed_variables), np.nan)
    else:
        model_fluxes = np.array([sed.compute_fnu(filter_) for filter_ in
                                 gbl_filters])
        model_variables = np.array([sed.info[name]
                                    for name in gbl_analysed_variables])

    gbl_model_redshifts[idx] = sed.info['universe.redshift']
    gbl_model_fluxes[idx, :] = model_fluxes
    gbl_model_variables[idx, :] = model_variables

    with gbl_n_computed.get_lock():
        gbl_n_computed.value += 1
        n_computed = gbl_n_computed.value
    if n_computed % 250 == 0 or n_computed == gbl_params.size:
        t_elapsed = time.time() - gbl_t_begin
        print("{}/{} models computed in {} seconds ({} models/s)".
              format(n_computed, gbl_params.size,
                     np.around(t_elapsed, decimals=1),
                     np.around(n_computed/t_elapsed, decimals=1)),
              end="\n" if n_computed == gbl_params.size else "\r")


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
    np.seterr(invalid='ignore')

    # Tolerance threshold under which any flux or error is considered as 0.
    global gbl_keys
    tolerance = 1e-12

    obs_fluxes = np.array([obs[name] for name in gbl_filters])
    obs_errors = np.array([obs[name + "_err"] for name in gbl_filters])
    nobs = np.where(np.isfinite(obs_fluxes))[0].size

    # We pick the the models with the closest redshift using a slice to work on
    # views of the arrays and not on copies to save on RAM.
    wz = slice(np.abs(obs['redshift'] - gbl_redshifts).argmin(), None,
               gbl_redshifts.size)

    chi2, scaling = compute_chi2(gbl_model_fluxes[wz, :], obs_fluxes,
                                 obs_errors, gbl_lim_flag)

    ##################################################################
    # Variable analysis                                              #
    ##################################################################

    # We select only models that have at least 0.1% of the probability of
    # the best model to reproduce the observations. It helps eliminating
    # very bad models.
    maxchi2 = st.chi2.isf(st.chi2.sf(np.min(chi2), nobs - 1) * 1e-3, nobs - 1)
    wlikely = np.where(chi2 < maxchi2)

    if wlikely[0].size == 0:
        # It sometimes happen because models are older than the Universe's age
        print("No suitable model found for the object {}. One possible origin "
              "is that models are older than the Universe.".format(obs['id']))
        gbl_analysed_averages[idx, :] = np.nan
        gbl_analysed_std[idx, :] = np.nan
        gbl_best_fluxes[idx, :] = np.nan
        gbl_best_parameters[idx, :] = np.nan
        gbl_best_chi2[idx] = np.nan
        gbl_best_chi2_red[idx] = np.nan
    else:
        # We use the exponential probability associated with the χ² as
        # likelihood function.
        likelihood = np.exp(-chi2[wlikely]/2)

        # We define the best fitting model for each observation as the one
        # with the least χ².
        best_index_z = chi2.argmin()  # index for models at given z
        best_index = wz.start + best_index_z * wz.step  # index for all models

        # We compute once again the best sed to obtain its info
        global gbl_previous_idx
        if gbl_previous_idx > -1:
            gbl_warehouse.partial_clear_cache(
                gbl_params.index_module_changed(gbl_previous_idx,
                                                best_index))
        gbl_previous_idx = best_index

        sed = gbl_warehouse.get_sed(gbl_params.modules,
                                    gbl_params.from_index(best_index))

        # We correct the mass-dependent parameters
        for key in sed.mass_proportional_info:
            sed.info[key] *= scaling[best_index_z]

        # We compute the weighted average and standard deviation using the
        # likelihood as weight.
        for i, variable in enumerate(gbl_analysed_variables):
            if variable in sed.mass_proportional_info:
                mean, std = weighted_param(gbl_model_variables[wz, i][wlikely]
                                           * scaling[wlikely], likelihood)
            else:
                mean, std = weighted_param(gbl_model_variables[wz, i][wlikely],
                                           likelihood)
            gbl_analysed_averages[idx, i] = mean
            gbl_analysed_std[idx, i] = max(0.05 * mean, std)

        gbl_best_fluxes[idx, :] = gbl_model_fluxes[best_index, :] \
            * scaling[best_index_z]

        if gbl_keys is None:
            gbl_keys = list(sed.info.keys())
            gbl_keys.sort()
        gbl_best_parameters[idx, :] = np.array([sed.info[k] for k in gbl_keys])
        gbl_best_chi2[idx] = chi2[best_index_z]
        gbl_best_chi2_red[idx] = chi2[best_index_z] / (nobs - 1)

        if gbl_save['best_sed']:
            save_best_sed(obs['id'], sed, scaling[best_index_z])
        if gbl_save['chi2']:
            for i, variable in enumerate(gbl_analysed_variables):
                if variable in sed.mass_proportional_info:
                    save_chi2(obs['id'], variable, gbl_model_variables[wz, i] *
                              scaling, chi2 / (nobs - 1))
                else:
                    save_chi2(obs['id'], variable, gbl_model_variables[wz, i],
                              chi2 / (nobs - 1))
        if gbl_save['pdf']:
            for i, variable in enumerate(gbl_analysed_variables):
                if variable in sed.mass_proportional_info:
                    save_pdf(obs['id'], variable,
                             gbl_model_variables[wz, i][wlikely] *
                             scaling[wlikely], likelihood)
                else:
                    save_pdf(obs['id'], variable,
                             gbl_model_variables[wz, i][wlikely], likelihood)

    with gbl_n_computed.get_lock():
        gbl_n_computed.value += 1
        n_computed = gbl_n_computed.value
    t_elapsed = time.time() - gbl_t_begin
    print("{}/{} objects analysed in {} seconds ({} objects/s)".
          format(n_computed, gbl_n_obs, np.around(t_elapsed, decimals=1),
                 np.around(n_computed/t_elapsed, decimals=2)),
          end="\n" if n_computed == gbl_n_obs else "\r")
