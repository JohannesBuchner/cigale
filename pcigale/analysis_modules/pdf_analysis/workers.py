# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013-2014 Institute of Astronomy
# Copyright (C) 2014 Laboratoire d'Astrophysique de Marseille, AMU
# Copyright (C) 2014 Yannick Roehlly <yannick@iaora.eu>
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly, Médéric Boquien & Denis Burgarella

from textwrap import wrap
import time

import numpy as np
from scipy import optimize
from scipy.special import erf
import scipy.stats as st

from .utils import (save_best_sed, save_pdf, save_chi2, dchi2_over_ds2)
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

    gbl_warehouse = SedWarehouse()

def init_analysis(params, filters, analysed, redshifts, fluxes, variables,
                  t_begin, n_computed, analysed_averages, analysed_std,
                  best_fluxes, best_parameters, best_chi2, best_chi2_red, save,
                  lim_flag, n_obs, phase):
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
    phase: int
        Phase of the analysis (data or mock).

    """
    init_sed(params, filters, analysed, redshifts, fluxes, variables,
             t_begin, n_computed)
    global gbl_redshifts, gbl_w_redshifts, gbl_analysed_averages
    global gbl_analysed_std, gbl_best_fluxes, gbl_best_parameters
    global gbl_best_chi2, gbl_best_chi2_red, gbl_save, gbl_n_obs
    global gbl_lim_flag, gbl_phase

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
    gbl_lim_flag = lim_flag

    gbl_n_obs = n_obs
    gbl_phase = phase

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
    # Tolerance threshold under which any flux or error is considered as 0.
    tolerance = 1e-12

    obs_fluxes = np.array([obs[name] for name in gbl_filters])
    obs_errors = np.array([obs[name + "_err"] for name in gbl_filters])

    wobs = np.where(obs_fluxes > tolerance)
    obs_fluxes = obs_fluxes[wobs]
    obs_errors = obs_errors[wobs]

    # We pick the indices of the models with closest redshift assuming we have
    # limited the number of decimals (usually set to 2 decimals).
    wz = np.where(gbl_w_redshifts[gbl_redshifts[np.abs(obs['redshift'] -
                                               gbl_redshifts).argmin()]])

    # We only keep model with fluxes >= -90. If not => no data
    # Probably because age > age of the universe (see function sed(idx) above).
    model_fluxes = gbl_model_fluxes[wz[0], :]

    model_fluxes = model_fluxes[:, wobs[0]]
    model_variables = gbl_model_variables[wz[0], :]

    wvalid = np.where(model_variables[:, 0] >= -90.)
    model_fluxes = model_fluxes[wvalid[0], :]
    model_variables = model_variables[wvalid[0], :]

    # Some observations may not have flux values in some filter(s), but
    # they can have upper limit(s). To process upper limits, the user
    # is asked to put the upper limit as flux value and an error value with
    # (obs_errors>=-9990. and obs_errors<0.).
    # Next, the user has two options:
    # 1) s/he puts True in the boolean lim_flag
    # and the limits are processed as upper limits below.
    # 2) s/he puts False in the boolean lim_flag
    # and the limits are processed as no-data below.

    lim_flag = gbl_lim_flag and np.any((obs_errors >= -9990.)&
                                       (obs_errors < tolerance))

    # Normalisation factor to be applied to a model fluxes to best fit
    # an observation fluxes. Normalised flux of the models. χ² and
    # likelihood of the fitting. Reduced χ² (divided by the number of
    # filters to do the fit).
    norm_facts = (
        np.sum(model_fluxes * obs_fluxes / (obs_errors * obs_errors), axis=1) /
        np.sum(model_fluxes * model_fluxes / (obs_errors * obs_errors), axis=1)
    )

    if lim_flag is True:
        norm_init = norm_facts
        for imod in range(len(model_fluxes)):
            norm_facts[imod] = optimize.newton(dchi2_over_ds2, norm_init[imod],
                                               tol=1e-16,
                                               args=(obs_fluxes, obs_errors,
                                                     model_fluxes[mod, :]))
    model_fluxes *= norm_facts[:, np.newaxis]

    # χ² of the comparison of each model to each observation.
    if lim_flag is True:
        # This mask selects the filter(s) for which measured fluxes are given
        # i.e., when (obs_flux is >=0. and obs_errors>=0.) and lim_flag=True
        mask_data = (obs_errors >= tolerance)
        # This mask selects the filter(s) for which upper limits are given
        # i.e., when (obs_flux is >=0. (and obs_errors>=-9990., obs_errors<0.))
        # and lim_flag=True
        mask_lim = np.logical_and(obs_errors >= -9990., obs_errors < tolerance)
        chi2_ = np.sum(np.square(
            (obs_fluxes[mask_data]-model_fluxes[:, mask_data]) /
            obs_errors[mask_data]), axis=1)

        chi2_ += -2. * np.sum(
            np.log(
                np.sqrt(np.pi/2.)*(-obs_errors[mask_lim])*(
                    1.+erf(
                        (obs_fluxes[mask_lim]-model_fluxes[:, mask_lim]) /
                        (np.sqrt(2)*(-obs_errors[mask_lim]))))), axis=1)
    else:
        mask_data = np.logical_and(obs_fluxes > tolerance,
                                   obs_errors > tolerance)
        chi2_ = np.sum(np.square(
            (obs_fluxes[mask_data] - model_fluxes[:, mask_data]) /
            obs_errors[mask_data]), axis=1)

    ##################################################################
    # Variable analysis                                              #
    ##################################################################

    # We define the best fitting model for each observation as the one
    # with the least χ².
    if chi2_.size == 0:
        # It sometimes happen because models are older than the Universe's age
        print("No suitable model found for the object {}. One possible origin "
              "is that models are older than the Universe.".format(obs['id']))
    else:
        # We select only models that have at least 0.1% of the probability of the
        # best model to reproduce the observations. It helps eliminating very bad
        # models.
        maxchi2 = st.chi2.isf(st.chi2.sf(np.min(chi2_), obs_fluxes.size-1)*1e-3,
                            obs_fluxes.size-1)
        wlikely = np.where(chi2_ < maxchi2)
        # We use the exponential probability associated with the χ² as
        # likelihood function.
        likelihood = np.exp(-chi2_[wlikely]/2)

        best_index = chi2_.argmin()        

        # We compute once again the best sed to obtain its info
        global gbl_previous_idx
        if gbl_previous_idx > -1:
            gbl_warehouse.partial_clear_cache(
                gbl_params.index_module_changed(gbl_previous_idx,
                                            wz[0][wvalid[0][best_index]]))
        gbl_previous_idx = wz[0][wvalid[0][best_index]]

        sed = gbl_warehouse.get_sed(gbl_params.modules,
                                gbl_params.from_index([wz[0][wvalid[0][best_index]]]))

        # We correct the mass-dependent parameters
        for key in sed.mass_proportional_info:
            sed.info[key] *= norm_facts[best_index]
        for index, variable in enumerate(gbl_analysed_variables):
            if variable in sed.mass_proportional_info:
                model_variables[:, index] *= norm_facts

        # We compute the weighted average and standard deviation using the
        # likelihood as weight.
        analysed_averages = np.empty(len(gbl_analysed_variables))
        analysed_std = np.empty_like(analysed_averages)

        # We check how many unique parameter values are analysed and if less
        # than Npdf (= 100), the PDF is initally built assuming a number of
        # bins equal to the number of unique values for a given parameter
        # (e.g., average_sfr, age, attenuation.uv_bump_amplitude,
        # dust.luminosity, attenuation.FUV, etc.).
        Npdf = 100
        var = np.empty((Npdf, len(analysed_averages)))
        pdf = np.empty((Npdf, len(analysed_averages)))
        min_hist = np.min(model_variables, axis=0)
        max_hist = np.max(model_variables, axis=0)

        for i, val in enumerate(analysed_averages):
            Nhist = min(Npdf, len(np.unique(model_variables[:, i])))

            if min_hist[i] == max_hist[i]:
                analysed_averages[i] = model_variables[0, i]
                analysed_std[i] = 0.

                var[:, i] = max_hist[i]
                pdf[:, i] = 1.
            else:
                pdf_prob, pdf_grid = np.histogram(model_variables[wlikely[0], i],
                                              Nhist,
                                              (min_hist[i], max_hist[i]),
                                              weights=likelihood, density=True)
                pdf_x = (pdf_grid[1:]+pdf_grid[:-1])/2
                pdf_y = pdf_x * pdf_prob
                analysed_averages[i] = np.sum(pdf_y) / np.sum(pdf_prob)
                analysed_std[i] = np.sqrt(
                                     np.sum(
                                    np.square(pdf_x-analysed_averages[i]) * pdf_prob
                                           ) / np.sum(pdf_prob)
                                         )
                analysed_std[i] = max(0.05*analysed_averages[i], analysed_std[i])

                var[:, i] = np.linspace(min_hist[i], max_hist[i], Npdf)
                pdf[:, i] = np.interp(var[:, i], pdf_x, pdf_prob)


    # TODO Merge with above computation after checking it is fine with a MA.
        gbl_analysed_averages[idx, :] = analysed_averages
        gbl_analysed_std[idx, :] = analysed_std

        gbl_best_fluxes[idx, :] = gbl_model_fluxes[wz[0][wvalid[0][best_index]], :] \
                      *norm_facts[best_index]
        gbl_best_parameters[idx, :] = list(sed.info.values())
        gbl_best_chi2[idx] = chi2_[best_index]
        gbl_best_chi2_red[idx] = chi2_[best_index] / obs_fluxes.size

    # If observed SED analysis
        if gbl_phase == 1:
            if gbl_save['best_sed']:
                save_best_sed(obs['id'], sed, norm_facts[best_index])
            if gbl_save['chi2']:
                save_chi2(obs['id'], gbl_analysed_variables, model_variables, chi2_ /
                      obs_fluxes.size)
            if gbl_save['pdf']:
                save_pdf(obs['id'], gbl_analysed_variables, var, pdf)

    with gbl_n_computed.get_lock():
        gbl_n_computed.value += 1
        n_computed = gbl_n_computed.value
    t_elapsed = time.time() - gbl_t_begin
    print("{}/{} objects analysed in {} seconds ({} objects/s)".
            format(n_computed, gbl_n_obs, np.around(t_elapsed, decimals=1),
                    np.around(n_computed/t_elapsed, decimals=2)),
            end="\r")

