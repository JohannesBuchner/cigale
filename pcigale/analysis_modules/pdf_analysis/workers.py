# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013-2014 Institute of Astronomy
# Copyright (C) 2014 Yannick Roehlly <yannick@iaora.eu>
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly & Médéric Boquien

import numpy as np
from .utils import save_best_sed, save_pdf, save_chi2

# Probability threshold: models with a lower probability are excluded from
# the moments computation.
MIN_PROBABILITY = 1e-20


def sed(warehouse, creation_modules, model_params, analysed_variables,
        filters):
    """Worker process to retrieve a SED and return the relevant data

    Parameters
    ----------
    warehouse: SedWarhouse object
        Used to retrieve a SED. Ideally a different warehouse should be used
        for each forked process but it does not seem to cause any problem so
        far
    creation_modules: list
        List of creation modules to build the SED
    model_params: list
        Parameters of the creation modules
    analysed_variables: list
        Names of the analysed variables
    filters: list
        Filters to compute the SED fluxes

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
    sed = warehouse.get_sed(creation_modules, model_params)
    if 'age' in sed.info and sed.info['age'] > sed.info['universe.age']:
        model_fluxes = -99. * np.ones(len(filters))
        model_variables = -99. * np.ones(len(analysed_variables))
    else:
        model_fluxes = np.array([sed.compute_fnu(filter_.trans_table,
                                                 filter_.effective_wavelength)
                                 for filter_ in filters.values()])
        model_variables = np.array([sed.info[name]
                                    for name in analysed_variables])
    redshift = sed.info['redshift']
    info = sed.info.values()

    return model_fluxes, model_variables, redshift, info


def analysis(obs, model_fluxes, model_variables, info, filters, sed,
             analysed_variables, creation_modules, creation_modules_params,
             save):
    """Worker process to analyse the PDF and estimate parameters values

    Parameters
    ----------
    obs: row
        Input data for an individual object
    model_fluxes: 2D array
        Fluxes for each model and for each filter
    model_variables: 2D array
        Variables values for each model
    info: list
        sed.info for each model
    filters: list
        Filters to compute the fluxes
    sed: SED object
        Used to retrieve which parameters are normalisation dependent
    analysed_variables: list
        Names of analysed variables
    creation_modules: list
        Creation modules named to recreate the best SED
    creation_modules_params: list
        Creation modules parameters to recreate the best SED
    save: set
        Booleans indicating whether to save the best SED, best χ² and PDF

    Returns
    -------
    The analysed parameters (values+errors), best raw and reduced χ², best
    normalisation factor, info of the best SED, fluxes of the best SED
    """
    obs_fluxes = np.array([obs[name] for name in filters])
    obs_errors = np.array([obs[name + "_err"] for name in filters])

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
    for index, variable in enumerate(analysed_variables):
        if variable in sed.mass_proportional_info:
            model_variables[:, index] *= norm_facts

    # We also add the galaxy mass to the analysed variables if relevant
    if sed.sfh is not None:
        analysed_variables.insert(0, "galaxy_mass")
        model_variables = np.dstack((norm_facts,
                                     model_variables))

    ##################################################################
    # Variable analysis                                              #
    ##################################################################

    # We compute the weighted average and standard deviation using the
    # likelihood as weight. We first build the weight array by
    # expanding the likelihood along a new axis corresponding to the
    # analysed variable.
    weights = likelihood[:, np.newaxis].repeat(len(analysed_variables), axis=1)

    # Analysed variables average and standard deviation arrays.
    analysed_averages = np.ma.average(model_variables, axis=0,
                                      weights=weights)

    analysed_std = np.ma.sqrt(np.ma.average(
        (model_variables - analysed_averages[np.newaxis, :])**2, axis=0,
        weights=weights))

    # We define the best fitting model for each observation as the one
    # with the least χ².
    best_index = chi2_.argmin()

    if save[0]:
        save_best_sed(obs['id'], creation_modules,
                      creation_modules_params[best_index],
                      norm_facts[best_index])
    if save[1]:
        save_chi2(obs['id'], analysed_variables, model_variables, chi2_red)
    if save[2]:
        save_pdf(obs['id'], analysed_variables, model_variables, likelihood)

    return (analysed_averages,
            analysed_std,
            chi2_[best_index],
            chi2_red[best_index],
            norm_facts[best_index],
            list(info[best_index]),
            model_fluxes[best_index, :])
