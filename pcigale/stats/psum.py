# -*- coding: utf-8 -*-
"""
Copyright (C) 2012 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

This file implements the statistical analysis as performed by the calcX2_psum
programme of the Fortran Cigale code.

"""

import numpy as np
from . import common


# Tolerance threshold under which any flux or error is considered as 0.
TOLERANCE = 1.e-12


class Module(common.AnalysisModule):
    """psum analysis

    TODO: Description of the PSUM method.
    """

    paramtre_list = {}


def adjust_errors(flux, error, default_error=0.1, systematic_deviation=0.1):
    """Adjust the errors replacing the 0 values by the default error and
    adding the systematic deviation.

    The systematic deviation change the error to:
    sqrt( error² + (flux * deviation)² )

    Parametres
    ----------
    flux : array of floats
        Fluxes.
    error : array of floats
        Observational error in the same unit as the fluxes.
    default_error : float
        Default error factor used when the provided error in under the
        tolerance threshold.
    systematic_deviation : float
        Systematic deviation added to the error.

    Returns
    -------
    error : array of floats
        The corrected errors.

    """

    # The arrays must have the same lengths.
    if len(flux) != len(error):
        raise ValueError("The flux and error arrays must have the same "
                         "length.")

    # We copy the error array not to modify the original one.
    error = np.copy(error)

    # Replace errors below tolerance by the default one.
    error[error < TOLERANCE] = (default_error * error[error < TOLERANCE])

    # Add the systematic error.
    error = np.sqrt(np.square(error) + np.square(flux * systematic_deviation))

    return error


def compute_chi2(model_fluxes, obs_fluxes, obs_errors):
    """Compute chi square value and normalisation factor for the comparison
    of a model fluxes to observational ones.

    Parametres
    ----------
    model_fluxes : array of floats
        Model fluxes.
    obs_fluxes : array of floats
        Observation fluxes for the same filters as the model ones and
        in the same unit.
    obs_errors : array of floats
        Error the observation flux. The error must be Gaussian for the
        chi-square to be meaning full.

    Returns
    -------
    chi2_reduced : float
        Reduced chi square value for the comparison. The maximum Chi square
        value returned is 99.
    normalisation_factor : float
        Normalisation factor that must be applied to the model to fit the
        observation.

    """

    # The three arrays must have the same length.
    if (len(model_fluxes) != len(obs_fluxes) or
            len(obs_fluxes) != len(obs_errors)):
        raise ValueError("The model fluxes, observation fluxes and "
                         "observation errors arrays must have the "
                         "same length.")

    # We copy the dictionaries not to modify the original ones.
    model_fluxes = np.copy(model_fluxes)
    obs_fluxes = np.copy(obs_fluxes)
    obs_errors = np.copy(obs_errors)

    # If no observed flux is over the tolerance threshold, or if any error,
    # for valid fluxes, is under the threshold then the observation is set
    # as not fitting at all.
    if (max(obs_fluxes) < TOLERANCE or
            min(obs_errors[obs_fluxes > TOLERANCE]) < TOLERANCE):
        reduced_chi2 = 99
        normalisation_factor = 1
    else:
        # We make the computation using only the filters for which the
        # observation error is over the tolerance threshold.
        (model_fluxes, obs_fluxes, obs_errors) = \
            (model_fluxes[obs_errors > TOLERANCE],
             obs_fluxes[obs_errors > TOLERANCE],
             obs_errors[obs_errors > TOLERANCE])

        #FIXME
        degrees_of_freedom = len(model_fluxes) - 1

        if degrees_of_freedom == 0:
            #FIXME
            reduced_chi2 = 0
            normalisation_factor = sum(obs_fluxes) / sum(model_fluxes)
        else:
            normalisation_factor = (sum(obs_fluxes * model_fluxes) /
                                    sum(model_fluxes * model_fluxes))
            norm_model_fluxes = normalisation_factor * model_fluxes
            reduced_chi2 = (sum(np.square((obs_fluxes - norm_model_fluxes) /
                            obs_errors))
                            / degrees_of_freedom)
            reduced_chi2 = min(reduced_chi2, 99)

    return reduced_chi2, normalisation_factor
