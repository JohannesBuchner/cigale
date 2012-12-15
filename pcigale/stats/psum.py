# -*- coding: utf-8 -*-
"""
Copyright (C) 2012 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

This file implements the statistical analysis as performed by the calcX2_psum
programme of the Fortran Cigale code.

The models corresponding to all possible combinations of parametres are
computed are the integrated flux in the same filters as the observations are
used to compute the χ² of the fitting. This χ² give a propability that is
associated with the model values for the parametres. At the end, for each
parametre, the (probability) weighted mean and standard deviation are computed
and the best fitting model (the one with the least reduced χ²) is given.

TODO: Factorise the way the analysis is done to have general methods in the
      AnalysisModule class that are defined in the specific class; in order
      to make the statistical analysis modular as the SED creation modules
      are.

"""

import atpy
import numpy as np
from . import common
from ..sed.warehouse import create_sed
from ..data import Database


# Tolerance threshold under which any flux or error is considered as 0.
TOLERANCE = 1.e-12


class Module(common.AnalysisModule):
    """psum analysis

    TODO: Description of the PSUM method.
    """

    paramtre_list = {}

    # TODO: Don't use the configuration from the pcigale.session but the
    # individual needed components.
    def process(self, data_file, column_list, sed_modules,
                sed_modules_params):
        """Process with the psum analysis.

        Parametres
        ----------
        configuration : dictionary
            Session configuration dictionnary resulting of
            pcigale.session.Configuration.configuration

        Returns
        -------
        results : list of tuples (pcigale.sed object, dict, float, float)
            There is one tuple per observed object: the first element is the
            best fitting SED for this object, the second dictionary of
            parametre used to produce i, the third is the reduced Chi-square
            of the fit and the fourth is the normalisation factor to be
            applied to the SED to fit the observation.

        """
        filter_list = [name for name in column_list
                       if not name.endswith('_err')]

        # We get the transmission table and effective wavelength for each
        # used filter.
        transmission = {}
        effective_wavelength = {}
        base = Database()
        for name in filter_list:
            filt = base.get_filter(name)
            transmission[name] = filt.trans_table
            effective_wavelength[name] = filt.effective_wavelength
        base.close()

        # Read the observation table
        obs_table = atpy.Table(data_file)

        # 2D array containing the chi-squares. The first axis in the
        # observation number, the second axis is the model number (i.e. the
        # index of the parametre dictionary in sed_modules_params.
        chi_square_table = 99 * np.ones((obs_table.data.shape[0],
                                         len(sed_modules_params)))

        # Same 2D array for the normalisation factor.
        norm_factor_table = np.zeros((obs_table.data.shape[0],
                                      len(sed_modules_params)))

        # We complete the observation data by adding error where none is
        # provided and by adding the systematic deviation.
        for name in filter_list:
            name_err = name + '_err'
            if name_err not in column_list:
                if name_err not in obs_table.columns:
                    obs_table.add_column(name_err,
                                         np.zeros(obs_table.data.shape),
                                         dtype=float)
                else:
                    obs_table[name_err] = np.zeros(obs_table.data.shape)

            obs_table[name_err] = adjust_errors(obs_table[name],
                                                obs_table[name_err])

        # We loop over all the possible theoretical SEDs
        for model_index, parametres in enumerate(sed_modules_params):

            sed = create_sed(sed_modules, parametres)

            # Theoretical fluxes
            theor_fluxes = [sed.compute_fnu(transmission[name],
                                            effective_wavelength[name])
                            for name in filter_list]

            # Compute the reduced Chi-square and normalisation factor for each
            # observed SEDs
            for obs_index in range(obs_table.data.shape[0]):
                obs_fluxes = [obs_table[name][obs_index]
                              for name in filter_list]
                obs_errors = [obs_table[name + '_err'][obs_index]
                              for name in filter_list]

                chi2, norm_factor = compute_chi2(theor_fluxes,
                                                 obs_fluxes,
                                                 obs_errors)
                chi_square_table[obs_index, model_index] = chi2
                norm_factor_table[obs_index, model_index] = norm_factor

        # Find the model corresponding to the least reduced Chi-square for
        # each observation.
        results = []
        for obs_index in range(obs_table.data.shape[0]):
            # If there more than one model with the minimal chi-square value
            # only the first is returned.
            best_chi2 = min(chi_square_table[obs_index])
            best_index = chi_square_table[obs_index].argmin()
            best_norm_factor = norm_factor_table[obs_index, best_index]
            best_params = sed_modules_params[best_index]
            best_sed = create_sed(sed_modules, best_params)
            results.append((best_sed,
                            best_params,
                            best_chi2,
                            best_norm_factor))

            return results


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
