# -*- coding: utf-8 -*-
"""
Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

This file implements the statistical analysis as performed by the calcX2_psum
programme of the Fortran Cigale code.

The models corresponding to all possible combinations of parameters are
computed are the integrated flux in the same filters as the observations are
used to compute the χ² of the fitting. This χ² give a propability that is
associated with the model values for the parameters. At the end, for each
parameter, the (probability) weighted mean and standard deviation are computed
and the best fitting model (the one with the least reduced χ²) is given.

"""

import os
import sys
import atpy
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
from . import common
from ..sed.warehouse import create_sed
from ..data import Database


# Tolerance threshold under which any flux or error is considered as 0.
TOLERANCE = 1.e-12
# Name of the fits file containing the results
RESULT_FILE = 'psum_results.xml'
# Directory where the output files are storeds
OUT_DIR = 'out/'


class Module(common.AnalysisModule):
    """psum analysis

    TODO: Description of the PSUM method.
    """

    parameter_list = {
        "analysed_variables": (
            "array of strings",
            "List of the variables (in the SEDs info dictionaries) for which "
            "the statistical analysis will be done.",
            ["sfr", "average_sfr"]
        ),
        "save_best_sed": (
            "boolean",
            "If true, save the best SED for each observation to a file.",
            False
        ),
        "plot_best_sed": (
            "boolean",
            "If true, for each observation save a plot of the best SED "
            "and the observed fluxes.",
            False
        ),
        "plot_chi2_distribution": (
            "boolean",
            "If true, for each observation and each analysed variable "
            "plot the value vs reduced chi-square distribution.",
            False
        )
    }

    def process(self, data_file, column_list, sed_modules,
                sed_modules_params, parameters):
        """Process with the psum analysis.

        The analysis is done in two nested loops: over each observation and
        over each theoretical SEDs. We first loop over the SEDs to limit the
        number of time the SEDs are created.

        Parameters
        ----------
        data_file: string
            Name of the file containing the observations to fit.
        column_list: list of strings
            Name of the columns from the data file to use for the analysis.
        sed_modules: list of strings
            List of the module names (in the right order) to use for creating
            the SEDs.
        sed_modules_params: list of dictionaries
            List of the parameter dictionaries for each module.
        parameters: dictionary
            Dictionnary containing the parameters.

        Returns
        -------
        best_sed: list of tuples (pcigale.sed object, dict, float, float)
            There is one tuple per observed object: the first element is the
            best fitting SED for this object, the second dictionary of
            parameter used to produce it, the third is the reduced Chi-square
            of the fit and the fourth is the normalisation factor to be
            applied to the SED to fit the observation.
        results: dictionary
            Dictionary associating to NAME and NAME_err the weighted average
            and standard deviation lists, where name is galaxy_mass or the
            content of the analysed variables list. Each key is associated
            to an array which index corresponds to the rows in the observation
            data file. This dictionary is also saved as a FITS file to the
            disk.

        """

        # Create the output directory and stop it exists.
        try:
            os.mkdir(OUT_DIR)
        except OSError:
            print("pcigale can't create the {} directory, maybe "
                  "it yet exists.".format(OUT_DIR))
            sys.exit()

        # Get the parameters
        analysed_variables = parameters["analysed_variables"]
        save_best_sed = parameters["save_best_sed"]
        plot_best_sed = parameters["plot_best_sed"]
        plot_chi2_distribution = parameters["plot_chi2_distribution"]

        best_sed_list = []
        results = {'galaxy_mass': [], 'galaxy_mass_err': []}
        for variable in analysed_variables:
            results[variable] = []
            results[variable + '_err'] = []

        # We get the transmission table and effective wavelength for each
        # used filter.
        filter_list = [name for name in column_list
                       if not name.endswith('_err')]
        transmission = {}
        effective_wavelength = {}
        base = Database()
        for name in filter_list:
            filt = base.get_filter(name)
            transmission[name] = filt.trans_table
            effective_wavelength[name] = filt.effective_wavelength
        base.close()

        # Read the observation table and complete it by adding error where
        # none is provided and by adding the systematic deviation.
        obs_table = atpy.Table(data_file)
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

        # As we loop fist on the models (to limit the number of times a model
        # is computed) we need to store the computation results to make the
        # analysis for each observation in a second time. For this, we use a
        # three axis numpy array: axis 1 is the model (based on the index of
        # the sed_module_params list), axis 2 is the observation (base on the
        # data table row index) and axis 3 is the considered variable (based
        # on the analysed variables list + chi2min, probability and
        # galaxy_mass at the beginning).
        comp_table = np.zeros((len(sed_modules_params),
                               obs_table.data.shape[0],
                               len(analysed_variables) + 3), dtype=float)
        comp_table[:, :, :] = np.nan

        # We loop over all the possible theoretical SEDs
        for model_index, parameters in enumerate(sed_modules_params):
            sed = create_sed(sed_modules, parameters)

            # Compute the reduced Chi-square, the galaxy mass (normalisation
            # factor) and probability for each observed SEDs. Add these and
            # the values for the analysed variable to the comp_table.
            for obs_index in range(obs_table.data.shape[0]):
                obs_redshift = obs_table['redshift'][obs_index]
                obs_fluxes = [obs_table[name][obs_index]
                              for name in filter_list]
                obs_errors = [obs_table[name + '_err'][obs_index]
                              for name in filter_list]

                # Theoretical fluxes
                theor_fluxes = [sed.compute_fnu(transmission[name],
                                                effective_wavelength[name],
                                                obs_redshift)
                                for name in filter_list]

                chi2min, galaxy_mass, probability = compute_chi2(theor_fluxes,
                                                                 obs_fluxes,
                                                                 obs_errors)
                comp_table[model_index, obs_index, 0] = chi2min
                comp_table[model_index, obs_index, 1] = probability
                comp_table[model_index, obs_index, 2] = galaxy_mass

                for index, variable in enumerate(analysed_variables):
                    if variable in sed.mass_proportional_info:
                        comp_table[model_index, obs_index, index + 3] = \
                            galaxy_mass * sed.info[variable]
                    else:
                        comp_table[model_index, obs_index, index + 3] = \
                            sed.info[variable]

        # Find the model corresponding to the least reduced Chi-square for
        # each observation.
        # Now we loop over the observations.
        for obs_index, obs_name in enumerate(obs_table['id']):
            # Find the model corresponding to the least reduced Chi-square;
            # if there more than one model with the minimal chi-square value
            # only the first is returned.
            best_index = comp_table[:, obs_index, 0].argmin()
            best_chi2 = comp_table[best_index, obs_index, 0]
            best_norm_factor = comp_table[best_index, obs_index, 2]
            best_params = sed_modules_params[best_index]
            best_sed = create_sed(sed_modules, best_params)
            best_sed_list.append((best_sed,
                                  best_params,
                                  best_chi2,
                                  best_norm_factor))

            # Save best SED
            # TODO: For now, we only save the lambda vs fnu table. Once
            # we develop a way to serialise the SED, we should save the
            # complete SED object.
            if save_best_sed:
                best_sed_lambda_fnu = best_sed.lambda_fnu(
                    redshift=obs_table['redshift'][obs_index])
                best_sed_table = atpy.Table()
                best_sed_table.add_column("wavelength",
                                          best_sed_lambda_fnu[0],
                                          "nm")
                best_sed_table.add_column("fnu",
                                          best_norm_factor
                                          * best_sed_lambda_fnu[1],
                                          "mJy")
                best_sed_table.write(OUT_DIR + obs_name + 'bestSED.xml')

            # Plot the best SED
            if plot_best_sed:
                best_sed_lambda_fnu = best_sed.lambda_fnu(
                    redshift=obs_table['redshift'][obs_index])
                figure = plt.figure()
                ax = figure.add_subplot(111)
                ax.loglog(best_sed_lambda_fnu[0],
                          best_norm_factor * best_sed_lambda_fnu[1],
                          '-b')
                ax.loglog([effective_wavelength[name] for name in filter_list],
                          [obs_table[name][obs_index] for name in filter_list],
                          'or')
                ax.set_xlabel('Wavelength [nm]')
                ax.set_ylabel('Flux [mJy]')
                ax.set_title(obs_name +
                             ' best fitting SED - chi2min=' +
                             str(best_chi2))
                figure.savefig(OUT_DIR + obs_name + '_bestSED.pdf')

            # Compute the statistics for the desired variables.
            for index, variable in enumerate(['galaxy_mass'] +
                                             analysed_variables):
                # The 'variable' axis in comp_table as chi2 and probability
                # values at the beginning.
                idx = index + 2

                mean, sigma = w_mean_sigma(comp_table[:, obs_index, idx],
                                           comp_table[:, obs_index, 1])

                results[variable].append(mean)
                results[variable + '_err'].append(sigma)

                # We plot all the (value, chi2min) tuples.
                if plot_chi2_distribution:
                    figure = plt.figure()
                    ax = figure.add_subplot(111)
                    ax.plot(values,
                            comp_table[:, obs_index, 0],
                            '.')
                    ax.set_xlabel('value')
                    ax.set_ylabel('reduced chi square')
                    ax.set_title(variable)
                    figure.savefig(OUT_DIR +
                                   obs_name + '_' + variable + '_chi2plot.pdf')

        # Write the results to the fits file
        result_table = atpy.Table()
        for variable in (['galaxy_mass'] + analysed_variables):
            result_table.add_column(variable, results[variable])
            result_table.add_column(variable + '_err',
                                    results[variable + '_err'])
        result_table.write(OUT_DIR + RESULT_FILE)

        return best_sed_list


def adjust_errors(flux, error, default_error=0.1, systematic_deviation=0.1):
    """Adjust the errors replacing the 0 values by the default error and
    adding the systematic deviation.

    The systematic deviation change the error to:
    sqrt( error² + (flux * deviation)² )

    Parameters
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

    Parameters
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
    probability: float
        Probability associated with the chi-square and the considered number
        of degrees of freedom.

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
        probability = 0
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
            probability = 1
        else:
            normalisation_factor = (sum(obs_fluxes * model_fluxes) /
                                    sum(model_fluxes * model_fluxes))
            norm_model_fluxes = normalisation_factor * model_fluxes
            chi2 = sum(np.square((obs_fluxes - norm_model_fluxes) /
                                 obs_errors))
            reduced_chi2 = chi2 / degrees_of_freedom
            reduced_chi2 = min(reduced_chi2, 99)

            probability = stats.chi2.sf(chi2, degrees_of_freedom)

    return reduced_chi2, normalisation_factor, probability


def w_mean_sigma(values, weights):
    """Return the weighted average and standard deviation

    Parameters
    ----------
    values : list of floats
        List of values.
    weights : list of floats
        List of weights, must have the same shape as value list.

    Returns
    -------
    mean: float
        Weighted average.
    sigma: float
        Standard deviation.
    """
    mean = np.average(values, weights=weights)
    variance = np.dot(weights, (values - mean) ** 2) / np.sum(weights)

    return (mean, np.sqrt(variance))
