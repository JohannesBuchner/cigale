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
from copy import deepcopy
from scipy import stats
from progressbar import ProgressBar
from matplotlib import pyplot as plt
from . import common
from ..warehouse import SedWarehouse
from ..sed.modules.common import get_module
from ..data import Database


# Tolerance threshold under which any flux or error is considered as 0.
TOLERANCE = 1.e-12
# Name of the fits file containing the results
RESULT_FILE = 'psum_results.fits'
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
        ),
        "save_pdf": (
            "boolean",
            "If true, for each observation and each analysed variable "
            "save the probability density function.",
            False
        ),
        "plot_pdf": (
            "boolean",
            "If true, for each observation and each analysed variable "
            "plot the probability density function.",
            False
        ),
        "pdf_max_bin_number": (
            "integer",
            "Maximum number of bins used to compute the probability density "
            "function. This is only used when saving or printing the PDF. "
            "If there are less values, the probability is given for each "
            "one.",
            50
        ),
        "storage_type": (
            "string",
            "Type of storage used to cache the generate SED.",
            "memory"
        )
    }

    def process(self, data_file, column_list, sed_modules,
                sed_modules_params, redshift_module_name,
                redshift_configuration, parameters):
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
        redshift_module_name : string
            Name of the module used to redshift the SED.
        redshift_configuration : dictionary
            Configuration dictionary for the module used to redshift the SED.
        parameters: dictionary
            Dictionnary containing the parameters.

        """

        # Create the output directory and stop it exists.
        try:
            os.mkdir(OUT_DIR)
        except OSError:
            print("pcigale can't create the {} directory, maybe "
                  "it yet exists.".format(OUT_DIR))
            sys.exit()

        # Open the warehouse
        # TODO Why is parameters["storage_type"] an array?
        sed_warehouse = SedWarehouse(
            cache_type=parameters["storage_type"][0])

        # Get the parameters
        analysed_variables = parameters["analysed_variables"]
        save_best_sed = parameters["save_best_sed"]
        plot_best_sed = parameters["plot_best_sed"]
        plot_chi2_distribution = parameters["plot_chi2_distribution"]
        save_pdf = parameters["save_pdf"]
        plot_pdf = parameters["plot_pdf"]
        pdf_max_bin_number = parameters["pdf_max_bin_number"]

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

        # We get the redshift module.
        redshift_module = get_module(redshift_module_name)
        redshift_module.parameters = redshift_configuration

        # Read the observation table and complete it by adding error where
        # none is provided and by adding the systematic deviation.
        obs_table = atpy.Table(data_file, verbose=False)
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
        # on the analysed variables list + reduced_chi2, probability and
        # galaxy_mass at the beginning).
        comp_table = np.zeros((len(sed_modules_params),
                               obs_table.data.shape[0],
                               len(analysed_variables) + 3), dtype=float)
        comp_table[:, :, :] = np.nan

        # We loop over all the possible theoretical SEDs
        progress_bar = ProgressBar(maxval=len(sed_modules_params)).start()
        for model_index, parameters in enumerate(sed_modules_params):
            sed = sed_warehouse.get_sed(sed_modules, parameters)

            # Compute the reduced Chi-square, the galaxy mass (normalisation
            # factor) and probability for each observed SEDs. Add these and
            # the values for the analysed variable to the comp_table.
            for obs_index in range(obs_table.data.shape[0]):
                obs_redshift = obs_table['redshift'][obs_index]
                obs_fluxes = [obs_table[name][obs_index]
                              for name in filter_list]
                obs_errors = [obs_table[name + '_err'][obs_index]
                              for name in filter_list]

                # We copy the SED before redshifting it.
                red_sed = deepcopy(sed)
                redshift_module.parameters["redshift"] = obs_redshift
                redshift_module.process(red_sed)

                # Theoretical fluxes
                theor_fluxes = [red_sed.compute_fnu(transmission[name],
                                                    effective_wavelength[name],
                                                    obs_redshift)
                                for name in filter_list]

                reduced_chi2, galaxy_mass, probability = compute_chi2(
                    theor_fluxes, obs_fluxes, obs_errors)
                comp_table[model_index, obs_index, 0] = reduced_chi2
                comp_table[model_index, obs_index, 1] = probability
                comp_table[model_index, obs_index, 2] = galaxy_mass

                for index, variable in enumerate(analysed_variables):
                    if variable in sed.mass_proportional_info:
                        comp_table[model_index, obs_index, index + 3] = \
                            galaxy_mass * sed.info[variable]
                    else:
                        comp_table[model_index, obs_index, index + 3] = \
                            sed.info[variable]

            progress_bar.update(model_index + 1)

        progress_bar.finish()

        # Loop over the observations to find the best fitting model and
        # compute the parametre statistics.
        for obs_index, obs_name in enumerate(obs_table['id']):
            # We will save the part of the computation table corresponding to
            # the model as a FITS file.
            fits_table = atpy.Table()
            fits_table.table_name = "Analysis computation table"
            fits_table.add_column("reduced_chi-square",
                                  comp_table[:, obs_index, 0])
            fits_table.add_column("probability",
                                  comp_table[:, obs_index, 1])

            # Find the model corresponding to the least reduced Chi-square;
            # if there more than one model with the minimal chi-square value
            # only the first is returned.
            best_index = comp_table[:, obs_index, 0].argmin()
            best_chi2 = comp_table[best_index, obs_index, 0]
            best_norm_factor = comp_table[best_index, obs_index, 2]
            best_params = sed_modules_params[best_index]
            best_sed = sed_warehouse.get_sed(sed_modules, best_params)

            # Save best SED
            # TODO: For now, we only save the lambda vs fnu table. Once
            # we develop a way to serialise the SED, we should save the
            # complete SED object.
            if save_best_sed:
                best_sed_lambda_fnu = best_sed.lambda_fnu(
                    redshift=obs_table['redshift'][obs_index])
                best_sed_table = atpy.Table()
                best_sed_table.table_name = "Best SED"
                best_sed_table.add_column("wavelength",
                                          best_sed_lambda_fnu[0],
                                          "nm")
                best_sed_table.add_column("fnu",
                                          best_norm_factor
                                          * best_sed_lambda_fnu[1],
                                          "mJy")
                best_sed_table.write(OUT_DIR + obs_name + 'bestSED.fits',
                                     verbose=False)

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
                             ' best fitting SED - reduced chi2:' +
                             str(best_chi2))
                figure.savefig(OUT_DIR + obs_name + '_bestSED.pdf')

            # Compute the statistics for the desired variables. First, we
            # build the probability distribution function (PDF) for the
            # variable, then we use it to compute the expected value and
            # the standard deviation.
            for index, variable in enumerate(['galaxy_mass'] +
                                             analysed_variables):
                # The variable index in the last axis of comp_table is
                # index+2 as chi2 and probability are at the beginning.
                # values at the beginning.
                idx = index + 2

                values = comp_table[:, obs_index, idx]
                probabilities = comp_table[:, obs_index, 1]

                fits_table.add_column(variable, values)

                # We compute the Probability Density Function by binning the
                # data into evenly populated bins. To estimate the binning
                # quality, we also provide the standard deviation of the
                # variable and probability values in the bin.
                pdf_values = []
                pdf_probs = []
                pdf_value_sigma = []
                pdf_prob_sigma = []

                # The maximum number of bins is the least between
                # pdf_max_bin_number and the number of distinct values
                # for the analyse variable.
                pdf_bin_number = min(pdf_max_bin_number,
                                     len(np.unique(values)))

                pdf_bin_boundaries, pdf_bins = bin_evenly(values,
                                                          pdf_bin_number)
                pdf_bin_sizes = [
                    pdf_bin_boundaries[i + 1] - pdf_bin_boundaries[i]
                    for i in range(pdf_bin_number)
                ]

                for bin in range(1, pdf_bin_number + 1):
                    # The bin probability is the average of the probabilities
                    # in the bin.
                    pdf_probs.append(
                        np.average(probabilities[pdf_bins == bin]))
                    pdf_prob_sigma.append(
                        np.std(probabilities[pdf_bins == bin])
                    )
                    pdf_values.append(
                        (pdf_bin_boundaries[bin] +
                         pdf_bin_boundaries[bin - 1]) / 2
                    )
                    pdf_value_sigma.append(
                        np.std(values[pdf_bins == bin])
                    )

                pdf_probs = np.array(pdf_probs)
                pdf_prob_sigma = np.array(pdf_prob_sigma)
                pdf_values = np.array(pdf_values)
                pdf_value_sigma = np.array(pdf_value_sigma)

                # We normalise the PDF to 1 over all the parameter values.
                pdf_probs = pdf_probs / np.sum(pdf_bin_sizes * pdf_probs)

                # From the PDF, compute the expected value and standard
                # deviation.
                expected_value = np.sum(pdf_values * pdf_bin_sizes * pdf_probs)

                sigma = np.sqrt(
                    np.sum((pdf_values - expected_value) ** 2 *
                           pdf_bin_sizes *
                           pdf_probs)
                )
                results[variable].append(expected_value)
                results[variable + '_err'].append(sigma)

                # We plot all the (value, reduced_chi2) tuples.
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

                if save_pdf:
                    pdf_table = atpy.Table()
                    pdf_table.table_name = "Probability Density Function"
                    pdf_table.add_column("bin_start",
                                         pdf_bin_boundaries[:-1])
                    pdf_table.add_column("bin_end",
                                         pdf_bin_boundaries[1:])
                    pdf_table.add_column("value", pdf_values)
                    pdf_table.add_column("value_sigma", pdf_value_sigma)
                    pdf_table.add_column("probability", pdf_probs)
                    pdf_table.add_column("probability_sigma", pdf_prob_sigma)
                    pdf_table.write(OUT_DIR + obs_name + "_" + variable +
                                    "_pdf.fits", verbose=False)

                if plot_pdf:
                    figure = plt.figure()
                    ax = figure.add_subplot(111)
                    ax.bar(
                        pdf_bin_boundaries[:-1],
                        pdf_probs,
                        pdf_bin_sizes
                    )
                    ax.axvline(expected_value, color="r")
                    ax.axvline(expected_value - sigma,
                               color="r",
                               linestyle="-.")
                    ax.axvline(expected_value + sigma,
                               color="r",
                               linestyle="-.")
                    ax.set_title(obs_name + ' ' + variable + ' PDF')
                    ax.set_xlabel(variable)
                    ax.set_ylabel('Probability')
                    figure.savefig(OUT_DIR + obs_name + "_" + variable +
                                   "_pdf.pdf")

        # Write the computation table FITS
        fits_table.write(OUT_DIR + obs_name + '_comptable.fits', verbose=False)

        # Write the results to the fits file
        result_table = atpy.Table()
        result_table.table_name = "Analysis results"
        result_table.add_column('id', obs_table['id'])
        for variable in (['galaxy_mass'] + analysed_variables):
            result_table.add_column(variable, results[variable])
            result_table.add_column(variable + '_err',
                                    results[variable + '_err'])
        result_table.write(OUT_DIR + RESULT_FILE, verbose=False)

        sed_warehouse.close()


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
        Probability associated with the chi-square.

    """

    # The three arrays must have the same length.
    if (len(model_fluxes) != len(obs_fluxes) or
            len(obs_fluxes) != len(obs_errors)):
        raise ValueError("The model fluxes, observation fluxes and "
                         "observation errors arrays must have the "
                         "same length.")

    # We copy the arrays not to modify the original ones.
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

        # We consider that we fit the observation to the SED, independently of
        # the number of parameters used to build it. So the number of degrees
        # of freedom depends only on the number of fluxes.
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

            # We use the probability from the chi-square law for the
            # considered degrees of freedom.
            probability = stats.chi2.sf(chi2, degrees_of_freedom)

    return reduced_chi2, normalisation_factor, probability


def bin_evenly(values, max_bins):
    """Divide some values in evenly populated bins

    Given a list of values and a desired number of bins, this method computes
    the bins boundaries to have bins with the same number of elements in each
    one and digitises the value list with these boundaries.

    Parameters
    ----------
    values : list of floats
        List of values to be binned.
    max_bins : integer
        Maximum number of bins. If there are less distinct value, every value
        is in it's own bin.

    Returns
    -------
    boundaries : array of floats
        The value of the boundaries of the bins.
    bins_digits : numpy array of integers
        Array of the same length as the value list giving for each value the
        bin number (between 1 and nb_of_bins) it belongs to.

    """
    # If there are less values than asked bins, raise an error.
    if max_bins > len(values):
        max_bins = len(values)

    # The bin boundaries are the nb_of_bins + 1 quantiles.
    quantiles = np.linspace(0, 1, max_bins + 1)
    boundaries = stats.mstats.mquantiles(values, quantiles)

    # Because of the way np.digitize works, we must have the last boundary
    # higher than the value maximum to have this maximum belong to the last
    # bin.
    digitize_boundaries = np.copy(boundaries)
    digitize_boundaries[-1] += 1
    bin_digits = np.digitize(values, digitize_boundaries)

    return (boundaries, bin_digits)
