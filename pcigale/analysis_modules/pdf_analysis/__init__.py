# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013-2014 Yannick Roehlly
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

"""
Probability Density Function analysis module
============================================

This module builds the probability density functions (PDF) of the SED
parameters to compute their moments.

The models corresponding to all possible combinations of parameters are
computed and their fluxes in the same filters as the observations are
integrated. These fluxes are compared to the observed ones to compute the
χ² value of the fitting. This χ² give a probability that is associated with
the model values for the parameters.

At the end, for each parameter, the probability-weighted mean and standard
deviation are computed and the best fitting model (the one with the least
reduced χ²) is given for each observation.

"""

import os
import numpy as np
import matplotlib
from numpy import newaxis
from collections import OrderedDict
from datetime import datetime
from progressbar import ProgressBar
from matplotlib import pyplot as plt
from astropy.table import Table, Column
from ...utils import read_table
from .. import AnalysisModule, complete_obs_table
from ...creation_modules import get_module as get_creation_module
from .utils import gen_compute_fluxes_at_redshift, gen_pdf, gen_best_sed_fig
from ...warehouse import SedWarehouse
from ...data import Database

# Tolerance threshold under which any flux or error is considered as 0.
TOLERANCE = 1.e-12
# Probability threshold: models with a lower probability are excluded from
# the moments computation.
MIN_PROBABILITY = 1.e-20
# Name of the file containing the analysis results
RESULT_FILE = "analysis_results.fits"
# Name of the file containing the best models information
BEST_MODEL_FILE = "best_models.fits"
# Directory where the output files are stored
OUT_DIR = "out/"
# Wavelength limits (restframe) when plotting the best SED.
PLOT_L_MIN = 91
PLOT_L_MAX = 1e6
# Number of points in the PDF
PDF_NB_POINTS = 1000


class PdfAnalysis(AnalysisModule):
    """PDF analysis module"""

    parameter_list = OrderedDict([
        ("analysed_variables", (
            "array of strings",
            "List of the variables (in the SEDs info dictionaries) for which "
            "the statistical analysis will be done.",
            ["sfr", "average_sfr"]
        )),
        ("use_observation_redshift", (
            "boolean",
            "If true, the redshift of each observation will be taken from the "
            "redshift column in the observation table (default) and you must "
            "not use a redshifting module. Set to false if you want to use a "
            "redshifting module to test various redshifts.",
            True
        )),
        ("save_best_sed", (
            "boolean",
            "If true, save the best SED for each observation to a file.",
            False
        )),
        ("plot_best_sed", (
            "boolean",
            "If true, for each observation save a plot of the best SED "
            "and the observed fluxes.",
            False
        )),
        ("plot_chi2_distribution", (
            "boolean",
            "If true, for each observation and each analysed variable "
            "plot the value vs reduced chi-square distribution.",
            False
        )),
        ("save_pdf", (
            "boolean",
            "If true, for each observation and each analysed variable "
            "save the probability density function.",
            False
        )),
        ("plot_pdf", (
            "boolean",
            "If true, for each observation and each analysed variable "
            "plot the probability density function.",
            False
        )),
        ("storage_type", (
            "string",
            "Type of storage used to cache the generate SED.",
            "memory"
        ))
    ])

    def process(self, data_file, column_list, creation_modules,
                creation_modules_params, parameters):
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
        creation_modules: list of strings
            List of the module names (in the right order) to use for creating
            the SEDs.
        creation_modules_params: list of dictionaries
            List of the parameter dictionaries for each module.
        parameters: dictionary
            Dictionary containing the parameters.

        """

        # To be sure matplotlib will not display the interactive window.
        matplotlib.interactive(0)

        # Rename the output directory if it exists
        if os.path.exists(OUT_DIR):
            new_name = datetime.now().strftime("%Y%m%d%H%M") + "_" + OUT_DIR
            os.rename(OUT_DIR, new_name)
            print("The existing {} directory was renamed to {}".format(
                OUT_DIR,
                new_name
            ))
        os.mkdir(OUT_DIR)

        # Get the parameters
        analysed_variables = parameters["analysed_variables"]
        use_observation_redshift = (parameters["use_observation_redshift"]
                                    .lower() == "true")
        save_best_sed = (parameters["save_best_sed"].lower() == "true")
        plot_best_sed = (parameters["plot_best_sed"].lower() == "true")
        plot_chi2_distribution = (
            parameters["plot_chi2_distribution"].lower() == "true")
        save_pdf = (parameters["save_pdf"].lower() == "true")
        plot_pdf = (parameters["plot_pdf"].lower() == "true")

        # Get the needed filters in the pcigale database. We use an ordered
        # dictionary because we need the keys to always be returned in the
        # same order.
        with Database() as base:
            filters = OrderedDict([(name, base.get_filter(name))
                                   for name in column_list
                                   if not name.endswith('_err')])

        # If needed, get the redshifting and IGM attenuation modules (we set
        # redshift=0 but the redshift value is adapted later).
        redshifting_module = None
        igm_module = None
        if use_observation_redshift:
            redshifting_module = get_creation_module("redshifting", redshift=0)
            igm_module = get_creation_module("igm_attenuation")

        # Read the observation table and complete it by adding error where
        # none is provided and by adding the systematic deviation.
        obs_table = complete_obs_table(
            read_table(data_file),
            column_list,
            filters,
            TOLERANCE
        )

        ##################################################################
        # Model computation                                              #
        ##################################################################

        print("Computing the models fluxes...")

        # First, we compute for all the possible theoretical models (one for
        # each parameter set in sed_module_parameters) the fluxes in all the
        # filters at the redshift of each observed galaxy.  These fluxes are
        # stored in:

        # model_fluxes:
        # - axis 0: model index
        # - axis 1: observation index
        # - axis 2: filter index

        # We use a numpy masked array to mask the fluxes of models that would
        # be older than the age of the Universe at the considered redshift

        # The values for the analysed variables are stored in:

        # model_variables:
        # - axis 0: the model index in sed_module_params
        # - axis 1: the variable index in analysed_variables

        model_fluxes = np.ma.zeros((len(creation_modules_params),
                                    len(obs_table),
                                    len(filters)))
        model_variables = np.ma.zeros((len(creation_modules_params),
                                       len(analysed_variables)))

        # We keep the information (i.e. the content of the sed.info
        # dictionary) for each model.
        model_info = []

        progress_bar = ProgressBar(maxval=len(creation_modules_params)).start()

        # The SED warehouse is used to retrieve SED corresponding to some
        # modules and parameters.
        with SedWarehouse(cache_type=parameters["storage_type"]) as \
                sed_warehouse:

            for model_index, model_params in enumerate(
                    creation_modules_params):

                sed = sed_warehouse.get_sed(creation_modules, model_params)

                # Cached function to compute the SED fluxes at a redshift
                gen_fluxes = gen_compute_fluxes_at_redshift(
                    sed, filters.values(), redshifting_module, igm_module)

                for obs_index, redshift in enumerate(obs_table["redshift"]):
                    model_fluxes[model_index, obs_index] = gen_fluxes(redshift)

                model_variables[model_index] = np.array(
                    [sed.info[name] for name in analysed_variables]
                )

                model_info.append(sed.info.values())

                progress_bar.update(model_index + 1)

        # Mask the invalid fluxes
        model_fluxes = np.ma.masked_less(model_fluxes, -90)

        progress_bar.finish()

        ##################################################################
        # Observations to models comparison                              #
        ##################################################################

        print("Comparing the observations to the models...")

        # To accelerate the computations, we put observations fluxes and
        # errors in multi-dimensional numpy array:
        # - axis 0: the observation index
        # - axis 1: the filter index
        obs_fluxes = np.array([
            obs_table[name] for name in filters]
        ).transpose()
        obs_errors = np.array([
            obs_table[name + "_err"] for name in filters]
        ).transpose()

        # Some observations may not have flux value in some filters, in that
        # case the user is asked to put -9999 as value. We mask these values.
        # Note, we must mask obs_fluxes after obs_errors.
        obs_errors = np.ma.masked_where(obs_fluxes < -9990., obs_errors)
        obs_fluxes = np.ma.masked_less(obs_fluxes, -9990.)

        # Normalisation factor to be applied to a model fluxes to best fit an
        # observation fluxes. Normalised flux of the models. χ² and
        # likelihood of the fitting. Reduced χ² (divided by the number of
        # filters to do the fit).
        # axis 0: model index
        # axis 1: observation index
        normalisation_factors = (
            np.sum(
                model_fluxes * obs_fluxes[newaxis, :, :] / (
                    obs_errors * obs_errors)[newaxis, :, :], axis=2
            ) / np.sum(
                model_fluxes * model_fluxes / (
                    obs_errors * obs_errors)[newaxis, :, :], axis=2)
        )
        norm_model_fluxes = model_fluxes * normalisation_factors[:, :, newaxis]

        # χ² of the comparison of each model to each observation. The
        # chi_squares array is:
        # axis 0: model index
        # axis 1: observation index
        chi_squares = np.sum(
            np.square(
                (obs_fluxes[newaxis, :, :] - norm_model_fluxes) /
                obs_errors[newaxis, :, :]
            ), axis=2
        )

        # We define the reduced χ² as the χ² divided by the number of fluxes
        # used for the fitting.
        reduced_chi_squares = chi_squares / obs_fluxes.count(1)[newaxis, :]

        # We use the exponential probability associated with the χ² as
        # likelihood function. The likelihood array is:
        # axis 0: model index
        # axis 1: observation index
        likelihood = np.exp(-chi_squares/2)
        # For the analysis, we consider that the computed models explain each
        # observation. We normalise the likelihood function to have a
        # total likelihood of 1 for each observation.
        likelihood /= np.sum(likelihood, axis=0)[newaxis, :]
        # We don't want to take into account the models with a probability
        # less that the threshold.
        likelihood = np.ma.masked_less(likelihood, MIN_PROBABILITY)
        # We re-normalise the likelihood.
        likelihood /= np.sum(likelihood, axis=0)[newaxis, :]

        # Some model variables depend on the galaxy mass (the normalisation
        # factor) so we expand the model_variables array to have a new axis
        # corresponding to the observation and multiply (when needed) the
        # value of the variables by the normalisation factor of the fitting
        # for each observation.
        # The new array will have be:
        # axis 0: model index
        # axis 1: observation index
        # axis 2: variable index
        model_variables = model_variables[:, newaxis, :].repeat(len(obs_table),
                                                                axis=1)
        # We take the mass-dependent variable list from the last computed sed.
        for index, variable in enumerate(analysed_variables):
            if variable in sed.mass_proportional_info:
                model_variables[:, :, index] *= normalisation_factors

        # We also add the galaxy mass to the analysed variables
        if sed.sfh is not None:
            analysed_variables.insert(0, "galaxy_mass")
            model_variables = np.dstack((normalisation_factors, model_variables))

        ##################################################################
        # Variable analysis                                              #
        ##################################################################

        print("Analysing the variables...")

        # We compute the weighted average and standard deviation using the
        # likelihood as weight. We first build the weight array by expanding
        # the likelihood along a new axis corresponding to the analysed
        # variable.
        weights = likelihood[:, :, newaxis].repeat(len(analysed_variables),
                                                   axis=2)

        # Analysed variables average and standard devisation arrays.
        # axis 0: observation index
        # axis 1: variable index
        analysed_averages = np.ma.average(model_variables,
                                          axis=0,
                                          weights=weights)

        analysed_variances = np.ma.average(
            (model_variables - analysed_averages[newaxis, :, :])**2,
            axis=0,
            weights=weights)

        analysed_std = np.ma.sqrt(analysed_variances)

        # Create and save the result table.
        result_table = Table()
        result_table.add_column(Column(
            obs_table["id"].data,
            name="observation_id"
        ))
        for index, variable in enumerate(analysed_variables):
            result_table.add_column(Column(
                analysed_averages[:, index],
                name=variable
            ))
            result_table.add_column(Column(
                analysed_std[:, index],
                name=variable+"_err"
            ))
        result_table.write(OUT_DIR + RESULT_FILE)

        ##################################################################
        # Best models                                                    #
        ##################################################################

        print("Analysing the best models...")

        # We define the best fitting model for each observation as the one
        # with the least χ².
        best_model_index = list(chi_squares.argmin(axis=0))

        # We take the list of information added to the SEDs from the last
        # computed one.
        model_info_names = sed.info.keys()

        best_model_table = Table()
        best_model_table.add_column(Column(
            obs_table["id"].data,
            name="observation_id"
        ))
        best_model_table.add_column(Column(
            chi_squares[best_model_index, range(len(best_model_index))],
            name="chi_square"
        ))
        best_model_table.add_column(Column(
            reduced_chi_squares[best_model_index, range(len(best_model_index))],
            name="reduced_chi_square"
        ))
        if sed.sfh is not None:
            best_model_table.add_column(Column(
                normalisation_factors[best_model_index,
                                    range(len(best_model_index))],
                name="galaxy_mass",
                unit="Msun"
            ))

        for index, name in enumerate(model_info_names):
            column = Column([list(model_info[model_idx])[index] for model_idx
                             in best_model_index], name=name)
            if name in sed.mass_proportional_info:
                column *= (normalisation_factors[best_model_index,
                                                range(len(best_model_index))])
            best_model_table.add_column(column)

        best_model_table.write(OUT_DIR + BEST_MODEL_FILE)

        if plot_best_sed or save_best_sed:

            print("Plotting/saving the best models...")

            with SedWarehouse(cache_type=parameters["storage_type"]) as \
                    sed_warehouse:
                for obs_index, obs_name in enumerate(obs_table["id"]):

                    obs_redshift = obs_table["redshift"][obs_index]
                    best_index = best_model_index[obs_index]

                    sed = sed_warehouse.get_sed(
                        creation_modules,
                        creation_modules_params[best_index]
                    )
                    if use_observation_redshift:
                        redshifting_module.parameters["redshift"] = \
                            obs_redshift
                        redshifting_module.process(sed)
                        igm_module.process(sed)

                    best_lambda = sed.wavelength_grid
                    best_fnu = sed.fnu * normalisation_factors[best_index,
                                                               obs_index]

                    if save_best_sed:
                        sed.to_votable(
                            OUT_DIR + "{}_best_model.xml".format(obs_name),
                            mass=normalisation_factors[best_index, obs_index]
                        )

                    if plot_best_sed:

                        plot_mask = (
                            (best_lambda >= PLOT_L_MIN * (1 + obs_redshift)) &
                            (best_lambda <= PLOT_L_MAX * (1 + obs_redshift))
                        )

                        figure = gen_best_sed_fig(
                            best_lambda[plot_mask],
                            best_fnu[plot_mask],
                            [f.effective_wavelength for f in filters.values()],
                            norm_model_fluxes[best_index, obs_index, :],
                            [obs_table[f][obs_index] for f in filters]
                        )

                        if figure is None:
                            print("Can not plot best model for observation "
                                  "{}!".format(obs_name))
                        else:
                            figure.suptitle(
                                u"Best model for {} - red-chi² = {}".format(
                                    obs_name,
                                    reduced_chi_squares[best_index, obs_index]
                                )
                            )
                            figure.savefig(OUT_DIR + "{}_best_model.pdf".format(
                                obs_name))
                            plt.close(figure)

        ##################################################################
        # Probability Density Functions                                  #
        ##################################################################

        # We estimate the probability density functions (PDF) of the
        # parameters using a weighted kernel density estimation. This part
        # should definitely be improved as we simulate the weigth by adding
        # as many value as their probability * 100.

        if save_pdf or plot_pdf:

            print("Computing the probability density functions...")

            for obs_index, obs_name in enumerate(obs_table["id"]):

                probabilities = likelihood[:, obs_index]

                for var_index, var_name in enumerate(analysed_variables):

                    values = model_variables[:, obs_index, var_index]

                    pdf_grid = np.linspace(values.min(), values.max(),
                                           PDF_NB_POINTS)
                    pdf_prob = gen_pdf(values, probabilities, pdf_grid)

                    if pdf_prob is None:
                        # TODO: use logging
                        print("Can not compute PDF for observation <{}> and "
                              "variable <{}>.".format(obs_name, var_name))

                    if save_pdf and pdf_prob is not None:
                        table = Table((
                            Column(pdf_grid, name=var_name),
                            Column(pdf_prob, name="probability density")
                        ))
                        table.write(OUT_DIR + "{}_{}_pdf.fits".format(
                            obs_name, var_name))

                    if plot_pdf and pdf_prob is not None:
                        figure = plt.figure()
                        ax = figure.add_subplot(111)
                        ax.plot(pdf_grid, pdf_prob)
                        ax.set_xlabel(var_name)
                        ax.set_ylabel("Probability density")
                        figure.savefig(OUT_DIR + "{}_{}_pdf.pdf".format(
                            obs_name, var_name))
                        plt.close(figure)

        ##################################################################
        # Reduced-chisquares plots                                       #
        ##################################################################
        if plot_chi2_distribution:

            print("Plotting the reduced chi squares distributions...")

            for obs_index, obs_name in enumerate(obs_table["id"]):

                obs_red_chisquares = reduced_chi_squares[:, obs_index]

                for var_index, var_name in enumerate(analysed_variables):

                    values = model_variables[:, obs_index, var_index]

                    figure = plt.figure()
                    ax = figure.add_subplot(111)
                    ax.plot(values, obs_red_chisquares, "ob")
                    ax.set_xlabel(var_name)
                    ax.set_ylabel("reduced chi-square")
                    figure.suptitle("Reduced chi-square distribution of {} "
                                    "values for {}".format(obs_index, var_name))
                    figure.savefig(OUT_DIR + "{}_{}_chisquares.pdf".format(
                            obs_name, var_name))
                    plt.close(figure)


# AnalysisModule to be returned by get_module
Module = PdfAnalysis
