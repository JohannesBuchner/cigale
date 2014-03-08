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
from numpy import newaxis
from collections import OrderedDict
from datetime import datetime
from progressbar import ProgressBar
from astropy.table import Table, Column
from ...utils import read_table
from .. import AnalysisModule, complete_obs_table
from .utils import gen_pdf
from ...warehouse import SedWarehouse
from ...data import Database

# Tolerance threshold under which any flux or error is considered as 0.
TOLERANCE = 1e-12
# Probability threshold: models with a lower probability are excluded from
# the moments computation.
MIN_PROBABILITY = 1e-20
# Limit the redshift to this number of decimals
REDSHIFT_DECIMALS = 2
# Name of the file containing the analysis results
RESULT_FILE = "analysis_results.fits"
# Name of the file containing the best models information
BEST_MODEL_FILE = "best_models.fits"
# Directory where the output files are stored
OUT_DIR = "out/"
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
        ("save_best_sed", (
            "boolean",
            "If true, save the best SED for each observation to a file.",
            False
        )),
        ("save_chi2", (
            "boolean",
            "If true, for each observation and each analysed variable save "
            "the reduced chi².",
            False
        )),
        ("save_pdf", (
            "boolean",
            "If true, for each observation and each analysed variable save "
            "the probability density function.",
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
        save_best_sed = (parameters["save_best_sed"].lower() == "true")
        save_chi2 = (parameters["save_chi2"].lower() == "true")
        save_pdf = (parameters["save_pdf"].lower() == "true")

        # Get the needed filters in the pcigale database. We use an ordered
        # dictionary because we need the keys to always be returned in the
        # same order.
        with Database() as base:
            filters = OrderedDict([(name, base.get_filter(name))
                                   for name in column_list
                                   if not name.endswith('_err')])

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
        # filters. These fluxes are stored in:

        # model_fluxes:
        # - axis 0: model index
        # - axis 1: filter index

        # We use a numpy masked array to mask the fluxes of models that would
        # be older than the age of the Universe at the considered redshift.

        # The values for the analysed variables are stored in:

        # model_variables:
        # - axis 0: the model index in sed_module_params
        # - axis 1: the variable index in analysed_variables

        # For convenience, the redshift of each model is stored in
        # model_redshift.

        model_fluxes = np.ma.empty((len(creation_modules_params),
                                    len(filters)))
        model_variables = np.ma.empty((len(creation_modules_params),
                                       len(analysed_variables)))

        model_redshift = np.empty(len(creation_modules_params))

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

                model_fluxes[model_index, :] = np.array(
                    [sed.compute_fnu(filter_.trans_table,
                                     filter_.effective_wavelength)
                     for filter_ in filters.values()])
                model_variables[model_index, :] = np.array(
                    [sed.info[name] for name in analysed_variables]
                )

                model_redshift[model_index] = sed.info['redshift']

                model_info.append(sed.info.values())

                progress_bar.update(model_index + 1)

        unique_redshifts = np.unique(model_redshift)

        # Mask the invalid fluxes
        model_fluxes = np.ma.masked_less(model_fluxes, -90)

        progress_bar.finish()

        ##################################################################
        # Observations to models comparison                              #
        ##################################################################

        print("Comparing the observations to the models...")

        # As we are looping over all the observations we store data for the
        # output tables in various arrays
        analysed_averages_all = np.empty((len(obs_table),
                                          len(analysed_variables)))
        analysed_std_all = np.empty_like(analysed_averages_all)

        best_idx_all = np.empty(len(obs_table))
        best_chi2_all = np.empty_like(best_idx_all)
        best_chi2_red_all = np.empty_like(best_idx_all)
        normalisation_factors_all = np.empty_like(best_idx_all)
        
        best_fluxes = np.empty((len(obs_table), len(filters)))

        best_variables_all = [None]*len(obs_table)

        for idx_obs, obs in enumerate(obs_table):
            obs_fluxes = np.array([obs[name] for name in filters])
            obs_errors = np.array([obs[name + "_err"] for name in filters])

            # Some observations may not have flux value in some filters, in
            # that case the user is asked to put -9999 as value. We mask these
            # values. Note, we must mask obs_fluxes after obs_errors.
            obs_errors = np.ma.masked_where(obs_fluxes < -9990., obs_errors)
            obs_fluxes = np.ma.masked_less(obs_fluxes, -9990.)

            # We compute the χ² only for models with the closest redshift. We
            # extract model fluxes and information into arrays dedicated to a
            # given observation.
            closest_redshift = unique_redshifts[np.abs(obs["redshift"] -
                                                unique_redshifts).argmin()]
            w_models = model_redshift == closest_redshift
            model_fluxes_obs = model_fluxes[w_models, :]
            model_info_obs = np.array(model_info)[w_models]
            model_variables_obs = model_variables[w_models]

            # Normalisation factor to be applied to a model fluxes to best fit
            # an observation fluxes. Normalised flux of the models. χ² and
            # likelihood of the fitting. Reduced χ² (divided by the number of
            # filters to do the fit).
            normalisation_factors = (
                np.sum(
                    model_fluxes_obs * obs_fluxes / (
                        obs_errors * obs_errors), axis=1
                ) / np.sum(
                    model_fluxes_obs * model_fluxes_obs / (
                        obs_errors * obs_errors), axis=1)
            )
            norm_model_fluxes = (model_fluxes_obs *
                                 normalisation_factors[:, np.newaxis])

            # χ² of the comparison of each model to each observation.
            chi_squares = np.sum(
                np.square((obs_fluxes - norm_model_fluxes) / obs_errors),
                axis=1)

            # We define the reduced χ² as the χ² divided by the number of
            # fluxes used for the fitting.
            reduced_chi_squares = chi_squares / obs_fluxes.count()

            # We use the exponential probability associated with the χ² as
            # likelihood function.
            likelihood = np.exp(-chi_squares/2)
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
                    model_variables_obs[:, index] *= normalisation_factors

            # We also add the galaxy mass to the analysed variables if relevant
            if sed.sfh is not None:
                analysed_variables.insert(0, "galaxy_mass")
                model_variables_obs = np.dstack((normalisation_factors,
                                                 model_variables_obs))

            ##################################################################
            # Variable analysis                                              #
            ##################################################################

            print("Analysing the variables...")

            # We compute the weighted average and standard deviation using the
            # likelihood as weight. We first build the weight array by
            # expanding the likelihood along a new axis corresponding to the
            # analysed variable.
            weights = likelihood[:, newaxis].repeat(len(analysed_variables),
                                                    axis=1)

            # Analysed variables average and standard deviation arrays.
            analysed_averages = np.ma.average(model_variables_obs,
                                              axis=0, weights=weights)

            analysed_std = np.ma.sqrt(np.ma.average(
                (model_variables_obs - analysed_averages[newaxis, :])**2,
                axis=0, weights=weights))

            # We record the estimated averages and standard deviations to
            # save in a table later on when this has been computed for all
            # objects.
            analysed_averages_all[idx_obs, :] = analysed_averages
            analysed_std_all[idx_obs, :] = analysed_std

            ##################################################################
            # Best models                                                    #
            ##################################################################

            print("Analysing the best models...")

            # We define the best fitting model for each observation as the one
            # with the least χ².
            best_index = chi_squares.argmin()

            # We save the relevant data related to the model with the lowest
            # χ²
            best_idx_all[idx_obs] = best_index
            normalisation_factors_all[idx_obs] = \
                normalisation_factors[best_index]
            best_chi2_all[idx_obs] = chi_squares[best_index]
            best_chi2_red_all[idx_obs] = reduced_chi_squares[best_index]
            best_variables_all[idx_obs] = list(model_info_obs[best_index])
            best_fluxes[idx_obs, :] = model_fluxes_obs[best_index, :]

            if save_best_sed:

                print("Saving the best models...")

                with SedWarehouse(cache_type=parameters["storage_type"]) as \
                        sed_warehouse:

                    sed = sed_warehouse.get_sed(
                        creation_modules,
                        np.array(creation_modules_params)[w_models][best_index]
                    )

                    sed.to_votable(
                        OUT_DIR + "{}_best_model.xml".format(obs['id']),
                        mass=normalisation_factors[best_index]
                    )

            ##################################################################
            # Probability Density Functions                                  #
            ##################################################################

            # We estimate the probability density functions (PDF) of the
            # parameters using a weighted kernel density estimation. This part
            # should definitely be improved as we simulate the weight by adding
            # as many value as their probability * 100.
            if save_pdf:

                print("Computing the probability density functions...")

                for var_index, var_name in enumerate(analysed_variables):

                    values = model_variables_obs[:, var_index]

                    pdf_grid = np.linspace(values.min(), values.max(),
                                           PDF_NB_POINTS)
                    pdf_prob = gen_pdf(values, likelihood, pdf_grid)

                    if pdf_prob is None:
                        # TODO: use logging
                        print("Can not compute PDF for observation <{}> and "
                              "variable <{}>.".format(obs['id'], var_name))
                    else:
                        table = Table((
                            Column(pdf_grid, name=var_name),
                            Column(pdf_prob, name="probability density")
                        ))
                        table.write(OUT_DIR + "{}_{}_pdf.fits".format(
                            obs['id'], var_name))

            if save_chi2:

                print("Saving the chi²...")

                for var_index, var_name in enumerate(analysed_variables):
                    table = Table((
                        Column(model_variables_obs[:, var_index],
                               name=var_name),
                        Column(reduced_chi_squares, name="chi2")))
                    table.write(OUT_DIR + "{}_{}_chi2.fits".format(obs['id'],
                                var_name))

        # Create and save the result table.
        result_table = Table()
        result_table.add_column(Column(
            obs_table["id"].data,
            name="observation_id"
        ))
        for index, variable in enumerate(analysed_variables):
            result_table.add_column(Column(
                analysed_averages_all[:, index],
                name=variable
            ))
            result_table.add_column(Column(
                analysed_std_all[:, index],
                name=variable+"_err"
            ))
        result_table.write(OUT_DIR + RESULT_FILE)

        best_model_table = Table()
        best_model_table.add_column(Column(
            obs_table["id"].data,
            name="observation_id"
        ))
        best_model_table.add_column(Column(
            best_chi2_all,
            name="chi_square"
        ))
        best_model_table.add_column(Column(
            best_chi2_red_all,
            name="reduced_chi_square"
        ))
        if sed.sfh is not None:
            best_model_table.add_column(Column(
                normalisation_factors_all,
                name="galaxy_mass",
                unit="Msun"
            ))

        for index, name in enumerate(sed.info.keys()):
            column = Column([best_variables[index]
                             for best_variables in best_variables_all],
                            name=name)
            if name in sed.mass_proportional_info:
                column *= normalisation_factors_all
            best_model_table.add_column(column)

        for index, name in enumerate(filters):
            column = Column(best_fluxes[:, index] * normalisation_factors_all,
                            name=name, unit='mJy')
            best_model_table.add_column(column)

        best_model_table.write(OUT_DIR + BEST_MODEL_FILE)


# AnalysisModule to be returned by get_module
Module = PdfAnalysis
