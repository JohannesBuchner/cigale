# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013-2014 Institute of Astronomy
# Copyright (C) 2013-2014 Yannick Roehlly <yannick@iaora.eu>
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly & Médéric Boquien

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
from collections import OrderedDict
from datetime import datetime
from itertools import repeat
import multiprocessing as mp
from ...utils import read_table
from .. import AnalysisModule, complete_obs_table
from .utils import save_table_analysis, save_table_best
from ...warehouse import SedWarehouse
from ...data import Database
from .workers import sed as worker_sed
from .workers import analysis as worker_analysis

# Tolerance threshold under which any flux or error is considered as 0.
TOLERANCE = 1e-12
# Limit the redshift to this number of decimals
REDSHIFT_DECIMALS = 2
# Directory where the output files are stored
OUT_DIR = "out/"


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
                creation_modules_params, parameters, cores):
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
        core: integer
            Number of cores to run the analysis on

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
        save = (parameters["save_best_sed"].lower() == "true",
                parameters["save_chi2"].lower() == "true",
                parameters["save_pdf"].lower() == "true")

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
        model_info = [None] * len(creation_modules_params)

        # The SED warehouse is used to retrieve SED corresponding to some
        # modules and parameters.
        with SedWarehouse(cache_type=parameters["storage_type"]) as warehouse,\
                mp.Pool(processes=cores) as pool:
            # First we get a dummy sed to obtain some information on the
            # parameters, such as which ones are dependent on the normalisation
            # factor of the fit.
            sed = warehouse.get_sed(creation_modules,
                                    creation_modules_params[0])
            items = pool.starmap(worker_sed, zip(repeat(warehouse),
                                                 repeat(creation_modules),
                                                 creation_modules_params,
                                                 repeat(analysed_variables),
                                                 repeat(filters)))
            pool.close()
            pool.join()
            for idx_item, item in enumerate(items):
                model_fluxes[idx_item, :] = item[0]
                model_variables[idx_item, :] = item[1]
                model_redshift[idx_item] = item[2]
                model_info[idx_item] = item[3]

            del items
        # Mask the invalid fluxes
        model_fluxes = np.ma.masked_less(model_fluxes, -90)

        ##################################################################
        # Observations to models comparison                              #
        ##################################################################

        # We compute the χ² only for models with the closest redshift. We
        # extract model fluxes and information into arrays dedicated to a
        # given observation.

        # This is a tricky part here. The data arrays we have computed are for
        # all redshifts. However, for memory efficiency concerns, we want to
        # transmit to the worker arrays containing only data corresponding to
        # the redshift of the observed object. The idea here is that we will
        # constructs these arrays thanks to generators. Thw basic algorithm is
        # the following:

        # 1) We compute the set of unique redshifts (as we do not have access
        # to the configuration file)
        # 2) We build a dictionary containing the indices of the models at a
        # given redshift
        # 3) We find for each observation, which is the closest redshift
        # 4) We build the generators slicing the input arrays with the indices
        # corresponding to the observation redshift
        redshifts = np.unique(model_redshift)
        w_redshifts = {redshift: model_redshift == redshift
                       for redshift in redshifts}
        closest_redshifts = [redshifts[np.abs(obs_redshift -
                                       redshifts).argmin()]
                             for obs_redshift in obs_table['redshift']]
        model_fluxes_obs = (model_fluxes[w_redshifts[redshift], :]
                            for redshift in closest_redshifts)
        model_info_obs = (np.array(model_info)[w_redshifts[redshift]]
                          for redshift in closest_redshifts)
        model_variables_obs = (model_variables[w_redshifts[redshift], :]
                               for redshift in closest_redshifts)

        with mp.Pool(processes=cores) as pool:
            items = pool.starmap(worker_analysis,
                                 zip(obs_table,
                                     model_fluxes_obs,
                                     model_variables_obs,
                                     model_info_obs,
                                     repeat(filters),
                                     repeat(sed),
                                     repeat(analysed_variables),
                                     repeat(creation_modules),
                                     repeat(creation_modules_params),
                                     repeat(save)))

            pool.close()
            pool.join()

        analysed_averages = np.empty((len(items), len(analysed_variables)))
        analysed_std = np.empty_like(analysed_averages)
        chi2_ = np.empty(len(obs_table))
        chi2_red = np.empty_like(chi2_)
        normalisation_factors = np.empty_like(chi2_)
        variables = [None] * len(items)
        fluxes = np.empty((len(items), len(filters)))

        for item_idx, item in enumerate(items):
            analysed_averages[item_idx, :] = item[0]
            analysed_std[item_idx, :] = item[1]
            chi2_[item_idx] = item[2]
            chi2_red[item_idx] = item[3]
            normalisation_factors[item_idx] = item[4]
            variables[item_idx] = item[5]
            fluxes[item_idx, :] = item[6]

        del items

        save_table_analysis(obs_table['id'], analysed_variables,
                            analysed_averages, analysed_std)
        save_table_best(obs_table['id'], chi2_, chi2_red,
                        normalisation_factors, variables, fluxes, filters, sed)

# AnalysisModule to be returned by get_module
Module = PdfAnalysis
