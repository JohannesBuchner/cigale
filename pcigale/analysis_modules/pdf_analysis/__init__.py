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

from collections import OrderedDict
import ctypes
import multiprocessing as mp
from multiprocessing.sharedctypes import RawArray
import time

import numpy as np

from ...utils import read_table
from .. import AnalysisModule, complete_obs_table
from .utils import save_table_analysis, save_table_best, backup_dir
from ...warehouse import SedWarehouse
from ...data import Database
from .workers import sed as worker_sed
from .workers import init_sed as init_worker_sed
from .workers import init_analysis as init_worker_analysis
from .workers import analysis as worker_analysis
from ..utils import find_changed_parameters
import pcigale.analysis_modules.myglobals as gbl

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

        The analysis is done in two steps which can both run on multiple
        processors to run faster. The first step is to compute all the fluxes
        associated with each model as well as ancillary data such as the SED
        information. The second step is to carry out the analysis of each
        object, considering all models at once.

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

        print("Initialising the analysis module... ")

        # Rename the output directory if it exists
        backup_dir(OUT_DIR)

        # Get the parameters. They are stored in a shared module so that this
        # can be accessed by subprocesses during the model generation and
        # analysis. This avoids transmitting the exact same data all over again
        gbl.analysed_variables = parameters["analysed_variables"]
        gbl.creation_modules = creation_modules
        gbl.creation_modules_params = creation_modules_params
        gbl.save_best_sed = parameters["save_best_sed"].lower() == "true"
        gbl.save_chi2 = parameters["save_chi2"].lower() == "true"
        gbl.save_pdf = parameters["save_pdf"].lower() == "true"
        gbl.n_models = len(creation_modules_params)
        gbl.n_variables = len(gbl.analysed_variables)

        # Get the needed filters in the pcigale database. We use an ordered
        # dictionary because we need the keys to always be returned in the
        # same order. We also put the filters in the shared modules as they
        # are needed to compute the fluxes during the models generation.
        with Database() as base:
            gbl.filters = OrderedDict([(name, base.get_filter(name))
                                       for name in column_list
                                       if not name.endswith('_err')])
        gbl.n_filters = len(gbl.filters)

        # Read the observation table and complete it by adding error where
        # none is provided and by adding the systematic deviation.
        obs_table = complete_obs_table(read_table(data_file), column_list,
                                       gbl.filters, TOLERANCE)
        gbl.n_obs = len(obs_table)

        print("Computing the models fluxes...")

        # First we get a dummy sed to obtain some information on the
        # parameters, such as which ones are dependent on the normalisation
        # factor of the fit.

        # The SED warehouse is used to retrieve SED corresponding to some
        # modules and parameters.
        gbl.warehouse = SedWarehouse(cache_type=parameters["storage_type"])

        sed = gbl.warehouse.get_sed(gbl.creation_modules,
                                    creation_modules_params[0])
        gbl.info_keys = list(sed.info.keys())
        gbl.n_keys = len(gbl.info_keys)
        gbl.mass_proportional_info = sed.mass_proportional_info
        gbl.has_sfh = sed.sfh is not None

        # Arrays where we store the data related to the models. For memory
        # efficiency reasons, we use RawArrays that will be passed in argument
        # to the pool. Each worker will fill a part of the RawArrays. It is
        # important that there is no conflict and that two different workers do
        # not write on the same section.
        # We put the shape in a tuple along with the RawArray because workers
        # need to know the shape to create the numpy array from the RawArray.
        model_redshifts = (RawArray(ctypes.c_double, gbl.n_models),
                           (gbl.n_models))
        model_fluxes = (RawArray(ctypes.c_double,
                                 gbl.n_models * gbl.n_filters),
                        (gbl.n_models, gbl.n_filters))
        model_variables = (RawArray(ctypes.c_double,
                                    gbl.n_models * gbl.n_variables),
                           (gbl.n_models, gbl.n_variables))

        # In order not to clog the warehouse memory with SEDs that will not be
        # used again during the SED generation, we identify which parameter has
        # changed, invalidating  models with this parameter
        changed_pars = find_changed_parameters(creation_modules_params)

        # Compute the SED fluxes and ancillary data in parallel
        initargs = (model_redshifts, model_fluxes, model_variables,
                    mp.Value('i', 0), time.time())
        if cores == 1:  # Do not create a new process
            init_worker_sed(*initargs)
            for changed_par, idx in zip(changed_pars, range(gbl.n_models)):
                worker_sed(changed_par, idx)
        else:
            with mp.Pool(processes=cores, initializer=init_worker_sed,
                         initargs=initargs) as pool:
                pool.starmap(worker_sed, zip(changed_pars,
                                             range(gbl.n_models)))

        print('\nAnalysing models...')

        # Analysis of each object in parallel.
        initargs = (model_redshifts, model_fluxes, model_variables,
                    mp.Value('i', 0), time.time())
        if cores == 1:  # Do not create a new process
            init_worker_analysis(*initargs)
            items = []
            for obs in obs_table:
                items.append(worker_analysis(obs))
        else:
            with mp.Pool(processes=cores, initializer=init_worker_analysis,
                         initargs=initargs) as pool:
                items = pool.starmap(worker_analysis, zip(obs_table))

        # Local arrays where to unpack the results of the analysis
        analysed_averages = np.empty((gbl.n_obs, gbl.n_variables))
        analysed_std = np.empty_like(analysed_averages)
        best_fluxes = np.empty((gbl.n_obs, gbl.n_filters))
        best_parameters = [None] * gbl.n_obs
        best_normalisation_factors = np.empty(gbl.n_obs)
        best_chi2 = np.empty_like(best_normalisation_factors)
        best_chi2_red = np.empty_like(best_normalisation_factors)

        for item_idx, item in enumerate(items):
            analysed_averages[item_idx, :] = item[0]
            analysed_std[item_idx, :] = item[1]
            best_fluxes[item_idx, :] = item[2]
            best_parameters[item_idx] = item[3]
            best_normalisation_factors[item_idx] = item[4]
            best_chi2[item_idx] = item[5]
            best_chi2_red[item_idx] = item[6]

        print("\nSaving results...")

        save_table_analysis(obs_table['id'], gbl.analysed_variables,
                            analysed_averages, analysed_std)
        save_table_best(obs_table['id'], best_chi2, best_chi2_red,
                        best_normalisation_factors, best_parameters,
                        best_fluxes, gbl.filters)

        print("Run completed!")

# AnalysisModule to be returned by get_module
Module = PdfAnalysis
