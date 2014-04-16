# -*- coding: utf-8 -*-
# Copyright (C) 2014 Laboratoire d'Astrophysique de Marseille, AMU
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013-2014 Institute of Astronomy
# Copyright (C) 2013-2014 Yannick Roehlly <yannick@iaora.eu>
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly, Médéric Boquien & Denis Burgarella

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
from .utils import save_table_analysis, save_table_best
from ...warehouse import SedWarehouse
from ...data import Database
from .workers import sed as worker_sed
from .workers import init_sed as init_worker_sed
from .workers import init_analysis as init_worker_analysis
from .workers import analysis as worker_analysis
from ..utils import ParametersHandler, backup_dir

# Tolerance threshold under which any flux or error is considered as 0.
TOLERANCE = 1e-12
# Limit the redshift to this number of decimals
REDSHIFT_DECIMALS = 2


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
            "the reduced chi2.",
            False
        )),
        ("save_pdf", (
            "boolean",
            "If true, for each observation and each analysed variable save "
            "the probability density function.",
            False
        )),
        ("lim_flag", (
            "boolean",
            "If true, for each object check whether upper limits are present "
            "and analyse them.",
            False
        )),
        ("storage_type", (
            "string",
            "Type of storage used to cache the generate SED.",
            "memory"
        ))
    ])

    def process(self, data_file, column_list, creation_modules,
                creation_modules_params, config, cores):
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
        config: dictionary
            Dictionary containing the configuration.
        core: integer
            Number of cores to run the analysis on

        """

        print("Initialising the analysis module... ")

        # Rename the output directory if it exists
        backup_dir()

        # Initalise variables from input arguments.
        analysed_variables = config["analysed_variables"]
        n_variables = len(analysed_variables)
        save = {key:config["save_{}".format(key)].lower() == "true"
                for key in ["best_sed", "chi2", "pdf"]}
        lim_flag = config["lim_flag"].lower() == "true"

        # Get the needed filters in the pcigale database. We use an ordered
        # dictionary because we need the keys to always be returned in the
        # same order. We also put the filters in the shared modules as they
        # are needed to compute the fluxes during the models generation.
        with Database() as base:
            filters = OrderedDict([(name, base.get_filter(name))
                                   for name in column_list
                                   if not name.endswith('_err')])
        n_filters = len(filters)

        # Read the observation table and complete it by adding error where
        # none is provided and by adding the systematic deviation.
        obs_table = complete_obs_table(read_table(data_file), column_list,
                                       filters, TOLERANCE, lim_flag)
        n_obs = len(obs_table)

        w_redshifting = creation_modules.index('redshifting')
        if creation_modules_params[w_redshifting]['redshift'] == ['']:
            z = np.unique(np.around(obs_table['redshift'],
                                    decimals=REDSHIFT_DECIMALS))
            creation_modules_params[w_redshifting]['redshift'] = z
            del z

        # The parameters handler allows us to retrieve the models parameters
        # from a 1D index. This is useful in that we do not have to create
        # a list of parameters as they are computed on-the-fly. It also has
        # nice goodies such as finding the index of the first parameter to
        # have changed between two indices or the number of models.
        params = ParametersHandler(creation_modules, creation_modules_params)
        n_params = params.size

        # Retrieve an arbitrary SED to obtain the list of output parameters
        warehouse = SedWarehouse(cache_type=config["storage_type"])
        sed = warehouse.get_sed(creation_modules, params.from_index(0))
        info = sed.info
        n_info = len(sed.info)
        del warehouse, sed

        print("Computing the models fluxes...")

        # Arrays where we store the data related to the models. For memory
        # efficiency reasons, we use RawArrays that will be passed in argument
        # to the pool. Each worker will fill a part of the RawArrays. It is
        # important that there is no conflict and that two different workers do
        # not write on the same section.
        # We put the shape in a tuple along with the RawArray because workers
        # need to know the shape to create the numpy array from the RawArray.
        model_redshifts = (RawArray(ctypes.c_double, n_params),
                           (n_params))
        model_fluxes = (RawArray(ctypes.c_double,
                                 n_params * n_filters),
                        (n_params, n_filters))
        model_variables = (RawArray(ctypes.c_double,
                                    n_params * n_variables),
                           (n_params, n_variables))

        initargs = (params, filters, analysed_variables, model_redshifts,
                    model_fluxes, model_variables, time.time(),
                    mp.Value('i', 0))
        if cores == 1:  # Do not create a new process
            init_worker_sed(*initargs)
            for idx in range(n_params):
                worker_sed(idx)
        else:  # Analyse observations in parallel
            with mp.Pool(processes=cores, initializer=init_worker_sed,
                         initargs=initargs) as pool:
                pool.map(worker_sed, range(n_params))

        print('\nAnalysing models...')

        # We use RawArrays for the same reason as previously
        analysed_averages = (RawArray(ctypes.c_double, n_obs * n_variables),
                           (n_obs, n_variables))
        analysed_std = (RawArray(ctypes.c_double, n_obs * n_variables),
                           (n_obs, n_variables))
        best_fluxes = (RawArray(ctypes.c_double, n_obs * n_filters),
                           (n_obs, n_filters))
        best_parameters = (RawArray(ctypes.c_double, n_obs * n_info),
                           (n_obs, n_info))
        best_chi2 = (RawArray(ctypes.c_double, n_obs), (n_obs))
        best_chi2_red = (RawArray(ctypes.c_double, n_obs), (n_obs))

        initargs = (params, filters, analysed_variables, model_redshifts,
                    model_fluxes, model_variables, time.time(),
                    mp.Value('i', 0), analysed_averages, analysed_std,
                    best_fluxes, best_parameters, best_chi2, best_chi2_red,
                    save, lim_flag, n_obs)
        if cores == 1:  # Do not create a new process
            init_worker_analysis(*initargs)
            for idx, obs in enumerate(obs_table):
                worker_analysis(idx, obs)
        else:  # Analyse observations in parallel
            with mp.Pool(processes=cores, initializer=init_worker_analysis,
                         initargs=initargs) as pool:
                pool.starmap(worker_analysis, enumerate(obs_table))

        print("\nSaving results...")

        save_table_analysis(obs_table['id'], analysed_variables,
                            analysed_averages, analysed_std)
        save_table_best(obs_table['id'], best_chi2, best_chi2_red,
                        best_parameters, best_fluxes, filters, info)

        print("Run completed!")

# AnalysisModule to be returned by get_module
Module = PdfAnalysis
