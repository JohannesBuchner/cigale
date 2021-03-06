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
from .utils import save_results, analyse_chi2
from ...warehouse import SedWarehouse
from .workers import sed as worker_sed
from .workers import init_sed as init_worker_sed
from .workers import init_analysis as init_worker_analysis
from .workers import analysis as worker_analysis
from ..utils import backup_dir
from ...handlers.parameters_handler import ParametersHandler


# Tolerance threshold under which any flux or error is considered as 0.
TOLERANCE = 1e-12


class PdfAnalysis(AnalysisModule):
    """PDF analysis module"""

    parameter_list = OrderedDict([
        ("variables", (
            "cigale_string_list()",
            "List of the physical properties to estimate. Leave empty to "
            "analyse all the physical properties (not recommended when there "
            "are many models).",
            ["sfh.sfr", "sfh.sfr10Myrs", "sfh.sfr100Myrs"]
        )),
        ("save_best_sed", (
            "boolean()",
            "If true, save the best SED for each observation to a file.",
            False
        )),
        ("save_chi2", (
            "boolean()",
            "If true, for each observation and each analysed variable save "
            "the reduced chi2.",
            False
        )),
        ("save_pdf", (
            "boolean()",
            "If true, for each observation and each analysed variable save "
            "the probability density function.",
            False
        )),
        ("lim_flag", (
            "boolean()",
            "If true, for each object check whether upper limits are present "
            "and analyse them.",
            False
        )),
        ("mock_flag", (
            "boolean()",
            "If true, for each object we create a mock object "
            "and analyse them.",
            False
        ))
    ])

    def process(self, conf):
        """Process with the psum analysis.

        The analysis is done in two steps which can both run on multiple
        processors to run faster. The first step is to compute all the fluxes
        associated with each model as well as ancillary data such as the SED
        information. The second step is to carry out the analysis of each
        object, considering all models at once.

        Parameters
        ----------
        conf: dictionary
            Contents of pcigale.ini in the form of a dictionary

        """
        np.seterr(invalid='ignore')

        print("Initialising the analysis module... ")

        # Rename the output directory if it exists
        backup_dir()

        # Initalise variables from input arguments.
        variables = conf['analysis_params']["variables"]
        variables_nolog = [variable[:-4] if variable.endswith('_log') else
                           variable for variable in variables]
        n_variables = len(variables)
        save = {key: conf['analysis_params']["save_{}".format(key)] for key in
                ["best_sed", "chi2", "pdf"]}
        lim_flag = conf['analysis_params']["lim_flag"]

        filters = [name for name in conf['bands'] if not
                   name.endswith('_err')]
        n_filters = len(filters)

        # Read the observation table and complete it by adding error where
        # none is provided and by adding the systematic deviation.
        obs_table = complete_obs_table(read_table(conf['data_file']),
                                       conf['bands'], filters, TOLERANCE,
                                       lim_flag)
        n_obs = len(obs_table)

        z = np.array(conf['sed_modules_params']['redshifting']['redshift'])

        # The parameters handler allows us to retrieve the models parameters
        # from a 1D index. This is useful in that we do not have to create
        # a list of parameters as they are computed on-the-fly. It also has
        # nice goodies such as finding the index of the first parameter to
        # have changed between two indices or the number of models.
        params = ParametersHandler(conf)
        n_params = params.size

        # Retrieve an arbitrary SED to obtain the list of output parameters
        warehouse = SedWarehouse()
        sed = warehouse.get_sed(conf['sed_modules'], params.from_index(0))
        info = list(sed.info.keys())
        info.sort()
        n_info = len(info)
        del warehouse, sed

        print("Computing the models fluxes...")

        # Arrays where we store the data related to the models. For memory
        # efficiency reasons, we use RawArrays that will be passed in argument
        # to the pool. Each worker will fill a part of the RawArrays. It is
        # important that there is no conflict and that two different workers do
        # not write on the same section.
        # We put the shape in a tuple along with the RawArray because workers
        # need to know the shape to create the numpy array from the RawArray.
        model_fluxes = (RawArray(ctypes.c_double, n_params * n_filters),
                        (n_params, n_filters))
        model_variables = (RawArray(ctypes.c_double, n_params * n_variables),
                           (n_params, n_variables))

        initargs = (params, filters, variables_nolog, model_fluxes,
                    model_variables, time.time(), mp.Value('i', 0))
        if conf['cores'] == 1:  # Do not create a new process
            init_worker_sed(*initargs)
            for idx in range(n_params):
                worker_sed(idx)
        else:  # Compute the models in parallel
            with mp.Pool(processes=conf['cores'], initializer=init_worker_sed,
                         initargs=initargs) as pool:
                pool.map(worker_sed, range(n_params))

        print("\nAnalysing models...")

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

        initargs = (params, filters, variables, z, model_fluxes,
                    model_variables, time.time(), mp.Value('i', 0),
                    analysed_averages, analysed_std, best_fluxes,
                    best_parameters, best_chi2, best_chi2_red, save, lim_flag,
                    n_obs)
        if conf['cores'] == 1:  # Do not create a new process
            init_worker_analysis(*initargs)
            for idx, obs in enumerate(obs_table):
                worker_analysis(idx, obs)
        else:  # Analyse observations in parallel
            with mp.Pool(processes=conf['cores'],
                         initializer=init_worker_analysis,
                         initargs=initargs) as pool:
                pool.starmap(worker_analysis, enumerate(obs_table))

        analyse_chi2(best_chi2_red)

        print("\nSaving results...")

        save_results("results", obs_table['id'], variables, analysed_averages,
                     analysed_std, best_chi2, best_chi2_red, best_parameters,
                     best_fluxes, filters, info)

        if conf['analysis_params']['mock_flag'] is True:

            print("\nMock analysis...")

            # For the mock analysis we do not save the ancillary files
            for k in save:
                save[k] = False

            obs_fluxes = np.array([obs_table[name] for name in filters]).T
            obs_errors = np.array([obs_table[name + "_err"] for name in
                                   filters]).T
            mock_fluxes = obs_fluxes.copy()
            bestmod_fluxes = np.ctypeslib.as_array(best_fluxes[0])
            bestmod_fluxes = bestmod_fluxes.reshape(best_fluxes[1])
            wdata = np.where((obs_fluxes > TOLERANCE) &
                             (obs_errors > TOLERANCE))
            mock_fluxes[wdata] = np.random.normal(bestmod_fluxes[wdata],
                                                  obs_errors[wdata])

            mock_table = obs_table.copy()
            for idx, name in enumerate(filters):
                mock_table[name] = mock_fluxes[:, idx]

            initargs = (params, filters, variables, z, model_fluxes,
                        model_variables, time.time(), mp.Value('i', 0),
                        analysed_averages, analysed_std, best_fluxes,
                        best_parameters, best_chi2, best_chi2_red, save,
                        lim_flag, n_obs)
            if conf['cores'] == 1:  # Do not create a new process
                init_worker_analysis(*initargs)
                for idx, mock in enumerate(mock_table):
                    worker_analysis(idx, mock)
            else:  # Analyse observations in parallel
                with mp.Pool(processes=conf['cores'],
                             initializer=init_worker_analysis,
                             initargs=initargs) as pool:
                    pool.starmap(worker_analysis, enumerate(mock_table))

            print("\nSaving results...")

            save_results("results_mock", mock_table['id'], variables,
                         analysed_averages, analysed_std, best_chi2,
                         best_chi2_red, best_parameters, best_fluxes, filters,
                         info)

        print("Run completed!")

# AnalysisModule to be returned by get_module
Module = PdfAnalysis
