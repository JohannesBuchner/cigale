# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

"""
Save fluxes analysis module
===========================

This module does not perform a statistical analysis. It computes and save the
fluxes in a set of filters for all the possible combinations of input SED
parameters.

The data file is used only to get the list of fluxes to be computed.

"""
from collections import OrderedDict
import ctypes
from datetime import datetime
from itertools import product, repeat
import os
import multiprocessing as mp
from multiprocessing.sharedctypes import RawArray
import time

import numpy as np

from .. import AnalysisModule
from ...data import Database
from ..utils import ParametersHandler, backup_dir, save_fluxes
from ...utils import read_table
from ...warehouse import SedWarehouse
from .workers import init_fluxes as init_worker_fluxes
from .workers import fluxes as worker_fluxes

# Limit the redshift to this number of decimals
REDSHIFT_DECIMALS = 2


class SaveFluxes(AnalysisModule):
    """Save fluxes analysis module

    This module saves a table containing all the parameters and desired fluxes
    for all the computed models.

    """

    parameter_list = OrderedDict([
        ("output_file", (
            "string",
            "Name of the output file that contains the parameters of the model(s)"
            "and the flux densities in the bands",
            "computed_fluxes.txt"
        )),
        ("save_sed", (
            "boolean",
            "If True, save the generated spectrum for each model.",
            "False"
        )),
        ("output_format", (
            "string",
            "Format of the output file. Any format supported by astropy.table "
            "e.g. votable or ascii.",
            "ascii"
        ))
    ])

    def process(self, data_file, column_list, creation_modules,
                creation_modules_params, parameters, cores):
        """Process with the savedfluxes analysis.

        All the possible theoretical SED are created and the fluxes in the
        filters from the column_list are computed and saved to a table,
        alongside the parameter values.

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
        cores: integer
            Number of cores to run the analysis on

        """

        # Rename the output directory if it exists
        backup_dir()

        out_file = parameters["output_file"]
        out_format = parameters["output_format"]
        save_sed = parameters["save_sed"].lower() == "true"

        # Get the needed filters in the pcigale database. We use an ordered
        # dictionary because we need the keys to always be returned in the
        # same order. We also put the filters in the shared modules as they
        # are needed to compute the fluxes during the models generation.
        with Database() as base:
            filters = OrderedDict([(name, base.get_filter(name))
                                   for name in column_list
                                   if not name.endswith('_err')])
        n_filters = len(filters)

        w_redshifting = creation_modules.index('redshifting')
        if creation_modules_params[w_redshifting]['redshift'] == ['']:
            obs_table = read_table(data_file)
            z = np.unique(np.around(obs_table['redshift'],
                                    decimals=REDSHIFT_DECIMALS))
            creation_modules_params[w_redshifting]['redshift'] = z
            del obs_table, z

        # The parameters handler allows us to retrieve the models parameters
        # from a 1D index. This is useful in that we do not have to create
        # a list of parameters as they are computed on-the-fly. It also has
        # nice goodies such as finding the index of the first parameter to
        # have changed between two indices or the number of models.
        params = ParametersHandler(creation_modules, creation_modules_params)
        n_params = params.size

        # Retrieve an arbitrary SED to obtain the list of output parameters
        warehouse = SedWarehouse()
        sed = warehouse.get_sed(creation_modules, params.from_index(0))
        info = sed.info
        n_info = len(sed.info)
        del warehouse, sed

        model_fluxes = (RawArray(ctypes.c_double,
                                 n_params * n_filters),
                        (n_params, n_filters))
        model_parameters = (RawArray(ctypes.c_double,
                                     n_params * n_info),
                        (n_params, n_info))

        initargs = (params, filters, save_sed, model_fluxes,
                    model_parameters, time.time(), mp.Value('i', 0))
        if cores == 1:  # Do not create a new process
            init_worker_fluxes(*initargs)
            for idx in range(n_params):
                worker_fluxes(idx)
        else:  # Analyse observations in parallel
            with mp.Pool(processes=cores, initializer=init_worker_fluxes,
                         initargs=initargs) as pool:
                pool.map(worker_fluxes, range(n_params))

        save_fluxes(model_fluxes, model_parameters, filters, info, out_file,
                    out_format=out_format)

# AnalysisModule to be returned by get_module
Module = SaveFluxes
