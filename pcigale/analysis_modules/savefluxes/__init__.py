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
import multiprocessing as mp
from multiprocessing.sharedctypes import RawArray
import time

import numpy as np

from .. import AnalysisModule
from ..utils import backup_dir, save_fluxes
from ...utils import read_table
from .workers import init_fluxes as init_worker_fluxes
from .workers import fluxes as worker_fluxes
from ...handlers.parameters_handler import ParametersHandler


class SaveFluxes(AnalysisModule):
    """Save fluxes analysis module

    This module saves a table containing all the parameters and desired fluxes
    for all the computed models.

    """

    parameter_list = OrderedDict([
        ("variables", (
            "cigale_string_list()",
            "List of the physical properties to save. Leave empty to save all "
            "the physical properties (not recommended when there are many "
            "models).",
            None
        )),
        ("output_file", (
            "string()",
            "Name of the output file that contains the parameters of the "
            "model(s) and the flux densities in the bands",
            "computed_fluxes.txt"
        )),
        ("save_sed", (
            "boolean()",
            "If True, save the generated spectrum for each model.",
            False
        )),
        ("output_format", (
            "string()",
            "Format of the output file. Any format supported by astropy.table "
            "e.g. votable or ascii.",
            "ascii"
        ))
    ])

    def process(self, conf):
        """Process with the savedfluxes analysis.

        All the possible theoretical SED are created and the fluxes in the
        filters from the column_list are computed and saved to a table,
        alongside the parameter values.

        Parameters
        ----------
        conf: dictionary
            Contents of pcigale.ini in the form of a dictionary
        """

        # Rename the output directory if it exists
        backup_dir()
        out_file = conf['analysis_params']['output_file']
        out_format = conf['analysis_params']['output_format']
        save_sed = conf['analysis_params']['save_sed']

        filters = [name for name in conf['column_list'] if not
                   name.endswith('_err')]
        n_filters = len(filters)

        # The parameters handler allows us to retrieve the models parameters
        # from a 1D index. This is useful in that we do not have to create
        # a list of parameters as they are computed on-the-fly. It also has
        # nice goodies such as finding the index of the first parameter to
        # have changed between two indices or the number of models.
        params = ParametersHandler(conf)
        n_params = params.size

        info = conf['analysis_params']['variables']
        n_info = len(info)

        model_fluxes = (RawArray(ctypes.c_double, n_params * n_filters),
                        (n_params, n_filters))
        model_parameters = (RawArray(ctypes.c_double, n_params * n_info),
                            (n_params, n_info))

        initargs = (params, filters, save_sed, info, model_fluxes,
                    model_parameters, time.time(), mp.Value('i', 0))
        if conf['cores'] == 1:  # Do not create a new process
            init_worker_fluxes(*initargs)
            for idx in range(n_params):
                worker_fluxes(idx)
        else:  # Analyse observations in parallel
            with mp.Pool(processes=conf['cores'],
                         initializer=init_worker_fluxes,
                         initargs=initargs) as pool:
                pool.map(worker_fluxes, range(n_params))

        save_fluxes(model_fluxes, model_parameters, filters, info, out_file,
                    out_format=out_format)

# AnalysisModule to be returned by get_module
Module = SaveFluxes
