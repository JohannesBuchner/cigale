# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013-2014 Institute of Astronomy
# Copyright (C) 2014 Yannick Roehlly <yannick@iaora.eu>
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly & Médéric Boquien

import time

import numpy as np

from ...warehouse import SedWarehouse
from ..utils import OUT_DIR


def init_fluxes(params, filters, save_sed, fluxes, info, t_begin, n_computed):
    """Initializer of the pool of processes. It is mostly used to convert
    RawArrays into numpy arrays. The latter are defined as global variables to
    be accessible from the workers.

    Parameters
    ----------
    params: ParametersHandler
        Handles the parameters from a 1D index.
    filters: List
        Contains the names of the filters to compute the fluxes.
    save_sed: boolean
        Indicates whether the SED should be saved.
    fluxes: RawArray and tuple containing the shape
        Fluxes of individual models. Shared among workers.
    n_computed: Value
        Number of computed models. Shared among workers.
    t_begin: float
        Time of the beginning of the computation.

    """
    global gbl_model_fluxes, gbl_model_info, gbl_n_computed, gbl_t_begin
    global gbl_params, gbl_previous_idx, gbl_filters, gbl_save_sed
    global gbl_warehouse, gbl_keys

    gbl_model_fluxes = np.ctypeslib.as_array(fluxes[0])
    gbl_model_fluxes = gbl_model_fluxes.reshape(fluxes[1])

    gbl_model_info = np.ctypeslib.as_array(info[0])
    gbl_model_info = gbl_model_info.reshape(info[1])

    gbl_n_computed = n_computed
    gbl_t_begin = t_begin

    gbl_params = params

    gbl_previous_idx = -1

    gbl_filters = filters

    gbl_save_sed = save_sed

    gbl_warehouse = SedWarehouse()

    gbl_keys = None

def fluxes(idx):
    """Worker process to retrieve a SED and affect the relevant data to shared
    RawArrays.

    Parameters
    ----------
    idx: int
        Index of the model to retrieve its parameters from the parameters
        handler.

    """
    global gbl_previous_idx, gbl_keys
    if gbl_previous_idx > -1:
        gbl_warehouse.partial_clear_cache(
            gbl_params.index_module_changed(gbl_previous_idx, idx))
    gbl_previous_idx = idx

    sed = gbl_warehouse.get_sed(gbl_params.modules,
                                gbl_params.from_index(idx))

    if gbl_save_sed is True:
        sed.to_votable(OUT_DIR + "{}_best_model.xml".format(idx))

    if 'sfh.age' in sed.info and sed.info['sfh.age'] > sed.info['universe.age']:
        gbl_model_fluxes[idx, :] = np.full(len(gbl_filters), np.nan)
    else:
        gbl_model_fluxes[idx, :] = np.array([sed.compute_fnu(filter_) for
                                             filter_ in gbl_filters])

    if gbl_keys is None:
        gbl_keys = list(sed.info.keys())
        gbl_keys.sort()
    gbl_model_info[idx, :] = np.array([sed.info[k] for k in gbl_keys])

    with gbl_n_computed.get_lock():
        gbl_n_computed.value += 1
        n_computed = gbl_n_computed.value
    if n_computed % 250 == 0 or n_computed == gbl_params.size:
        t_elapsed = time.time() - gbl_t_begin
        print("{}/{} models computed in {} seconds ({} models/s)".
              format(n_computed, gbl_params.size,
                     np.around(t_elapsed, decimals=1),
                     np.around(n_computed/t_elapsed, decimals=1)),
              end="\n" if n_computed == gbl_params.size else "\r")
