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

import os
from itertools import product
from collections import OrderedDict
from datetime import datetime
from astropy.table import Table
from progressbar import ProgressBar
from . import AnalysisModule
from ..warehouse import SedWarehouse
from ..data import Database


class SaveFluxes(AnalysisModule):
    """Save fluxes analysis module

    This module saves a table containing all the parameters and desired fluxes
    for all the computed models.

    """

    parameter_list = OrderedDict([
        ("output_file", (
            "string",
            "Name of the output file.",
            "computed_fluxes.xml"
        )),
        ("output_format", (
            "string",
            "Format of the output file. Any format supported by astropy.table "
            "e.g. votable or ascii.",
            "votable"
        )),
        ("storage_type", (
            "string",
            "Type of storage used to cache the generate SED.",
            "memory"
        ))
    ])

    def process(self, data_file, column_list, creation_modules,
                creation_modules_params, parameters):
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

        """

        out_file = parameters["output_file"]
        out_format = parameters["output_format"]

        # If the output file already exists make a copy.
        if os.path.isfile(out_file):
            new_name = datetime.now().strftime("%Y%m%d%H%M") + "_" + out_file
            os.rename(out_file, new_name)
            print("The existing {} file was renamed to {}".format(
                out_file,
                new_name
            ))

        # Get the filters in the database
        filter_names = [name for name in column_list
                        if not name.endswith('_err')]
        with Database() as base:
            filter_list = [base.get_filter(name) for name in filter_names]

        # Content of the output table.
        # In the output table, we put the content of the sed.info dictionary
        # plus the flux in all the filters. As all the SEDs are made with the
        #  same pipeline, they should have the same sed.info dictionary keys.
        output = []

        # Open the warehouse
        sed_warehouse = SedWarehouse(
            cache_type=parameters["storage_type"])

        # We loop over all the possible theoretical SEDs
        progress_bar = ProgressBar(maxval=len(creation_modules_params)).start()
        for model_index, parameters in enumerate(creation_modules_params):
            sed = sed_warehouse.get_sed(creation_modules, parameters)

            # Take the content of the sed info dictionary.
            row = list(sed.info.values())

            # Add the flux in each filter to the row
            row += [sed.compute_fnu(filter_.trans_table,
                                    filter_.effective_wavelength)
                    for filter_ in filter_list]

            output.append(row)

            progress_bar.update(model_index + 1)

        progress_bar.finish()

        # We take the names of the columns from the last computed SED.
        out_columns = list(sed.info.keys()) + filter_names

        # The zip call is to convert the list of rows to a list of columns.
        out_table = Table(list(zip(*output)), names=out_columns)
        out_table.write(out_file, format=out_format)

# AnalysisModule to be returned by get_module
Module = SaveFluxes
