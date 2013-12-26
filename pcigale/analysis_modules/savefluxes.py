# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly <yannick.roehlly@oamp.fr>

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


class Module(AnalysisModule):
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
                creation_modules_params, redshift_module,
                redshift_configuration, parameters):
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
        redshift_module_name : string
            Name of the module used to redshift the SED.
        redshift_configuration : dictionary
            Configuration dictionary for the module used to redshift the SED.
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
        base = Database()
        filter_list = [base.get_filter(name) for name in filter_names]
        base.close()

        # Columns of the output table
        out_columns = []
        for module_param_list in zip(creation_modules,
                                     creation_modules_params[0]):
            for module_param in product([module_param_list[0]],
                                        module_param_list[1].keys()):
                out_columns.append(".".join(module_param))
        out_columns += filter_names

        # Content of the output table
        out_rows = []

        # Open the warehouse
        sed_warehouse = SedWarehouse(
            cache_type=parameters["storage_type"])

        # We loop over all the possible theoretical SEDs
        progress_bar = ProgressBar(maxval=len(creation_modules_params)).start()
        for model_index, parameters in enumerate(creation_modules_params):
            sed = sed_warehouse.get_sed(creation_modules, parameters)

            row = []

            # Add the parameter values to the row. Some parameters are array
            # so we must join their content.
            for module_param in parameters:
                for value in module_param.values():
                    if type(value) == list:
                        value = ".".join(value)
                    row.append(value)

            # Add the flux in each filter to the row
            row += [sed.compute_fnu(filter.trans_table,
                                    filter.effective_wavelength,
                                    0)
                    for filter in filter_list]

            out_rows.append(row)

            progress_bar.update(model_index + 1)

        progress_bar.finish()

        # The zip call is to convert the list of rows to a list of columns.
        out_table = Table(zip(*out_rows), names=out_columns)
        out_table.write(out_file, format=out_format)
