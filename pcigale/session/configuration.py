# -*- coding: utf-8 -*-
"""
Copyright (C) 2012 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""


import atpy
import configobj
import pkg_resources
import pkgutil
import numpy as np
from textwrap import wrap
from ..data import Database
from ..sed.modules import common as modules
from ..stats import common as analysis


def list_modules(package_name):
    """Lists the modules available in a package

    Parametres
    ----------
    package_name : string
        Name of the package (e.g. pcigale.sed.modules).

    Returns
    -------
    module_name : array of strings
        List of the available modules.

    """
    directory = pkg_resources.resource_filename(package_name, '')
    module_names = [name for _, name, _ in pkgutil.iter_modules([directory])]
    if 'common' in module_names:
        module_names.remove('common')

    return module_names


def read_array_description(description):
    """Read an array description from configuration file

    Interpret a string (from the configuration file) as a Numpy array
    definition. The string can be either a Numpy command beginning with 'np.'
    and evaluated (may be dangerous), or a range definition beginning with
    'range' and followed by the start value, the stop value and the step (the
    stop value is included if reached), or a list of values separated with
    spaces.

    The returned array is always an array of floats.

    Parametres
    ----------
    description : string
        The description to be evaluated.

    Returns
    -------
    array : array of floats
        The evaluated array.

    """

    if description[:3] == 'np.':
        array = eval(description)
    else:
        value_list = description.split(' ')

        if value_list[0] == 'range':
            start = float(value_list[1])
            step = float(value_list[3])
            stop = float(value_list[2]) + step  # to include stop value
            array = np.arange(start, stop, step, float)
        else:
            array = np.array(value_list, dtype=float)

    return array


class Configuration(object):
    """This class manages the configuration of pcigale.
    """

    def __init__(self, filename="pcigale.ini"):
        """Initialise a pcigale configuration.

        Parametres
        ----------
        filename : string
            Name of the configuration file (pcigale.conf by default).

        """
        self.config = configobj.ConfigObj(filename,
                                          write_empty_values=True,
                                          indent_type='  ')

    def create_blank_conf(self):
        """Create the initial configuration file

        Write the initial pcigale configuration file where the user can state
        which data file to use, which modules to use for the SED creation, as
        well as the method selected for statistical analysis.

        """

        self.config['data_file'] = ""
        self.config.comments['data_file'] = wrap(
            "File containing the observation data to be fitted. Each flux "
            "column must have the name of the corresponding filter, the "
            "error columns are suffixed with '_err'. The values must be "
            "in mJy.")

        self.config['sed_modules'] = []
        self.config.comments['sed_modules'] = [""] + wrap(
            "Order of the modules use for SED creation. Available modules : "
            + ', '.join(list_modules('pcigale.sed.modules')) + ".")

        self.config['analysis_method'] = ""
        self.config.comments['analysis_method'] = [""] + wrap(
            "Method used for statistical analysis. Available methods: "
            + ', '.join(list_modules('pcigale.stats')) + ".")

        self.config.write()

    def generate_conf(self):
        """Generate the full configuration file

        Reads the user entries in the initial configuration file and add the
        configuration options of all selected modules as well as the filter
        selection based on the filters identified in the data table file.

        """

        # Getting the list of the filters available in pcigale database
        base = Database()
        filter_list = base.get_filter_list()[0]
        base.close()

        # Finding the known filters in the data table
        obs_table = atpy.Table(self.config['data_file'])
        column_list = []
        for column in obs_table.columns:
            filter_name = column[:-4] if column.endswith('_err') else column
            if filter_name in filter_list:
                column_list.append(column)

        # Check that we don't have an error column without the associated flux
        for column in column_list:
            if column.endswith('_err') and (column[:-4] not in column_list):
                raise StandardError("The observation table as a {} column "
                                    "but no {} column.".format(column,
                                                               column[:-4]))

        self.config['column_list'] = column_list
        self.config.comments['column_list'] = [""] + wrap(
            "List of the columns in the observation data file to use for "
            "the fitting.")

        # SED creation modules configurations. For each module, we generate
        # the configuration section from its parametre list.
        self.config['sed_creation_modules'] = {}
        self.config.comments['sed_creation_modules'] = ["", ""] + wrap(
            "Configuration of the SED creation modules.")

        for module_name in self.config['sed_modules']:
            self.config["sed_creation_modules"][module_name] = {}
            sub_config = self.config["sed_creation_modules"][module_name]

            for name, (typ, unit, description, default) in \
                    modules.get_module(module_name).parametre_list.items():
                if default is None:
                    default = ''
                sub_config[name] = default
                sub_config.comments[name] = wrap(description)

            self.config['sed_creation_modules'].comments[module_name] = [
                modules.get_module(module_name).comments]

        # Configuration for the analysis method
        self.config['analysis_configuration'] = {}
        self.config.comments['analysis_configuration'] = ["", ""] + wrap(
            "Configuration of the statistical analysis method.")
        module_name = self.config['analysis_method']
        for name, (typ, unit, desc, default) in \
                analysis.get_module(module_name).parametre_list.items():
            if default is None:
                default = ''
            self.config['analysis_configuration'][name] = default
            self.config['analysis_configuration'].comments[name] = wrap(desc)

        self.config.write()
