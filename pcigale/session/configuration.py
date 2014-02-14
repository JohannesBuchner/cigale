# -*- coding: utf-8 -*-
# Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

from astropy.table import Table
import configobj
import pkg_resources
import pkgutil
import collections
import itertools
import numpy as np
from glob import glob # To allow the use of glob() in "eval..."
from textwrap import wrap
from .tools import param_dict_combine
from ..data import Database
from .. import creation_modules
from .. import analysis_modules


def list_modules(package_name):
    """Lists the modules available in a package

    Parameters
    ----------
    package_name : string
        Name of the package (e.g. pcigale.creation_modules).

    Returns
    -------
    module_name : array of strings
        List of the available modules.

    """
    directory = pkg_resources.resource_filename(package_name, '')
    module_names = [name for _, name, _ in pkgutil.iter_modules([directory])]

    return module_names


def evaluate_description(description):
    """Evaluate a description from the config file as a list.

    The description is read from the config file by configobj that transforms
    coma separated value in a list. From this description, this function try
    to evaluate the desired list of values:
    - If the description is a string beginning with 'eval ', then its content
      (without 'eval ') is evaluated as Python code and its result returned.
      An array is expected.
    - If the description is a string beginning by 'range', the start, step and
      stop values are then expected and the range is evaluated (stop included
      if reached.
    - Then the function tries to evaluate the description as a Numpy array of
      float and returns the mere list if this fails.

    Parameters
    ----------
    description : string or list
        The description to be evaluated.

    Returns
    -------
     results : list
        The evaluated list of values.

    """
    results = description
    if type(description) == str:
        if description.startswith('eval '):
            results = eval(description[4:])
            # If the evaluation lead to a single value, we put it in a list.
            if not isinstance(results, collections.Iterable):
                results = [results]
        elif description.startswith('range '):
            start, stop, step = [float(item) for item
                                 in description[5:].split()]
            results = np.arange(start, stop+step, step)
        else:
            # We need to return a list to combine the list of possible values
            # for each parameter.
            results = [results]

    # We prefer to evaluate the parameter as a numpy array of floats if
    # possible.
    try:
        results = np.array(results, float)
    except ValueError:
        pass

    return results


class Configuration(object):
    """This class manages the configuration of pcigale.
    """

    def __init__(self, filename="pcigale.ini"):
        """Initialise a pcigale configuration.

        Parameters
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

        self.config['creation_modules'] = []
        self.config.comments['creation_modules'] = [""] + wrap(
            "Order of the modules use for SED creation. Available modules : "
            + ', '.join(list_modules('pcigale.creation_modules')) + ".")

        self.config['analysis_method'] = ""
        self.config.comments['analysis_method'] = [""] + wrap(
            "Method used for statistical analysis. Available methods: "
            + ', '.join(list_modules('pcigale.analysis_modules')) + ".")

        self.config.write()

    def generate_conf(self):
        """Generate the full configuration file

        Reads the user entries in the initial configuration file and add the
        configuration options of all selected modules as well as the filter
        selection based on the filters identified in the data table file.

        """

        # Getting the list of the filters available in pcigale database
        with Database() as base:
            filter_list = base.get_filter_list()[0]

        # Finding the known filters in the data table
        if self.config['data_file'][-4:] == 'fits':
            obs_table = Table.read(self.config['data_file'], format='fits')
        elif self.config['data_file'][-3:] == 'vot':
            obs_table = Table.read(self.config['data_file'], format='votable')
        else:
            obs_table = Table.read(self.config['data_file'], format='ascii')
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
        # the configuration section from its parameter list.
        self.config['sed_creation_modules'] = {}
        self.config.comments['sed_creation_modules'] = ["", ""] + wrap(
            "Configuration of the SED creation modules.")

        for module_name in self.config['creation_modules']:
            self.config["sed_creation_modules"][module_name] = {}
            sub_config = self.config["sed_creation_modules"][module_name]

            for name, (typ, description, default) in \
                    creation_modules.get_module(
                        module_name,
                        blank=True).parameter_list.items():
                if default is None:
                    default = ''
                sub_config[name] = default
                sub_config.comments[name] = wrap(description)

            self.config['sed_creation_modules'].comments[module_name] = [
                creation_modules.get_module(module_name, blank=True).comments]

        # Configuration for the analysis method
        self.config['analysis_configuration'] = {}
        self.config.comments['analysis_configuration'] = ["", ""] + wrap(
            "Configuration of the statistical analysis method.")
        module_name = self.config['analysis_method']
        for name, (typ, desc, default) in \
                analysis_modules.get_module(module_name).parameter_list.items():
            if default is None:
                default = ''
            self.config['analysis_configuration'][name] = default
            self.config['analysis_configuration'].comments[name] = wrap(desc)

        self.config.write()

    @property
    def configuration(self):
        """Returns a dictionary for the session configuration.

        Returns
        -------
        configuration['data_file'] : string
            File containing the observations to fit.
        configuration['column_list'] : list of strings
            List of the columns of data_file to use in the fitting.
        configuration['creation_modules'] : list of strings
            List of the modules (in the right order) used to create the SEDs.
        configuration['creation_modules_params'] : list of dictionaries
            Configuration parameters for each module. To each parameter, the
            dictionary associates a list of possible values (possibly only
            one).
        configuration['analysis_method'] : string
            Statistical analysis module used to fit the data.
        configuration['analysis_method_params'] : dictionary
            Parameters for the statistical analysis module. To each parameter
            is associated a list of possible values.
        """
        configuration = {}

        for section in ['data_file', 'column_list', 'creation_modules',
                        'analysis_method']:
            configuration[section] = self.config[section]

        # Parsing the SED modules parameters
        configuration['creation_modules_params'] = []
        for module in self.config['creation_modules']:
            module_params = {}
            for key, value in \
                    self.config['sed_creation_modules'][module].items():
                module_params[key] = evaluate_description(value)
            configuration['creation_modules_params'].append(module_params)

        # Analysis method parameters
        configuration['analysis_method_params'] = \
            self.config['analysis_configuration']

        return configuration

    @property
    def creation_modules_conf_array(self):
        """Return the array of all the possible parameter sets from the
        SED creation modules.

        TODO: Maybe it would be more optimal to create an iterator that would
              iterate over the whole parameter combinations instead of
              creating the array.

        Returns
        -------
        result : array of arrays of dictionaries
            The inner arrays contains the various parameter dictionaries
            for the modules listed in configuration['creation_modules'].

        """

        # First, for each module, we transform the dictionary containing all
        # the possible value for each parameter in a list of dictionaries
        # containing one value for each parameter. We put this list in a list
        # corresponding to the SED modules one.
        tmp_list = [param_dict_combine(dictionary) for dictionary in
                    self.configuration['creation_modules_params']]

        # The we use itertools to create an array of all possible
        # combinations.
        return [x for x in itertools.product(*tmp_list)]
