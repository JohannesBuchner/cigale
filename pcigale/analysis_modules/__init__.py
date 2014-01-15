# -*- coding: utf-8 -*-
# Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

from importlib import import_module
from collections import OrderedDict


class AnalysisModule(object):
    """Abstract class, the pCigale analysis modules are based on.
    """

    # parameter_list is a ordered dictionary containing all the parameters
    # used by the module. Each parameter name is associate to a tuple
    # (variable type, description [string], default value). Each module must
    # define its parameter list, unless it does not use any parameter. Using
    # None means that there is no description, unit or default value. If None
    # should be the default value, use the 'None' string instead.
    parameter_list = OrderedDict()

    def __init__(self, **kwargs):
        """Instantiate a analysis module

        The module parameters values can be passed as keyword parameters.
        """
        # parameters is a dictionary containing the actual values for each
        # module parameter.
        self.parameters = kwargs

    def _process(self, data_file, column_list, creation_modules,
                 creation_modules_params, redshift_module,
                 redshift_configuration, parameters):
        """Do the actual analysis

        This method is responsible for the fitting / analysis process
        and must be implemented by each real module.

        Parameters
        ----------
        data_file : string
            Name of the file containing the observations to be fitted.
        column_list : array of strings
            Names of the columns from the data file to use in the analysis.
        creation_modules : array of strings
            Names (in the right order) of the modules to use to build the SED.
        creation_modules_params : array of array of dictionaries
            Array containing all the possible combinations of configurations
            for the creation_modules. Each 'inner' array has the same length as
            the creation_modules array and contains the configuration dictionary
            for the corresponding module.
        redshift_module : string
            Name of the module used to redshift the SED.
        redshift_configuration : dictionary
            Configuration dictionary for the module used to redshift the SED.
        parameters : dictionary
            Configuration for the module.

        Returns
        -------
        The process results are saved to disk by the analysis module.

        """
        raise NotImplementedError()

    def process(self, data_file, column_list, creation_modules,
                creation_modules_params, redshift_module,
                redshift_configuration, parameters):
        """Process with the analysis

        This method is responsible for checking the module parameters before
        doing the actual processing (_process method). If a parameter is not
        given but exists in the parameter_list with a default value, this
        value is used.

        Parameters
        ----------
        data_file : string
            Name of the file containing the observations to be fitted.
        column_list : array of strings
            Names of the columns from the data file to use in the analysis.
        creation_modules : array of strings
            Names (in the right order) of the modules to use to build the SED.
        creation_modules_params : array of array of dictionaries
            Array containing all the possible combinations of configurations
            for the creation_modules. Each 'inner' array has the same length as
            the creation_modules array and contains the configuration
            dictionary for the corresponding module.
        redshift_module : string
            Name of the module used to redshift the SED.
        redshift_configuration : dictionary
            Configuration dictionary for the module used to redshift the SED.
        parameters : dictionary
            Configuration for the module.

        Returns
        -------
        The process results are saved to disk by the analysis module

        Raises
        ------
        KeyError : when not all the needed parameters are given.

        """
        # For parameters that are present on the parameter_list with a default
        # value and that are not in the parameters dictionary, we add them
        # with their default value.
        for key in self.parameter_list:
            if (not key in parameters) and (
                    self.parameter_list[key][2] is not None):
                parameters[key] = self.parameter_list[key][2]

        # If the keys of the parameters dictionary are different from the one
        # of the parameter_list dictionary, we raises a KeyError. That means
        # that a parameter is missing (and has no default value) or that an
        # unexpected one was given.
        if not set(parameters) == set(self.parameter_list):
            missing_parameters = (set(self.parameter_list) - set(parameters))
            unexpected_parameters = (set(parameters) - set(self.parameter_list))
            message = ""
            if missing_parameters:
                message += ("Missing parameters: " +
                            ", ".join(missing_parameters) +
                            ".")
            if unexpected_parameters:
                message += ("Unexpected parameters: " +
                            ", ".join(unexpected_parameters) +
                            ".")
            raise KeyError("The parameters passed are different from the "
                           "expected one." + message)

        #We do the actual processing
        self._process(data_file, column_list, creation_modules,
                      creation_modules_params, redshift_module,
                      redshift_configuration, parameters)


def get_module(module_name):
    """Return the main class of the module provided

    Parameters
    ----------
    module_name : string
        The name of the module we want to get the class.

    Returns
    -------
    module_class : class
    """

    try:
        module = import_module('.' + module_name, 'pcigale.analysis_modules')
        return module.Module()
    except ImportError:
        print('Module ' + module_name + ' does not exists!')
        raise
