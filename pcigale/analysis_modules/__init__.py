# -*- coding: utf-8 -*-
# Copyright (C) 2014 Laboratoire d'Astrophysique de Marseille, AMU
# Copyright (C) 2012, 2014 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly & Denis Burgarella

from importlib import import_module
import numpy as np
from astropy.table import Column


class AnalysisModule(object):
    """Abstract class, the pCigale analysis modules are based on.
    """

    # parameter_list is a dictionary containing all the parameters
    # used by the module. Each parameter name is associate to a tuple
    # (variable type, description [string], default value). Each module must
    # define its parameter list, unless it does not use any parameter. Using
    # None means that there is no description, unit or default value. If None
    # should be the default value, use the 'None' string instead.
    parameter_list = dict()

    def __init__(self, **kwargs):
        """Instantiate a analysis module

        The module parameters values can be passed as keyword parameters.
        """
        # parameters is a dictionary containing the actual values for each
        # module parameter.
        self.parameters = kwargs

    def _process(self, data_file, column_list, creation_modules,
                 creation_modules_params, parameters):
        """Do the actual analysis

        This method is responsible for the fitting / analysis process
        and must be implemented by each real module.

        Parameters
        ----------
        data_file: string
            Name of the file containing the observations to be fitted.
        column_list: array of strings
            Names of the columns from the data file to use in the analysis.
        creation_modules: array of strings
            Names (in the right order) of the modules to use to build the SED.
        creation_modules_params: array of array of dictionaries
            Array containing all the possible combinations of configurations
            for the creation_modules. Each 'inner' array has the same length as
            the creation_modules array and contains the configuration
            dictionary for the corresponding module.
        parameters: dictionary
            Configuration for the module.

        Returns
        -------
        The process results are saved to disk by the analysis module.

        """
        raise NotImplementedError()

    def process(self, data_file, column_list, creation_modules,
                creation_modules_params, parameters):
        """Process with the analysis

        This method is responsible for checking the module parameters before
        doing the actual processing (_process method). If a parameter is not
        given but exists in the parameter_list with a default value, this
        value is used.

        Parameters
        ----------
        data_file: string
            Name of the file containing the observations to be fitted.
        column_list: array of strings
            Names of the columns from the data file to use in the analysis.
        creation_modules: array of strings
            Names (in the right order) of the modules to use to build the SED.
        creation_modules_params: array of array of dictionaries
            Array containing all the possible combinations of configurations
            for the creation_modules. Each 'inner' array has the same length as
            the creation_modules array and contains the configuration
            dictionary for the corresponding module.
        parameters: dictionary
            Configuration for the module.

        Returns
        -------
        The process results are saved to disk by the analysis module

        Raises
        ------
        KeyError: when not all the needed parameters are given.

        """
        # For parameters that are present on the parameter_list with a default
        # value and that are not in the parameters dictionary, we add them
        # with their default value.
        for key in self.parameter_list:
            if (key not in parameters) and (
                    self.parameter_list[key][2] is not None):
                parameters[key] = self.parameter_list[key][2]

        # If the keys of the parameters dictionary are different from the one
        # of the parameter_list dictionary, we raises a KeyError. That means
        # that a parameter is missing (and has no default value) or that an
        # unexpected one was given.
        if not set(parameters) == set(self.parameter_list):
            missing_parameters = (set(self.parameter_list) - set(parameters))
            unexpected_parameters = (set(parameters) -
                                     set(self.parameter_list))
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

        # We do the actual processing
        self._process(data_file, column_list, creation_modules,
                      creation_modules_params, parameters)


def get_module(module_name):
    """Return the main class of the module provided

    Parameters
    ----------
    module_name: string
        The name of the module we want to get the class.

    Returns
    -------
    module_class: class
    """

    try:
        module = import_module('.' + module_name, 'pcigale.analysis_modules')
        return module.Module()
    except ImportError:
        print('Module ' + module_name + ' does not exists!')
        raise


def adjust_errors(flux, error, tolerance, lim_flag, default_error=0.1,
                  systematic_deviation=0.1):
    """Adjust the errors replacing the 0 values by the default error and
    adding the systematic deviation.

    The systematic deviation change the error to:
    sqrt( error² + (flux * deviation)² )

    Parameters
    ----------
    flux: array of floats
        Fluxes.
    error: array of floats
        Observational error in the same unit as the fluxes.
    tolerance: float
        Tolerance threshold under flux error is considered as 0.
    lim_flag: boolean
        Do we process upper limits (True) or treat them as no-data (False)?
    default_error: float
        Default error factor used when the provided error in under the
        tolerance threshold.
    systematic_deviation: float
        Systematic deviation added to the error.

    Returns
    -------
    error: array of floats
        The corrected errors.

    """
    # The arrays must have the same lengths.
    if len(flux) != len(error):
        raise ValueError("The flux and error arrays must have the same "
                         "length.")

    # We copy the error array not to modify the original one.
    error = np.copy(error)

    # We define:
    # 1) upper limit == (lim_flag==True) and
    #                [(flux > tolerance) and (-9990. < error < tolerance)]
    # 2) no data == (flux < -9990.) and (error < -9990.)
    # Note that the upper-limit test must be performed before the no-data test
    # because if lim_flag is False, we process upper limits as no-data.
    #
    # Replace errors below tolerance by the default one.
    mask_noerror = np.logical_and(flux > tolerance, error < -9990.)
    error[mask_noerror] = (default_error * flux[mask_noerror])

    mask_limflag = np.logical_and.reduce(
                       (flux > tolerance, error >= -9990., error < tolerance))

    # Replace upper limits by no data if lim_flag==False
    if not lim_flag:
        flux[mask_limflag] = -9999.
        error[mask_limflag] = -9999.

    mask_ok = np.logical_and(flux > tolerance, error > tolerance)

    # Add the systematic error.
    error[mask_ok] = np.sqrt(error[mask_ok]**2 +
                             (flux[mask_ok]*systematic_deviation)**2)
    return error


def complete_obs_table(obs_table, used_columns, filter_list, tolerance,
                       lim_flag, default_error=0.1, systematic_deviation=0.1):
    """Complete the observation table

    For each filter:
    * If the corresponding error is not present in the used column list or in
      the table columns, add (or replace) an error column with the default
      error.
    * Adjust the error value.

    Parameters
    ----------
    obs_table: astropy.table.Table
        The observation table.
    used_columns: list of strings
        The list of columns to use in the observation table.
    filter_list: list of strings
        The list of filters used in the analysis.
    tolerance: float
        Tolerance threshold under flux error is considered as 0.
    lim_flag: boolean
        Do we process upper limits (True) or treat them as no-data (False)?
    default_error: float
        Default error factor used when the provided error in under the
        tolerance threshold.
    systematic_deviation: float
        Systematic deviation added to the error.

    Returns
    -------
    obs_table = astropy.table.Table
        The completed observation table

    Raises
    ------
    Exception: When a filter is not present in the observation table.

    """
    # TODO Print or log a warning when an error column is in the used column
    # list but is not present in the observation table.
    for name in filter_list:
        if name not in obs_table.columns:
            raise Exception("The filter <{}> (at least) is not present in "
                                "the observation table.".format(name))

        name_err = name + "_err"
        if name_err not in used_columns:
            if name_err not in obs_table.columns:
                obs_table.add_column(Column(
                    name=name_err,
                    data=np.ones(len(obs_table), dtype=float))*-9999.,
                    index=obs_table.colnames.index(name)+1
                )
            else:
                obs_table[name_err] = np.zeros(len(obs_table))

        obs_table[name_err] = adjust_errors(obs_table[name],
                                            obs_table[name_err],
                                            tolerance,
                                            lim_flag,
                                            default_error,
                                            systematic_deviation)
    return obs_table
