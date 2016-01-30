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

    def _process(self, configuration):
        """Do the actual analysis

        This method is responsible for the fitting / analysis process
        and must be implemented by each real module.

        Parameters
        ----------
        configuration: dictionary
            Configuration file

        Returns
        -------
        The process results are saved to disk by the analysis module.

        """
        raise NotImplementedError()

    def process(self, configuration):
        """Process with the analysis

        This method is responsible for checking the module parameters before
        doing the actual processing (_process method). If a parameter is not
        given but exists in the parameter_list with a default value, this
        value is used.

        Parameters
        ----------
        configuration: dictionary
            Contents of pcigale.ini in the form of a dictionary

        Returns
        -------
        The process results are saved to disk by the analysis module

        Raises
        ------
        KeyError: when not all the needed parameters are given.

        """
        parameters = configuration['analysis_params']
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
        self._process(configuration)


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


def adjust_data(fluxes, errors, tolerance, lim_flag, default_error=0.1,
                systematic_deviation=0.1):
    """Adjust the fluxes and errors replacing the invalid values by NaN, and
    adding the systematic deviation. The systematic deviation changes the
    errors to: sqrt(errors² + (fluxes*deviation)²)

    Parameters
    ----------
    fluxes: array of floats
        Observed fluxes.
    errors: array of floats
        Observational errors in the same unit as the fluxes.
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
    if len(fluxes) != len(errors):
        raise ValueError("The flux and error arrays must have the same "
                         "length.")

    # We copy the arrays not to modify the original ones.
    fluxes = fluxes.copy()
    errors = errors.copy()

    # We set invalid data to NaN
    mask_invalid = np.where((fluxes <= tolerance) | (errors < -9990.))
    fluxes[mask_invalid] = np.nan
    errors[mask_invalid] = np.nan

    # Replace missing errors by the default ones.
    mask_noerror = np.where((fluxes > tolerance) & ~np.isfinite(errors))
    errors[mask_noerror] = (default_error * fluxes[mask_noerror])

    # Replace upper limits by no data if lim_flag==False
    if not lim_flag:
        mask_limflag = np.where((fluxes > tolerance) & (errors < tolerance))
        fluxes[mask_limflag] = np.nan
        errors[mask_limflag] = np.nan

    # Add the systematic error.
    mask_ok = np.where((fluxes > tolerance) & (errors > tolerance))
    errors[mask_ok] = np.sqrt(errors[mask_ok]**2 +
                             (fluxes[mask_ok]*systematic_deviation)**2)

    return fluxes, errors


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
        if name_err not in obs_table.columns:
            obs_table.add_column(Column(name=name_err,
                                        data=np.full(len(obs_table), np.nan)),
            index=obs_table.colnames.index(name)+1)
        elif name_err not in used_columns:
            obs_table[name_err] = np.full(len(obs_table), np.nan)

        obs_table[name], obs_table[name_err] = adjust_data(obs_table[name],
                                                           obs_table[name_err],
                                                           tolerance,
                                                           lim_flag,
                                                           default_error,
                                                           systematic_deviation)
    return obs_table
