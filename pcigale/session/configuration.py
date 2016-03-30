# -*- coding: utf-8 -*-
# Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

import pkgutil
from collections import Iterable, OrderedDict
import multiprocessing as mp
from textwrap import wrap

import configobj
from glob import glob  # To allow the use of glob() in "eval..."
import pkg_resources
import numpy as np
import validate

from ..handlers.parameters_handler import ParametersHandler
from ..data import Database
from ..utils import read_table
from .. import sed_modules
from .. import analysis_modules
from ..warehouse import SedWarehouse
from . import validation


# Limit the redshift to this number of decimals
REDSHIFT_DECIMALS = 2


class Configuration(object):
    """This class manages the configuration of pcigale.
    """

    def __init__(self, filename="pcigale.ini"):
        """Initialise a pcigale configuration.

        Parameters
        ----------
        filename: string
            Name of the configuration file (pcigale.conf by default).

        """
        self.spec = configobj.ConfigObj(filename+'.spec',
                                        write_empty_values=True,
                                        indent_type='  ',
                                        encoding='UTF8',
                                        list_values=False,
                                        _inspec=True)
        self.config = configobj.ConfigObj(filename,
                                          write_empty_values=True,
                                          indent_type='  ',
                                          encoding='UTF8',
                                          configspec=self.spec)

        # We validate the configuration so that the variables are converted to
        # the expected that. We do not handle errors at the point but only when
        # we actually return the configuration file from the property() method.
        self.config.validate(validate.Validator(validation.functions))

    def create_blank_conf(self):
        """Create the initial configuration file

        Write the initial pcigale configuration file where the user can state
        which data file to use, which modules to use for the SED creation, as
        well as the method selected for statistical analysis.

        """

        self.config['data_file'] = ""
        self.config.comments['data_file'] = wrap(
            "File containing the input data. The columns are 'id' (name of the"
            " object), 'redshift' (if 0 the distance is assumed to be 10 pc), "
            "the filter names for the fluxes, and the filter names with the "
            "'_err' suffix for the uncertainties. The fluxes and the "
            "uncertainties must be in mJy. This file is optional to generate "
            "the configuration file, in particular for the savefluxes module.")
        self.spec['data_file'] = "string"

        self.config['parameters_file'] = ""
        self.config.comments['parameters_file'] = [""] + wrap(
            "Optional file containing the list of physical parameters. Each "
            "column must be in the form module_name.parameter_name, with each "
            "line behind a different model. The columns must be in the order "
            "the modules will be called. The redshift column must be the last "
            "one. Finally, if this parameters is not left empty, cigale will "
            "not interpret the configuration parameters given in pcigale.ini. "
            "They will be given only for information.")
        self.spec['parameters_file'] = "string()"

        self.config['sed_modules'] = []
        self.config.comments['sed_modules'] = ([""] +
            ["Order of the modules use for SED creation. Available modules:"] +
            ["SFH: sfh2exp, sfhdelayed, sfhfromfile, sfhperiodic"] +
            ["SSP: bc03, m2005"] +
            ["Nebular emission: nebular"] +
            ["Dust attenuation: dustatt_calzleit, dustatt_powerlaw, "
             "dustatt_2powerlaws"] +
            ["Dust emission: casey2012, dale2014, dl2007, dl2014"] +
            ["AGN: dale2014, fritz2006"] +
            ["Radio: radio"] +
            ["Redshift: redshifting (mandatory!)"])
        self.spec['sed_modules'] = "cigale_string_list()"

        self.config['analysis_method'] = ""
        self.config.comments['analysis_method'] = [""] + wrap(
            "Method used for statistical analysis. Available methods: "
            "pdf_analysis, savefluxes.")
        self.spec['analysis_method'] = "string()"

        self.config['cores'] = ""
        self.config.comments['cores'] = [""] + wrap(
            "Number of CPU cores available. This computer has {} cores."
            .format(mp.cpu_count()))
        self.spec['cores'] = "integer(min=1)"

        self.config.write()
        self.spec.write()

    def generate_conf(self):
        """Generate the full configuration file

        Reads the user entries in the initial configuration file and add the
        configuration options of all selected modules as well as the filter
        selection based on the filters identified in the data table file.

        """

        # Getting the list of the filters available in pcigale database
        with Database() as base:
            filter_list = base.get_filter_names()

        if self.config['data_file'] != '':
            obs_table = read_table(self.config['data_file'])

            # Check that the id and redshift columns are present in the input
            # file
            if 'id' not in obs_table.columns:
                raise Exception("Column id not present in input file")
            if 'redshift' not in obs_table.columns:
                raise Exception("Column redshift not present in input file")

            # Finding the known filters in the data table
            bands = []
            for band in obs_table.columns:
                filter_name = band[:-4] if band.endswith('_err') else band
                if filter_name in filter_list:
                    bands.append(band)

            # Check that we don't have an band error without the associated
            # band
            for band in bands:
                if band.endswith('_err') and (band[:-4] not in bands):
                    raise Exception("The observation table as a {} column "
                                    "but no {} column.".format(band,
                                                               band[:-4]))

            self.config['bands'] = bands
        else:
            self.config['bands'] = ''
        self.config.comments['bands'] = [""] + wrap("Bands to consider. To "
            "consider uncertainties too, the name of the band must be "
            "indicated with the _err suffix. For instance: FUV, FUV_err.")
        self.spec['bands'] = "cigale_string_list()"

        # SED creation modules configurations. For each module, we generate
        # the configuration section from its parameter list.
        self.config['sed_modules_params'] = {}
        self.config.comments['sed_modules_params'] = ["", ""] + wrap(
            "Configuration of the SED creation modules.")
        self.spec['sed_modules_params'] = {}

        for module_name in self.config['sed_modules']:
            self.config['sed_modules_params'][module_name] = {}
            self.spec['sed_modules_params'][module_name] = {}
            sub_config = self.config['sed_modules_params'][module_name]
            sub_spec = self.spec['sed_modules_params'][module_name]

            for name, (typ, description, default) in \
                    sed_modules.get_module(
                        module_name,
                        blank=True).parameter_list.items():
                if default is None:
                    default = ''
                sub_config[name] = default
                sub_config.comments[name] = wrap(description)
                sub_spec[name] = typ
            self.config['sed_modules_params'].comments[module_name] = [
                sed_modules.get_module(module_name, blank=True).comments]

        self.check_modules()

        # Configuration for the analysis method
        self.config['analysis_params'] = {}
        self.config.comments['analysis_params'] = ["", ""] + wrap(
            "Configuration of the statistical analysis method.")
        self.spec['analysis_params'] = {}

        module_name = self.config['analysis_method']
        for name, (typ, desc, default) in \
                analysis_modules.get_module(module_name).parameter_list.items():
            if default is None:
                default = ''
            self.config['analysis_params'][name] = default
            self.config['analysis_params'].comments[name] = wrap(desc)
            self.spec['analysis_params'][name] = typ

        self.config.write()
        self.spec.write()

    @property
    def configuration(self):
        """Returns a dictionary for the session configuration if it is valid.
        Otherwise, print the erroneous keys.

        Returns
        -------
        configuration: dictionary
            Dictionary containing the information provided in pcigale.ini.
        """
        self.complete_redshifts()
        self.complete_analysed_parameters()

        vdt = validate.Validator(validation.functions)
        validity = self.config.validate(vdt, preserve_errors=True)

        if validity is not True:
            print("The following issues have been found in pcigale.ini:")
            for module, param, message in configobj.flatten_errors(self.config,
                                                                   validity):
                if len(module) > 0:
                    print("Module {}, parameter {}: {}".format('/'.join(module),
                                                               param, message))
                else:
                    print("Parameter {}: {}".format(param, message))
            print("Run the same command after having fixed pcigale.ini.")

            return None

        return self.config.dict()

    def check_modules(self):
        """Make a basic check to ensure that some required modules are present.
        Otherwise we emit a warning so the user knows their list of modules is
        suspicious. We do not emit an exception as they may be using an
        unofficial module that is not in our list
        """

        modules = OrderedDict((('SFH', ['sfh2exp', 'sfhdelayed', 'sfhfromfile',
                                        'sfhperiodic']),
                               ('SSP', ['bc03', 'm2005']),
                               ('nebular', ['nebular']),
                               ('dust attenuation', ['dustatt_calzleit',
                                                     'dustatt_powerlaw',
                                                     'dustatt_2powerlaws']),
                               ('dust emission', ['casey2012', 'dale2014',
                                                  'dl2007', 'dl2014']),
                               ('AGN', ['dale2014', 'fritz2006']),
                               ('radio', ['radio']),
                               ('redshift', ['redshifting'])))

        comments = {'SFH': "ERROR! Choosing one SFH module is mandatory.",
                    'SSP': "ERROR! Choosing one SSP module is mandatory.",
                    'nebular': "WARNING! Choosing the nebular module is "
                               "recommended. Without it the Lyman continuum "
                               "is left untouched.",
                    'dust attenuation': "No dust attenuation module found.",
                    'dust emission': "No dust emission module found.",
                    'AGN': "No AGN module found.",
                    'radio': "No radio module found.",
                    'redshift': "ERROR! No redshifting module found."}

        for module in modules:
            if all([user_module not in modules[module] for user_module in
                    self.config['sed_modules']]):
                print("{} Options are: {}.".
                      format(comments[module], ', '.join(modules[module])))

    def complete_redshifts(self):
        """Complete the configuration when the redshifts are missing from the
        configuration file and must be extracted from the input flux file.
        """

        z_mod = self.config['sed_modules_params']['redshifting']['redshift']
        if type(z_mod) is str and not z_mod:
            if self.config['data_file']:
                obs_table = read_table(self.config['data_file'])
                z = list(np.unique(np.around(obs_table['redshift'],
                                        decimals=REDSHIFT_DECIMALS)))
                self.config['sed_modules_params']['redshifting']['redshift'] = z
            else:
                raise Exception("No flux file and no redshift indicated. "
                                "The spectra cannot be computed. Aborting.")

    def complete_analysed_parameters(self):
        """Complete the configuration when the variables are missing from the
        configuration file and must be extract from a dummy run."""
        if not self.config['analysis_params']['variables']:
            warehouse = SedWarehouse()
            params = ParametersHandler(self.config.dict())
            sed = warehouse.get_sed(params.modules, params.from_index(0))
            info = list(sed.info.keys())
            info.sort()
            self.config['analysis_params']['variables'] = info
