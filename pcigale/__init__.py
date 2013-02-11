# -*- coding: utf-8 -*-
"""
Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""

import argparse
from .session.configuration import Configuration
from .stats.common import get_module as get_stats_module


def init(config):
    "Create a blank configuration file."
    config.create_blank_conf()
    print("The initial configuration file was created. Please complete it "
          "with the data file name and the pcigale modules to use.")


def genconf(config):
    "Generate the full configuration."
    config.generate_conf()
    print("The configuration file has been updated. Please complete the "
          "various module parameters and the data file columns to use in "
          "the analysis.")


def check(config):
    "Check the configuration."
    # TODO : Check if all the parameters that don't have default values are
    # given for each module.
    print ("With this configuration, pcigale must compute {} "
           "SEDs.".format(len(config.sed_modules_conf_array)))


def run(config):
    "Run the analysis."
    data_file = config.configuration['data_file']
    column_list = config.configuration['column_list']
    sed_modules = config.configuration['sed_modules']
    sed_modules_params = config.sed_modules_conf_array
    stat_module = get_stats_module(config.configuration['analysis_method'])
    stat_module_params = config.configuration['analysis_method_params']

    stat_module.process(data_file, column_list, sed_modules,
                        sed_modules_params, stat_module_params)


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--conf-file', dest='config_file',
                        help="Alternative configuration file to use.")

    subparsers = parser.add_subparsers(help="List of commands")

    init_parser = subparsers.add_parser('init', help=init.__doc__)
    init_parser.set_defaults(parser='init')

    genconf_parser = subparsers.add_parser('genconf', help=genconf.__doc__)
    genconf_parser.set_defaults(parser='genconf')

    check_parser = subparsers.add_parser('check', help=check.__doc__)
    check_parser.set_defaults(parser='check')

    run_parser = subparsers.add_parser('run', help=run.__doc__)
    run_parser.set_defaults(parser='run')

    args = parser.parse_args()

    if args.config_file:
        config = Configuration(args.config_file)
    else:
        config = Configuration()

    if args.parser == 'init':
        init(config)
    elif args.parser == 'genconf':
        genconf(config)
    elif args.parser == 'check':
        check(config)
    elif args.parser == 'run':
        run(config)
