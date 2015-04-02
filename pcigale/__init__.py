# -*- coding: utf-8 -*-
# Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

__version__ = "0.1-alpha"

import argparse
import multiprocessing as mp
import sys

from .session.configuration import Configuration
from .analysis_modules import get_module as get_analysis_module
from .analysis_modules.utils import ParametersHandler


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
    # TODO: Check if all the parameters that don't have default values are
    # given for each module.
    print("With this configuration, pcigale must compute {} "
          "SEDs.".format(ParametersHandler(
                             config.configuration['creation_modules'],
                             config.configuration['creation_modules_params']
                             ).size))


def run(config):
    "Run the analysis."
    data_file = config.configuration['data_file']
    column_list = config.configuration['column_list']
    creation_modules = config.configuration['creation_modules']
    creation_modules_params = config.configuration['creation_modules_params']
    analysis_module = get_analysis_module(config.configuration[
        'analysis_method'])
    analysis_module_params = config.configuration['analysis_method_params']
    cores = config.configuration['cores']

    analysis_module.process(data_file, column_list, creation_modules,
                            creation_modules_params, analysis_module_params,
                            cores)


def main():
    # We set the sub processes start method to spawn because it solves
    # deadlocks when a library cannot handle being used on two sides of a
    # forked process. This happens on modern Macs with the Accelerate library
    # for instance. Unfortunately this only comes with python≥3.4. People using
    # older versions should upgrade if they encounter deadlocks.
    if sys.version_info[:2] >= (3, 4):
        mp.set_start_method('spawn')
    else:
        print("Could not set the multiprocessing start method to spawn. If "
              "you encounter a deadlock, please upgrade to Python≥3.4.")

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
