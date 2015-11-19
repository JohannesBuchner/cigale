# -*- coding: utf-8 -*-
# Copyright (C) 2015 Laboratoire d'Astrophysique de Marseille
# Copyright (C) 2015 Universidad de Antofagasta
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Médéric Boquien & Denis Burgarella

import argparse
import multiprocessing as mp
import os
import sys

from astropy.table import Table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pkg_resources
from scipy import stats

from pcigale.session.configuration import Configuration

__version__ = "0.1-alpha"

# Name of the file containing the best models information
BEST_MODEL_FILE = "best_models.txt"
MOCK_OUTPUT_FILE = "analysis_mock_results.txt"
# Directory where the output files are stored
OUT_DIR = "out/"


def worker(exact, estimated, param, nologo):
    """Plot the exact and estimated values of a parameter for the mock analysis

    Parameters
    ----------
    exact: Table column
        Exact values of the parameter.
    estimated: Table column
        Estimated values of the parameter.
    param: string
        Name of the parameter
    nologo: boolean
        Do not add the logo when set to true.

    """

    range_exact = np.linspace(np.min(exact), np.max(exact), 100)

    # We compute the linear regression
    if (np.min(exact) < np.max(exact)):
        slope, intercept, r_value, p_value, std_err = stats.linregress(exact,
                                                                       estimated)
    else:
        slope = 0.0
        intercept = 1.0
        r_value = 0.0

    plt.errorbar(exact, estimated, marker='.', label=param, color='k',
                 linestyle='None', capsize=0.)
    plt.plot(range_exact, range_exact, color='r', label='1-to-1')
    plt.plot(range_exact, slope * range_exact + intercept, color='b',
             label='exact-fit $r^2$ = {:.2f}'.format(r_value**2))
    plt.xlabel('Exact')
    plt.ylabel('Estimated')
    plt.title(param)
    plt.legend(loc='best', fancybox=True, framealpha=0.5, numpoints=1)
    plt.minorticks_on()
    if nologo is False:
        image = plt.imread(pkg_resources.resource_filename(__name__,
                                                           "data/CIGALE.png"))
        plt.figimage(image, 510, 55, origin='upper', zorder=10, alpha=1)

    plt.tight_layout()
    plt.savefig(OUT_DIR + 'mock_{}.pdf'.format(param))

    if np.all(exact > 0.) and range_exact[-1]/range_exact[0] > 20.:
        plt.loglog()
        plt.savefig(OUT_DIR + 'mock_log_{}.pdf'.format(param))
    plt.close()


def mock(config, nologo):
    """Plot the comparison of input/output values of analysed variables in
    parallel.

    Parameters
    ----------
    config: configuration object
        Contains the pcigale.ini configuration data
    nologo: boolean
        True when the CIGALE logo should not be drawn

    """

    try:
        exact = Table.read(OUT_DIR + BEST_MODEL_FILE, format='ascii')
    except FileNotFoundError:
        print("Best models file {} not found.".format(OUT_DIR +
                                                      BEST_MODEL_FILE))
        sys.exit(1)

    try:
        estimated = Table.read(OUT_DIR + MOCK_OUTPUT_FILE, format='ascii')
    except FileNotFoundError:
        print("Mock models file {} not found.".format(OUT_DIR +
                                                      MOCK_OUTPUT_FILE))
        sys.exit(1)

    params = config.configuration['analysis_method_params']['analysed_variables']
    arguments = ((exact[param], estimated[param], param, nologo) for param in
                 params)
    with mp.Pool(processes=config.configuration['cores']) as pool:
        pool.starmap(worker, arguments)
        pool.close()
        pool.join()


def main():

    if sys.version_info[:2] >= (3, 4):
        mp.set_start_method('spawn')
    else:
        print("Could not set the multiprocessing start method to spawn. If "
              "you encounter a deadlock, please upgrade to Python≥3.4.")

    parser = argparse.ArgumentParser(description="Create diagnostic plots for "
                                     "CIGALE")

    parser.add_argument('--nologo', action="store_true",
                        help='Do not plot the CIGALE logo')

    args = parser.parse_args()
    config = Configuration()

    mock(config, args.nologo)
