# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013-2014 Yannick Roehlly
# Copyright (C) 2013 Institute of Astronomy
# Copyright (C) 2014 Laboratoire d'Astrophysique de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly, Médéric Boquien & Denis Burgarella

__version__ = "0.1-alpha"

import argparse
from astropy.table import Table
from itertools import product, repeat
from collections import OrderedDict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import os
from pcigale.data import Database
from pcigale.utils import read_table
from pcigale.session.configuration import Configuration
import matplotlib.gridspec as gridspec

# Name of the file containing the best models information
BEST_MODEL_FILE = "best_models.txt"
# Directory where the output files are stored
OUT_DIR = "out/"
# Wavelength limits (restframe) when plotting the best SED.
PLOT_L_MIN = 91.
PLOT_L_MAX = 1e6


def _chi2_worker(obj_name, var_name):
    """Plot the reduced χ² associated with a given analysed variable

    Parameters
    ----------
    obj_name: string
        Name of the object.
    var_name: string
        Name of the analysed variable..

    """
    if os.path.isfile(OUT_DIR + "{}_{}_chi2.fits".format(obj_name, var_name)):
        chi2 = Table.read(OUT_DIR + "{}_{}_chi2.fits".format(obj_name,
                                                             var_name))
        figure = plt.figure()
        ax = figure.add_subplot(111)
        ax.scatter(chi2[var_name], chi2['chi2'], color='k', s=.1)
        ax.set_xlabel(var_name)
        ax.set_ylabel("Reduced $\chi^2$")
        ax.set_ylim(0., )
        ax.minorticks_on()
        figure.suptitle("Reduced $\chi^2$ distribution of {} for {}."
                        .format(var_name, obj_name))
        figure.savefig(OUT_DIR + "{}_{}_chi2.pdf".format(obj_name, var_name))
        plt.close(figure)
    else:
        print("No chi² found for {}. No plot created.".format(obj_name))


def _pdf_worker(obj_name, var_name):
    """Plot the PDF associated with a given analysed variable

    Parameters
    ----------
    obj_name: string
        Name of the object.
    var_name: string
        Name of the analysed variable..

    """
    if os.path.isfile(OUT_DIR + "{}_{}_pdf.fits".format(obj_name, var_name)):
        pdf = Table.read(OUT_DIR + "{}_{}_pdf.fits".format(obj_name, var_name))
        figure = plt.figure()
        ax = figure.add_subplot(111)
        ax.plot(pdf[var_name], pdf['probability density'], color='k')
        ax.set_xlabel(var_name)
        ax.set_ylabel("Probability density")
        ax.minorticks_on()
        figure.suptitle("Probability distribution function of {} for {}"
                        .format(var_name, obj_name))
        figure.savefig(OUT_DIR + "{}_{}_pdf.pdf".format(obj_name, var_name))
        plt.close(figure)
    else:
        print("No PDF found for {}. No plot created.".format(obj_name))


def _sed_worker(obs, mod, filters):
    """Plot the best SED with the associated fluxes in bands

    Parameters
    ----------
    obs: Table row
        Data from the input file regarding one object.
    mod: Table row
        Data from the best model of one object.
    filters: ordered dictionary of Filter objects
        The observed fluxes in each filter.

    """
    if os.path.isfile(OUT_DIR + "{}_best_model.xml".format(obs['id'])):
        sed = Table.read(OUT_DIR + "{}_best_model.xml".format(obs['id']),
                         table_id="Fnu")
        filters_wl = np.array([filt.effective_wavelength
                               for filt in filters.values()])
        obs_fluxes = np.array([obs[filt] for filt in filters.keys()])
        obs_fluxes_err = np.array([obs[filt+'_err']
                                   for filt in filters.keys()])
        mod_fluxes = np.array([mod[filt] for filt in filters.keys()])
        xmin = PLOT_L_MIN * (1. + obs['redshift'])
        xmax = PLOT_L_MAX * (1. + obs['redshift'])
        wsed = np.where((sed['wavelength'] > xmin)&(sed['wavelength'] < xmax))
        figure = plt.figure()
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        if (sed['F_nu'][wsed] > 0.).any():

          ax1 = plt.subplot(gs[0])
          ax2 = plt.subplot(gs[1])
          ax1.loglog(sed['wavelength'][wsed], sed['F_nu'][wsed],
                    label="Model spectrum", color='k', nonposy='clip')
          ax1.set_autoscale_on(False)
          ax1.scatter(filters_wl, mod_fluxes, marker='o', color='r',
                     label="Model fluxes")
          ax1.errorbar(filters_wl, obs_fluxes, yerr=obs_fluxes_err*3, ls='',
                       marker='_', label='Observed fluxes', color='b',
                       capsize=0.)
          mask = np.where(obs_fluxes > 0.)
          ax2.errorbar(filters_wl[mask],
                       (obs_fluxes[mask]-mod_fluxes[mask])/obs_fluxes[mask],
                       yerr=obs_fluxes_err[mask]/obs_fluxes[mask]*3, marker='_',
                       label="(Obs-Mod)/Obs", color='k', capsize=0.)
          ax2.plot([xmin, xmax], [0., 0.], ls='--', color='k')
          ax2.set_xscale('log')
          ax2.minorticks_on()

          figure.subplots_adjust(hspace=0.,wspace=0.)

          ax1.set_xlim(xmin, xmax)
          ax2.set_xlim(xmin, xmax)
          ax2.set_ylim(-1.0, 1.0)
          ax2.set_xlabel("Observed wavelength [nm]")
          ax1.set_ylabel("Flux [mJy]")
          ax2.set_ylabel("Relative residual flux")
          ax1.legend(loc='upper left', fontsize=10)
          ax2.legend(loc='upper right', fontsize=10)
          plt.setp(ax1.get_xticklabels(), visible=False)
          plt.setp(ax1.get_yticklabels()[1], visible=False)
          figure.suptitle("Best model for {} at z = {}. Reduced $\chi^2$={}".
                      format(
                          obs['id'], np.round(obs['redshift'], decimals=3),
                          np.round(mod['reduced_chi_square'], decimals=2))
                            )
          figure.savefig(OUT_DIR + "{}_best_model.pdf".format(obs['id']))
          plt.close(figure)
        else:
          print("No valid best SED found for {}. No plot created.".format(obs['id']))
    else:
        print("No SED found for {}. No plot created.".format(obs['id']))


def chi2(config):
    """Plot the χ² values of analysed variables.
    """
    input_data = read_table(config.configuration['data_file'])
    chi2_vars = (config.configuration['analysis_method_params']
                 ['analysed_variables'])

    with mp.Pool(processes=config.configuration['cores']) as pool:
        items = product(input_data['id'], chi2_vars)
        pool.starmap(_chi2_worker, items)
        pool.close()
        pool.join()


def chi2(config):
    """Plot the χ² values of analysed variables.
    """
    input_data = read_table(config.configuration['data_file'])
    chi2_vars = (config.configuration['analysis_method_params']
                 ['analysed_variables'])

    with mp.Pool(processes=config.configuration['cores']) as pool:
        items = product(input_data['id'], chi2_vars)
        pool.starmap(_chi2_worker, items)
        pool.close()
        pool.join()


def pdf(config):
    """Plot the PDF of analysed variables.
    """
    input_data = read_table(config.configuration['data_file'])
    pdf_vars = (config.configuration['analysis_method_params']
                ['analysed_variables'])

    with mp.Pool(processes=config.configuration['cores']) as pool:
        items = product(input_data['id'], pdf_vars)
        pool.starmap(_pdf_worker, items)
        pool.close()
        pool.join()


def sed(config):
    """Plot the best SED with associated observed and modelled fluxes.
    """
    obs = read_table(config.configuration['data_file'])
    mod = Table.read(OUT_DIR + BEST_MODEL_FILE, format='ascii')
    with Database() as base:
        filters = OrderedDict([(name, base.get_filter(name))
                               for name in config.configuration['column_list']
                               if not name.endswith('_err')])

    with mp.Pool(processes=config.configuration['cores']) as pool:
        pool.starmap(_sed_worker, zip(obs, mod, repeat(filters)))
        pool.close()
        pool.join()


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--conf-file', dest='config_file',
                        help="Alternative configuration file to use.")

    subparsers = parser.add_subparsers(help="List of commands")

    pdf_parser = subparsers.add_parser('pdf', help=pdf.__doc__)
    pdf_parser.set_defaults(parser='pdf')

    chi2_parser = subparsers.add_parser('chi2', help=chi2.__doc__)
    chi2_parser.set_defaults(parser='chi2')

    sed_parser = subparsers.add_parser('sed', help=sed.__doc__)
    sed_parser.set_defaults(parser='sed')

    args = parser.parse_args()

    if args.config_file:
        config = Configuration(args.config_file)
    else:
        config = Configuration()

    if args.parser == 'chi2':
        chi2(config)
    elif args.parser == 'pdf':
        pdf(config)
    elif args.parser == 'sed':
        sed(config)
