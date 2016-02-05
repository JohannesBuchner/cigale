# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013-2014 Yannick Roehlly
# Copyright (C) 2013 Institute of Astronomy
# Copyright (C) 2014 Laboratoire d'Astrophysique de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly, Médéric Boquien & Denis Burgarella

import argparse
from itertools import product, repeat
from collections import OrderedDict
import sys

from astropy.table import Table
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
import os
import pkg_resources
from scipy.constants import c
from scipy import stats
from pcigale.data import Database
from pcigale.utils import read_table
from pcigale.session.configuration import Configuration
import matplotlib.gridspec as gridspec

__version__ = "0.1-alpha"


# Name of the file containing the best models information
BEST_RESULTS = "results.fits"
MOCK_RESULTS = "results_mock.fits"
# Directory where the output files are stored
OUT_DIR = "out/"
# Wavelength limits (restframe) when plotting the best SED.
PLOT_L_MIN = 0.1
PLOT_L_MAX = 5e5


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


def _sed_worker(obs, mod, filters, sed_type, nologo):
    """Plot the best SED with the associated fluxes in bands

    Parameters
    ----------
    obs: Table row
        Data from the input file regarding one object.
    mod: Table row
        Data from the best model of one object.
    filters: ordered dictionary of Filter objects
        The observed fluxes in each filter.
    sed_type: string
        Type of SED to plot. It can either be "mJy" (flux in mJy and observed
        frame) or "lum" (luminosity in W and rest frame)
    nologo: boolean
        Do not add the logo when set to true.

    """

    if os.path.isfile(OUT_DIR + "{}_best_model.xml".format(obs['id'])):

        sed = Table.read(OUT_DIR + "{}_best_model.xml".format(obs['id']),
                         table_id="Llambda")

        filters_wl = np.array([filt.effective_wavelength
                               for filt in filters.values()])
        wavelength_spec = sed['wavelength']
        obs_fluxes = np.array([obs[filt] for filt in filters.keys()])
        obs_fluxes_err = np.array([obs[filt+'_err']
                                   for filt in filters.keys()])
        mod_fluxes = np.array([mod["best."+filt] for filt in filters.keys()])
        z = np.around(obs['redshift'], decimals=2)
        DL = mod['best.universe.luminosity_distance']

        if sed_type == 'lum':
            xmin = PLOT_L_MIN
            xmax = PLOT_L_MAX

            k_corr_SED = 1e-29 * (4.*np.pi*DL*DL) * c / (filters_wl*1e-9)
            obs_fluxes *= k_corr_SED
            obs_fluxes_err *= k_corr_SED
            mod_fluxes *= k_corr_SED

            for cname in sed.colnames[1:]:
                sed[cname] *= wavelength_spec

            filters_wl /= 1. + z
            wavelength_spec /= 1. + z
        elif sed_type == 'mJy':
            xmin = PLOT_L_MIN * (1. + z)
            xmax = PLOT_L_MAX * (1. + z)
            k_corr_SED = 1.

            for cname in sed.colnames[1:]:
                sed[cname] *= (wavelength_spec * 1e29 /
                               (c / (wavelength_spec * 1e-9)) /
                               (4. * np.pi * DL * DL))
        else:
            print("Unknown plot type")

        filters_wl /= 1000.
        wavelength_spec /= 1000.

        wsed = np.where((wavelength_spec > xmin) & (wavelength_spec < xmax))
        figure = plt.figure()
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
        if (sed.columns[1][wsed] > 0.).any():
            ax1 = plt.subplot(gs[0])
            ax2 = plt.subplot(gs[1])

            # Stellar emission
            if 'nebular.absorption_young' in sed.columns:
                ax1.loglog(wavelength_spec[wsed],
                           (sed['stellar.young'][wsed] +
                            sed['attenuation.stellar.young'][wsed] +
                            sed['nebular.absorption_young'][wsed] +
                            sed['stellar.old'][wsed] +
                            sed['attenuation.stellar.old'][wsed] +
                            sed['nebular.absorption_old'][wsed]),
                           label="Stellar attenuated ", color='orange',
                           marker=None, nonposy='clip', linestyle='-',
                           linewidth=0.5)
            else:
                ax1.loglog(wavelength_spec[wsed],
                           (sed['stellar.young'][wsed] +
                            sed['attenuation.stellar.young'][wsed] +
                            sed['stellar.old'][wsed] +
                            sed['attenuation.stellar.old'][wsed]),
                           label="Stellar attenuated ", color='orange',
                           marker=None, nonposy='clip', linestyle='-',
                           linewidth=0.5)
            ax1.loglog(wavelength_spec[wsed],
                       (sed['stellar.old'][wsed] +
                        sed['stellar.young'][wsed]),
                       label="Stellar unattenuated", color='b', marker=None,
                       nonposy='clip', linestyle='--', linewidth=0.5)
            # Nebular emission
            if 'nebular.lines_young' in sed.columns:
                ax1.loglog(wavelength_spec[wsed],
                           (sed['nebular.lines_young'][wsed] +
                            sed['nebular.lines_old'][wsed] +
                            sed['nebular.continuum_young'][wsed] +
                            sed['nebular.continuum_old'][wsed] +
                            sed['attenuation.nebular.lines_young'][wsed] +
                            sed['attenuation.nebular.lines_old'][wsed] +
                            sed['attenuation.nebular.continuum_young'][wsed] +
                            sed['attenuation.nebular.continuum_old'][wsed]),
                           label="Nebular emission", color='y', marker=None,
                           nonposy='clip', linewidth=.5)
            # Dust emission Draine & Li
            if 'dust.Umin_Umin' in sed.columns:
                ax1.loglog(wavelength_spec[wsed],
                           (sed['dust.Umin_Umin'][wsed] +
                            sed['dust.Umin_Umax'][wsed]),
                           label="Dust emission", color='r', marker=None,
                           nonposy='clip', linestyle='-', linewidth=0.5)
            # Dust emission Dale
            if 'dust' in sed.columns:
                ax1.loglog(wavelength_spec[wsed], sed['dust'][wsed],
                           label="Dust emission", color='r', marker=None,
                           nonposy='clip', linestyle='-', linewidth=0.5)
            # AGN emission Fritz
            if 'agn.fritz2006_therm' in sed.columns:
                ax1.loglog(wavelength_spec[wsed],
                           (sed['agn.fritz2006_therm'][wsed] +
                            sed['agn.fritz2006_scatt'][wsed] +
                            sed['agn.fritz2006_agn'][wsed]),
                           label="AGN emission", color='g', marker=None,
                           nonposy='clip', linestyle='-', linewidth=0.5)
            # Radio emission
            if 'radio_nonthermal' in sed.columns:
                ax1.loglog(wavelength_spec[wsed],
                           sed['radio_nonthermal'][wsed],
                           label="Radio nonthermal", color='brown',
                           marker=None, nonposy='clip', linestyle='-',
                           linewidth=0.5)

            ax1.loglog(wavelength_spec[wsed], sed['L_lambda_total'][wsed],
                       label="Model spectrum", color='k', nonposy='clip',
                       linestyle='-', linewidth=1.5)

            ax1.set_autoscale_on(False)
            ax1.scatter(filters_wl, mod_fluxes, marker='o', color='r', s=8,
                        zorder=3, label="Model fluxes")
            mask_ok = np.logical_and(obs_fluxes > 0., obs_fluxes_err > 0.)
            ax1.errorbar(filters_wl[mask_ok], obs_fluxes[mask_ok],
                         yerr=obs_fluxes_err[mask_ok]*3, ls='', marker='s',
                         label='Observed fluxes', markerfacecolor='None',
                         markersize=6, markeredgecolor='b', capsize=0.)
            mask_uplim = np.logical_and(np.logical_and(obs_fluxes > 0.,
                                                       obs_fluxes_err < 0.),
                                        obs_fluxes_err > -9990. * k_corr_SED)
            if not mask_uplim.any() == False:
                ax1.errorbar(filters_wl[mask_uplim], obs_fluxes[mask_uplim],
                             yerr=obs_fluxes_err[mask_uplim]*3, ls='',
                             marker='v', label='Observed upper limits',
                             markerfacecolor='None', markersize=6,
                             markeredgecolor='g', capsize=0.)
            mask_noerr = np.logical_and(obs_fluxes > 0.,
                                        obs_fluxes_err < -9990. * k_corr_SED)
            if not mask_noerr.any() == False:
                ax1.errorbar(filters_wl[mask_noerr], obs_fluxes[mask_noerr],
                             ls='', marker='s', markerfacecolor='None',
                             markersize=6, markeredgecolor='r',
                             label='Observed fluxes, no errors', capsize=0.)
            mask = np.where(obs_fluxes > 0.)
            ax2.errorbar(filters_wl[mask],
                         (obs_fluxes[mask]-mod_fluxes[mask])/obs_fluxes[mask],
                         yerr=obs_fluxes_err[mask]/obs_fluxes[mask]*3,
                         marker='_', label="(Obs-Mod)/Obs", color='k',
                         capsize=0.)
            ax2.plot([xmin, xmax], [0., 0.], ls='--', color='k')
            ax2.set_xscale('log')
            ax2.minorticks_on()

            figure.subplots_adjust(hspace=0., wspace=0.)

            ax1.set_xlim(xmin, xmax)
            ymin = min(np.min(obs_fluxes[mask_ok]),
                       np.min(mod_fluxes[mask_ok]))
            if not mask_uplim.any() == False:
                ymax = max(max(np.max(obs_fluxes[mask_ok]),
                               np.max(obs_fluxes[mask_uplim])),
                           max(np.max(mod_fluxes[mask_ok]),
                               np.max(mod_fluxes[mask_uplim])))
            else:
                ymax = max(np.max(obs_fluxes[mask_ok]),
                           np.max(mod_fluxes[mask_ok]))
            ax1.set_ylim(1e-1*ymin, 1e1*ymax)
            ax2.set_xlim(xmin, xmax)
            ax2.set_ylim(-1.0, 1.0)
            if sed_type == 'lum':
                ax2.set_xlabel("Rest-frame wavelength [$\mu$m]")
                ax1.set_ylabel("Luminosity [W]")
                ax2.set_ylabel("Relative residual luminosity")
            else:
                ax2.set_xlabel("Observed wavelength [$\mu$m]")
                ax1.set_ylabel("Flux [mJy]")
                ax2.set_ylabel("Relative residual flux")
            ax1.legend(fontsize=6, loc='best', fancybox=True, framealpha=0.5)
            ax2.legend(fontsize=6, loc='best', fancybox=True, framealpha=0.5)
            plt.setp(ax1.get_xticklabels(), visible=False)
            plt.setp(ax1.get_yticklabels()[1], visible=False)
            figure.suptitle("Best model for {} at z = {}. Reduced $\chi^2$={}".
                            format(obs['id'], np.round(obs['redshift'],
                                   decimals=3),
                                   np.round(mod['best.reduced_chi_square'],
                                            decimals=2)))
            if nologo is False:
                image = plt.imread(pkg_resources.resource_filename(__name__,
                                   "data/CIGALE.png"))
                figure.figimage(image, 75, 330, origin='upper', zorder=10,
                                alpha=1)
            figure.savefig(OUT_DIR + "{}_best_model.pdf".format(obs['id']))
            plt.close(figure)
        else:
            print("No valid best SED found for {}. No plot created.".
                  format(obs['id']))
    else:
        print("No SED found for {}. No plot created.".format(obs['id']))


def _mock_worker(exact, estimated, param, nologo):
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

    plt.close()


def chi2(config):
    """Plot the χ² values of analysed variables.
    """
    input_data = read_table(config.configuration['data_file'])
    chi2_vars = config.configuration['analysis_params']['variables']

    with mp.Pool(processes=config.configuration['cores']) as pool:
        items = product(input_data['id'], chi2_vars)
        pool.starmap(_chi2_worker, items)
        pool.close()
        pool.join()


def pdf(config):
    """Plot the PDF of analysed variables.
    """
    input_data = read_table(config.configuration['data_file'])
    pdf_vars = config.configuration['analysis_params']['variables']

    with mp.Pool(processes=config.configuration['cores']) as pool:
        items = product(input_data['id'], pdf_vars)
        pool.starmap(_pdf_worker, items)
        pool.close()
        pool.join()


def sed(config, sed_type, nologo):
    """Plot the best SED with associated observed and modelled fluxes.
    """
    obs = read_table(config.configuration['data_file'])
    mod = Table.read(OUT_DIR + BEST_RESULTS)

    with Database() as base:
        filters = OrderedDict([(name, base.get_filter(name))
                               for name in config.configuration['bands']
                               if not name.endswith('_err')])

    with mp.Pool(processes=config.configuration['cores']) as pool:
        pool.starmap(_sed_worker, zip(obs, mod, repeat(filters),
                                      repeat(sed_type), repeat(nologo)))
        pool.close()
        pool.join()


def mock(config, nologo):
    """Plot the comparison of input/output values of analysed variables.
    """

    try:
        exact = Table.read(OUT_DIR + BEST_RESULTS)
    except FileNotFoundError:
        print("Best models file {} not found.".format(OUT_DIR + BEST_RESULTS))
        sys.exit(1)

    try:
        estimated = Table.read(OUT_DIR + MOCK_RESULTS)
    except FileNotFoundError:
        print("Mock models file {} not found.".format(OUT_DIR + MOCK_RESULTS))
        sys.exit(1)

    params = config.configuration['analysis_params']['variables']

    for param in params:
        if param.endswith('_log'):
            param = "best."+param
            exact[param] = np.log10(exact[param[:-4]])

    arguments = ((exact["best."+param], estimated["bayes."+param], param, nologo)
                 for param in params)

    with mp.Pool(processes=config.configuration['cores']) as pool:
        pool.starmap(_mock_worker, arguments)
        pool.close()
        pool.join()


def main():

    if sys.version_info[:2] >= (3, 4):
        mp.set_start_method('spawn')
    else:
        print("Could not set the multiprocessing start method to spawn. If "
              "you encounter a deadlock, please upgrade to Python≥3.4.")

    parser = argparse.ArgumentParser()

    parser.add_argument('-c', '--conf-file', dest='config_file',
                        help="Alternative configuration file to use.")

    subparsers = parser.add_subparsers(help="List of commands")

    pdf_parser = subparsers.add_parser('pdf', help=pdf.__doc__)
    pdf_parser.set_defaults(parser='pdf')

    chi2_parser = subparsers.add_parser('chi2', help=chi2.__doc__)
    chi2_parser.set_defaults(parser='chi2')

    sed_parser = subparsers.add_parser('sed', help=sed.__doc__)
    sed_parser.add_argument('--type', default='mJy')
    sed_parser.add_argument('--nologo', action="store_true")
    sed_parser.set_defaults(parser='sed')

    sed_parser = subparsers.add_parser('mock', help=mock.__doc__)
    sed_parser.add_argument('--nologo', action="store_true")
    sed_parser.set_defaults(parser='mock')

    args = parser.parse_args()

    if args.config_file:
        config = Configuration(args.config_file)
    else:
        config = Configuration()

    if len(sys.argv) == 1:
        parser.print_usage()
    else:
        if args.parser == 'chi2':
            chi2(config)
        elif args.parser == 'pdf':
            pdf(config)
        elif args.parser == 'sed':
            sed(config, args.type, args.nologo)
        elif args.parser == 'mock':
            mock(config, args.nologo)
