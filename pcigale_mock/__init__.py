# -*- coding: utf-8 -*-
# Copyright (C) 2015 Laboratoire d'Astrophysique de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Denis Burgarella

import argparse
from itertools import product, repeat
from collections import OrderedDict

from astropy.table import Table, join
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages

#matplotlib.use('Agg')
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
BEST_MODEL_FILE = "best_models.txt"
MOCK_OUTPUT_FILE = "analysis_mock_results.txt"
# Directory where the output files are stored
OUT_DIR = "out/"

def worker(best_model, mock_output, mock_params, NRows, NCols, nologo):
    """Plot the diagnostic diagrammes from the mock analysis

    Parameters
    ----------
    best_model: Table row
        Data from the input best_model file.
    mock_output: Table row
        Data from the out analysis of the mock catalogue
    mock_params: Table row
        Name of the parameter analysed that we want to plot
    NRows: Integer
        How many rows per page?
    NCols: Integer
        How many columns per page?
    nologo: boolean
        Do not add the logo when set to true.

    """

    id_best = best_model['observation_id']
    id_mock = mock_output['observation_id']

    if (len(id_best)>0 and len(id_mock)>0):
        joined_table = join(left=best_model,right=mock_output,
                   keys='observation_id',
                   table_names=['Best', 'Mock'],
                   uniq_col_name='{table_name}_{col_name}')
    else:
        print("*** WARNING ***: No object found in {} or {}."
              .format(OUT_DIR + BEST_MODEL_FILE, OUT_DIR + MOCK_OUTPUT_FILE))

    NParams = len(mock_params)
    if NParams <=0:
        print("*** WARNING ***: No analysed parameter found ")

    NPages = int(np.ceil(NParams / float(NRows*NCols)))

    #print(NRows, NCols, NPages, NParams)

    with PdfPages(OUT_DIR + 'mock.pdf') as pdf:

     for Page in range(NPages):
       fig = plt.figure()
       plt.suptitle('Page '+str(Page))
       gs = gridspec.GridSpec(NRows, NCols)
       minParam = Page*NRows*NCols
       maxParam = Page*NRows*NCols+NRows*NCols
       for Param in enumerate(mock_params[minParam: maxParam]):

          # We exclude outliers
          mean_x = np.mean(best_model[Param[1]])
          std_x = np.std(best_model[Param[1]])
          if std_x > 0.:
             mask = np.logical_and(mock_output[Param[1]] > mean_x-5.*std_x,
                                   mock_output[Param[1]] < mean_x+5.*std_x
                                   )
             PercentExcluded = 100*(1.-np.sum(mask)/len(best_model[Param[1]]))
          else:
             mask = True
             PercentExcluded = 0.0

          # We compute the linear regression (excludind outliers)
          if (np.min(best_model[Param[1]]) < np.max(best_model[Param[1]])):
             slope, intercept, r_value, p_value, std_err = stats.linregress(
                            best_model[Param[1]][mask], mock_output[Param[1]][mask])
             if p_value < 0.001:
                p_value = 0.001
          else:
             slope=0.0
             intercept=1.0
             r_value=0.0
             p_value=0.0

          ax = fig.add_subplot(gs[Param[0]])
          ax.set_xmargin(0.2)
          ax.set_ymargin(0.2)
          ax.errorbar(best_model[Param[1]][mask], mock_output[Param[1]][mask],
            yerr=mock_output[Param[1]+'_err'][mask],
            marker='.',
            label=Param[1]+'\n'+ str("%.1f" % PercentExcluded+'% of outliers ($> 5\sigma$)'),
            color='k', linestyle='None', capsize=0.)
          ax.errorbar(best_model[Param[1]], best_model[Param[1]],
            color='r', linestyle='-', label='1-to-1')
          ax.errorbar(best_model[Param[1]], slope*best_model[Param[1]]+intercept,
            color='b', linestyle='-',
            label='best-fit'+' $r^2$ = '+ str("%.2f" % np.square(r_value))+', outliers not used')
          ax.set_xlabel('Input', fontsize=10)
          ax.set_ylabel('Output', fontsize=10)
          ax.yaxis.labelpad = 0
          if NRows > 2:
             font_size = 4
          else:
             font_size = 6
          ax.legend(fontsize=font_size, loc='best', fancybox=True, framealpha=0.5, numpoints=1)
          ax.tick_params(axis='x', labelsize=8)
          ax.tick_params(axis='y', labelsize=8)
          if nologo is False:
             image = plt.imread(pkg_resources.resource_filename(__name__,
                               "data/CIGALE.png"))
             # Where do we plot CIGALE's logo?
             if NRows == 1:
                x0 = 225
                y0 = 380
             elif NRows == 2:
                x0 = 300
                y0 = 525
             else:
                x0 = 320
                y0 = 525
                            
             fig.figimage(image, x0, y0, origin='upper', zorder=10,
                             alpha=1)

       pdf.savefig()  # saves the current figure into a pdf page
    fig.tight_layout()
    plt.close(fig)

def mock(config, NRowsNCols, nologo):
    """Plot the comparison of input/output values of analysed variables.
    """
    if os.path.isfile(OUT_DIR + BEST_MODEL_FILE):
        input = Table.read(OUT_DIR + BEST_MODEL_FILE, format='ascii')
    else:
        print("*** WARNING ***: No best model file found {}: error."
              .format(OUT_DIR + BEST_MODEL_FILE))

    if os.path.isfile(OUT_DIR + MOCK_OUTPUT_FILE):
        output = Table.read(OUT_DIR + MOCK_OUTPUT_FILE, format='ascii')
    else:
        print("*** WARNING ***: No best model file found {}."
              .format(OUT_DIR + MOCK_OUTPUT_FILE))

    analysed_params = config.configuration['analysis_method_params']
    mock_params = (analysed_params['analysed_variables'])

    NRows = NRowsNCols
    NCols = NRowsNCols
    
    worker(input, output, mock_params, NRows, NCols, nologo)

def main():

    parser = argparse.ArgumentParser(description='Create diagnotic plots for CIGALE')

    parser.add_argument('--NRowsNCols', type=int,  default='2',
                         help='How many columns/rows (1/1 or 2/2) per page in the mosaic?')
    parser.add_argument('--nologo', action="store_true",
                         help='if you do not want the CIGALE logo in the figure')

    args = parser.parse_args()
    config = Configuration()

    mock(config, args.NRowsNCols, args.nologo)
