# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013-2014 Institute of Astronomy
# Copyright (C) 2014 Yannick Roehlly <yannick@iaora.eu>
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly & Médéric Boquien

from datetime import datetime
import os
from astropy.table import Table, Column
import numpy as np
from scipy.stats import gaussian_kde
from scipy.linalg import LinAlgError
from ...warehouse import SedWarehouse
import pcigale.analysis_modules.myglobals as gbl

# Directory where the output files are stored
OUT_DIR = "out/"
# Number of points in the PDF
PDF_NB_POINTS = 1000
# Name of the file containing the analysis results
RESULT_FILE = "analysis_results.fits"
# Name of the file containing the best models information
BEST_MODEL_FILE = "best_models.fits"


def gen_pdf(values, probabilities, grid):
    """Generate a probability density function

    For a list of values and associated probabilities, this function
    generates a probability density using a weighted gaussian kernel
    density estimation.

    This part should definitely be improved as it is done in the simplest
    way: each value is repeated (probability * 100) times and the standard
    scipy gaussian KDE is used.

    Parameters
    ----------
    values: array like of floats
        The values of the variable.
    probabilities: array like of floats
        The probability associated with each value
    grid: array like of float
        The list of values to which the probability will be evaluated.

    Returns
    -------
    The list of probabilities evaluated at each value of the grid. If there
    is only one input value with a probability superior to 0, the scipy KDE
    algorithm will fail and None is returned.

    """

    # We use masked arrays because in the analysis module the arrays are
    # already masked and the mask is important.
    values = np.ma.array(values, dtype=float)
    probabilities = np.ma.array(probabilities, dtype=float)

    probabilities = np.ma.array(np.around(100. * probabilities),
                                dtype=int, copy=True)

    combined_values = []

    # We must convert the values masked array to list and test each value
    # against None because 0. is a valid value. For the probabilities,
    # we don't have this problem because if the probability is 0 we won't add
    # the value.
    for val_idx, val in enumerate(values.tolist()):
        if val is not None and probabilities[val_idx]:
            combined_values += [val] * probabilities[val_idx]

    try:
        result = gaussian_kde(combined_values)(grid)
    except (LinAlgError, ValueError):
        result = None

    return result


def save_best_sed(obsid, sed, norm):
    """Save the best SED to a VO table.

    Parameters
    ----------
    obsid: string
        Name of the object. Used to prepend the output file name
    sed: SED object
        Best SED
    norm: float
        Normalisation factor to scale the scale to the observations

    """
    sed.to_votable(OUT_DIR + "{}_best_model.xml".format(obsid), mass=norm)


def save_pdf(obsid, analysed_variables, model_variables, likelihood):
    """Save the PDF to a FITS file

    We estimate the probability density functions (PDF) of the parameters using
    a weighted kernel density estimation. This part should definitely be
    improved as we simulate the weight by adding as many value as their
    probability * 100.

    Parameters
    ----------
    obsid: string
        Name of the object. Used to prepend the output file name
    analysed_variables: list
        Analysed variables names
    model_variables: 2D array
        Analysed variables values for all models
    likelihood: array
        Likelihood for all models

    """
    for var_index, var_name in enumerate(analysed_variables):
        values = model_variables[:, var_index]
        pdf_grid = np.linspace(values.min(), values.max(), PDF_NB_POINTS)
        pdf_prob = gen_pdf(values, likelihood, pdf_grid)

        if pdf_prob is None:
            # TODO: use logging
            print("Can not compute PDF for observation <{}> and "
                  "variable <{}>.".format(obsid, var_name))
        else:
            table = Table((
                Column(pdf_grid, name=var_name),
                Column(pdf_prob, name="probability density")
            ))
            table.write(OUT_DIR + "{}_{}_pdf.fits".format(obsid, var_name))


def save_chi2(obsid, analysed_variables, model_variables, reduced_chi2):
    """Save the best reduced χ² versus the analysed variables

    Parameters
    ----------
    obsid: string
        Name of the object. Used to prepend the output file name
    analysed_variables: list
        Analysed variable names
    model_variables: 2D array
        Analysed variables values for all models
    reduced_chi2:
        Reduced χ²

    """
    for var_index, var_name in enumerate(analysed_variables):
        table = Table((
            Column(model_variables[:, var_index],
                   name=var_name),
            Column(reduced_chi2, name="chi2")))
        table.write(OUT_DIR + "{}_{}_chi2.fits".format(obsid, var_name))


def save_table_analysis(obsid, analysed_variables, analysed_averages,
                        analysed_std):
    """Save the estimated values derived from the analysis of the PDF

    Parameters
    ----------
    obsid: table column
        Names of the objects
    analysed_variables: list
        Analysed variable names
    analysed_averages: array
        Analysed variables values estimates
    analysed_std: array
        Analysed variables errors estimates

    """
    result_table = Table()
    result_table.add_column(Column(obsid.data, name="observation_id"))
    for index, variable in enumerate(analysed_variables):
        result_table.add_column(Column(
            analysed_averages[:, index],
            name=variable
        ))
        result_table.add_column(Column(
            analysed_std[:, index],
            name=variable+"_err"
        ))
    result_table.write(OUT_DIR + RESULT_FILE)


def save_table_best(obsid, chi2, chi2_red, norm, variables, fluxes, filters):
    """Save the values corresponding to the best fit

    Parameters
    ----------
    obsid: table column
        Names of the objects
    chi2: array
        Best χ² for each object
    chi2_red: array
        Best reduced χ² for each object
    norm: array
        Normalisation factor for each object
    variables: list
        All variables corresponding to a SED
    fluxes: 2D array
        Fluxes in all bands for each object
    filters: list
        Filters used to compute the fluxes

    """
    best_model_table = Table()
    best_model_table.add_column(Column(obsid.data, name="observation_id"))
    best_model_table.add_column(Column(chi2, name="chi_square"))
    best_model_table.add_column(Column(chi2_red, name="reduced_chi_square"))

    for index, name in enumerate(gbl.info_keys):
        column = Column([variable[index] for variable in variables], name=name)
        if name in gbl.mass_proportional_info:
            column *= norm
        best_model_table.add_column(column)

    for index, name in enumerate(filters):
        column = Column(fluxes[:, index] * norm, name=name, unit='mJy')
        best_model_table.add_column(column)

    best_model_table.write(OUT_DIR + BEST_MODEL_FILE)


def backup_dir(directory):
    if os.path.exists(directory):
        new_name = datetime.now().strftime("%Y%m%d%H%M") + "_" + directory
        os.rename(directory, new_name)
        print("The existing {} directory was renamed to {}".format(
            OUT_DIR,
            new_name
        ))
    os.mkdir(OUT_DIR)
