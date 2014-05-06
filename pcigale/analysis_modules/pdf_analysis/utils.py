# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013-2014 Institute of Astronomy
# Copyright (C) 2014 Yannick Roehlly <yannick@iaora.eu>
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly & Médéric Boquien

from astropy.table import Table, Column
import numpy as np
from scipy.stats import scoreatpercentile

from ..utils import OUT_DIR

# Number of points in the PDF
PDF_NB_POINTS = 1000
# Name of the file containing the analysis results
RESULT_FILE = "analysis_results.txt"
# Name of the file containing the best models information
BEST_MODEL_FILE = "best_models.txt"


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
    likelihood: 2D array
        Likelihood for all models

    """
    for var_index, var_name in enumerate(analysed_variables):
        pdf_grid = model_variables[:, var_index]
        pdf_prob = likelihood[:, var_index]
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
    """Save the best reduced Ç² versus the analysed variables

    Parameters
    ----------
    obsid: string
        Name of the object. Used to prepend the output file name
    analysed_variables: list
        Analysed variable names
    model_variables: 2D array
        Analysed variables values for all models
    reduced_chi2:
        Reduced Ç²

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
    analysed_averages: RawArray
        Analysed variables values estimates
    analysed_std: RawArray
        Analysed variables errors estimates

    """
    np_analysed_averages = np.ctypeslib.as_array(analysed_averages[0])
    np_analysed_averages = np_analysed_averages.reshape(analysed_averages[1])

    np_analysed_std = np.ctypeslib.as_array(analysed_std[0])
    np_analysed_std = np_analysed_std.reshape(analysed_std[1])

    result_table = Table()
    result_table.add_column(Column(obsid.data, name="observation_id"))
    for index, variable in enumerate(analysed_variables):
        result_table.add_column(Column(
            np_analysed_averages[:, index],
            name=variable
        ))
        result_table.add_column(Column(
            np_analysed_std[:, index],
            name=variable+"_err"
        ))
    result_table.write(OUT_DIR + RESULT_FILE, format='ascii.commented_header')


def save_table_best(obsid, chi2, chi2_red, variables, fluxes, filters,
                    info_keys):
    """Save the values corresponding to the best fit

    Parameters
    ----------
    obsid: table column
        Names of the objects
    chi2: RawArray
        Best χ² for each object
    chi2_red: RawArray
        Best reduced χ² for each object
    variables: RawArray
        All variables corresponding to a SED
    fluxes: RawArray
        Fluxes in all bands for each object
    filters: OrderedDict
        Filters used to compute the fluxes
    info_keys: list
        Parameters names

    """
    np_fluxes = np.ctypeslib.as_array(fluxes[0])
    np_fluxes = np_fluxes.reshape(fluxes[1])

    np_variables = np.ctypeslib.as_array(variables[0])
    np_variables = np_variables.reshape(variables[1])

    np_chi2 = np.ctypeslib.as_array(chi2[0])

    np_chi2_red = np.ctypeslib.as_array(chi2_red[0])

    best_model_table = Table()
    best_model_table.add_column(Column(obsid.data, name="observation_id"))
    best_model_table.add_column(Column(np_chi2, name="chi_square"))
    best_model_table.add_column(Column(np_chi2_red, name="reduced_chi_square"))

    for index, name in enumerate(info_keys):
        column = Column(np_variables[:, index], name=name)
        best_model_table.add_column(column)

    for index, name in enumerate(filters):
        column = Column(np_fluxes[:, index], name=name, unit='mJy')
        best_model_table.add_column(column)

    best_model_table.write(OUT_DIR + BEST_MODEL_FILE,
                           format='ascii.commented_header')


def FDbinSize(values):
    """
    To define the size of the bin (parameter x), we use the Freedman-Diaconis
    rule : bin size = 2 * IQR(x) * N^(-1/3) where IQR is the InterQuartile
    Range containing 50% of sample. We do not use it here but, for a normal
    distribution IQR = 1.349 * sigma. Note that the actual rule, there is a
    factor 2. and not 1 like here.

    Parameters
    ----------
    values: array like of floats
        The values of the variable.

    Returns
    -------
    h:  float
      The Freedman-Diaconis bin size

    """
    # First Calculate the interquartile range
    values = np.sort(values)
    upperQuartile = scoreatpercentile(values, 75.)
    lowerQuartile = scoreatpercentile(values, 25.)
    IQR = upperQuartile - lowerQuartile

    # Find the Freedman-Diaconis bin size
    h = 2. * IQR * len(values)**(-1./3.)

    return h
