# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013-2014 Institute of Astronomy
# Copyright (C) 2014 Yannick Roehlly <yannick@iaora.eu>
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly & Médéric Boquien

from astropy import log
log.setLevel('ERROR')
from astropy.table import Table, Column
import numpy as np
from scipy.stats import scoreatpercentile

from ..utils import OUT_DIR

# Number of points in the PDF
PDF_NB_POINTS = 1000

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


def save_table_analysis(filename, obsid, analysed_variables, analysed_averages,
                        analysed_std):
    """Save the estimated values derived from the analysis of the PDF

    Parameters
    ----------
    filename: name of the file to save
        Name of the output file
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
    result_table.write(OUT_DIR + filename, format='ascii.fixed_width',
                       delimiter=None)


def save_table_best(filename, obsid, chi2, chi2_red, variables, fluxes, filters,
                    info_keys):
    """Save the values corresponding to the best fit

    Parameters
    ----------
    filename: name of the file to save
        Name of the output file
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

    best_model_table.write(OUT_DIR + filename,format='ascii.fixed_width',
                           delimiter=None)


def dchi2_over_ds2(s):
    """Function used to estimate the normalization factor in the SED fitting
    process when upper limits are included in the dataset to fit (from Eq. A11
    in Sawicki M. 2012, PASA, 124, 1008).

    Parameters
    ----------
    s: Float
        Contains value onto which we perform minimization = normalization
        factor
    obs_fluxes: RawArray
        Contains observed fluxes for each filter.
    obs_errors: RawArray
        Contains observed errors for each filter.
    model_fluxes: RawArray
        Contains modeled fluxes for each filter.
    lim_flag: Boolean
        Tell whether we use upper limits (True) or not (False).

   Returns
    -------
    func: Float
        Eq. A11 in Sawicki M. 2012, PASA, 124, 1008).

    """
    # We enter into this function if lim_flag = True.

    # The mask "data" selects the filter(s) for which measured fluxes are given
    # i.e., when obs_fluxes is >=0. and obs_errors >=0.
    # The mask "lim" selects the filter(s) for which upper limits are given
    # i.e., when obs_fluxes is >=0. and obs_errors = 9990 <= obs_errors < 0.

    wlim = np.where((gbl_obs_errors >= -9990.)&(gbl_obs_errors < 0.))
    wdata = np.where(gbl_obs_errors>=0.)

    mod_fluxes_data = gbl_mod_fluxes[wdata]
    mod_fluxes_lim = gbl_mod_fluxes[wlim]

    obs_fluxes_data = gbl_obs_fluxes[wdata]
    obs_fluxes_lim = gbl_obs_fluxes[wlim]

    obs_errors_data = gbl_obs_errors[wdata]
    obs_errors_lim = -gbl_obs_errors[wlim]

    dchi2_over_ds_data = np.sum(
        (obs_fluxes_data-s*mod_fluxes_data) *
        mod_fluxes_data/(obs_errors_data*obs_errors_data))

    dchi2_over_ds_lim = np.sqrt(2./np.pi)*np.sum(
        mod_fluxes_lim*np.exp(
            -np.square(
                (obs_fluxes_lim-s*mod_fluxes_lim)/(np.sqrt(2)*obs_errors_lim)
                      )
                             )/(
            obs_errors_lim*(
                1.+erf(
                   (obs_fluxes_lim-s*mod_fluxes_lim)/(np.sqrt(2)*obs_errors_lim)
                      )
                           )
                               )
                                                )

    func = dchi2_over_ds_data - dchi2_over_ds_lim

    return func

def analyse_chi2(chi2):
    """Function to analyse the best chi^2 and find out whether what fraction of
    objects seem to be overconstrainted.

    Parameters
    ----------
    chi2: RawArray
        Contains the reduced chi^2

    """
    chi2_red = np.ctypeslib.as_array(chi2[0])
    # If low values of reduced chi^2, it means that the data are overfitted
    # Errors might be under-estimated or not enough valid data.
    print("\n{}% of the objects have chi^2_red~0 and {}% chi^2_red<0.5"
        .format(np.round((chi2_red < 1e-12).sum()/chi2_red.size, 1),
                np.round((chi2_red < 0.5).sum()/chi2_red.size, 1)))
