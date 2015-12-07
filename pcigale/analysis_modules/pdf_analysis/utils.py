# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013-2014 Institute of Astronomy
# Copyright (C) 2014 Yannick Roehlly <yannick@iaora.eu>
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly & Médéric Boquien

from astropy import log
from astropy.table import Table, Column
import numpy as np
from scipy import optimize
from scipy.special import erf

from ..utils import OUT_DIR

log.setLevel('ERROR')


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


def save_pdf(obsid, name, model_variable, likelihood):
    """Compute and save the PDF to a FITS file

    We estimate the probability density functions (PDF) of the parameter from
    a likelihood-weighted histogram.

    Parameters
    ----------
    obsid: string
        Name of the object. Used to prepend the output file name
    name: string
        Analysed variable name
    model_variable: array
        Values of the model variable
    likelihood: 1D array
        Likelihood of the "likely" models

    """

    # We check how many unique parameter values are analysed and if
    # less than Npdf (= 100), the PDF is initally built assuming a
    # number of bins equal to the number of unique values for a given
    # parameter
    Npdf = 100
    min_hist = np.min(model_variable)
    max_hist = np.max(model_variable)
    Nhist = min(Npdf, len(np.unique(model_variable)))

    if min_hist == max_hist:
        pdf_grid = np.array([min_hist, max_hist])
        pdf_prob = np.array([1., 1.])
    else:
        pdf_prob, pdf_grid = np.histogram(model_variable, Nhist,
                                          (min_hist, max_hist),
                                          weights=likelihood, density=True)
        pdf_x = (pdf_grid[1:]+pdf_grid[:-1]) / 2.

        pdf_grid = np.linspace(min_hist, max_hist, Npdf)
        pdf_prob = np.interp(pdf_grid, pdf_x, pdf_prob)

    if pdf_prob is None:
        print("Can not compute PDF for observation <{}> and variable <{}>."
              "".format(obsid, name))
    else:
        table = Table((
            Column(pdf_grid, name=name),
            Column(pdf_prob, name="probability density")
        ))
        table.write(OUT_DIR + "{}_{}_pdf.fits".format(obsid, name))


def _save_chi2(obsid, name, model_variable, chi2):
    """Save the best reduced χ² versus an analysed variable

    Parameters
    ----------
    obsid: string
        Name of the object. Used to prepend the output file name
    name: string
        Analysed variable name
    model_variable: array
        Values of the model variable
    chi2:
        Reduced χ²

    """
    table = Table((Column(model_variable, name=name),
                   Column(chi2, name="chi2")))
    table.write(OUT_DIR + "{}_{}_chi2.fits".format(obsid, name))


def save_chi2(obsid, names, mass_proportional, model_variables, scaling, chi2):
    """Save the best reduced χ² versus analysed variables

    Parameters
    ----------
    obsid: string
        Name of the object. Used to prepend the output file name
    name: list of strings
        Analysed variables names
    model_variables: array
        Values of the model variables
    scaling: array
        Scaling factors of the models
    chi2:
        Reduced χ²
    """
    for i, name in enumerate(names):
        if name.endswith('_log'):
            if name[:-4] in mass_proportional:
                model_variable = np.log10(model_variables[:, i] * scaling)
            else:
                model_variable = np.log10(model_variables[:, i])
        else:
            if name in mass_proportional:
                model_variable = model_variables[:, i] * scaling
            else:
                model_variable = model_variables

        _save_chi2(obsid, name, model_variable, chi2)


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


def save_table_best(filename, obsid, chi2, chi2_red, variables, fluxes,
                    filters, info_keys):
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
    filters: list
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

    best_model_table.write(OUT_DIR + filename, format='ascii.fixed_width',
                           delimiter=None)


def dchi2_over_ds2(s, obs_fluxes, obs_errors, mod_fluxes):
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
    # i.e., when obs_errors < 0

    wlim = np.where(np.isfinite(obs_errors) & (obs_errors < 0.))
    wdata = np.where(obs_errors >= 0.)

    mod_fluxes_data = mod_fluxes[wdata]
    mod_fluxes_lim = mod_fluxes[wlim]

    obs_fluxes_data = obs_fluxes[wdata]
    obs_fluxes_lim = obs_fluxes[wlim]

    obs_errors_data = obs_errors[wdata]
    obs_errors_lim = -obs_errors[wlim]

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


def _compute_scaling(model_fluxes, obs_fluxes, obs_errors):
    """Compute the scaling factor to be applied to the model fluxes to best fit
    the observations. Note that we look over the bands to avoid the creation of
    an array of the same size as the model_fluxes array. Because we loop on the
    bands and not on the models, the impact on the performance should be small.

    Parameters
    ----------
    model_fluxes: array
        Fluxes of the models
    obs_fluxes: array
        Observed fluxes
    obs_errors: array
        Observed errors

    Returns
    -------
    scaling: array
        Scaling factors minimising the χ²
    """
    num = np.zeros(model_fluxes.shape[0])
    denom = np.zeros(model_fluxes.shape[0])
    for i in range(obs_fluxes.size):
        if np.isfinite(obs_fluxes[i]):
            num += model_fluxes[:, i] * (obs_fluxes[i] / (obs_errors[i] *
                                                          obs_errors[i]))
            denom += np.square(model_fluxes[:, i] / obs_errors[i])

    return num/denom


def compute_chi2(model_fluxes, obs_fluxes, obs_errors, lim_flag):
    """Compute the χ² of observed fluxes with respect to the grid of models. We
    take into account upper limits if need be. Note that we look over the bands
    to avoid the creation of an array of the same size as the model_fluxes
    array. Because we loop on the bands and not on the models, the impact on
    the performance should be small.

    Parameters
    ----------
    model_fluxes: array
        2D grid containing the fluxes of the models
    obs_fluxes: array
        Fluxes of the observed object
    obs_errors: array
        Uncertainties on the fluxes of the observed object
    lim_flag: boolean
        Boolean indicating whether upper limits should be treated (True) or
        discarded (False)

    Returns
    -------
    chi2: array
        χ² for all the models in the grid
    scaling: array
        scaling of the models to obtain the minimum χ²
    """
    scaling = _compute_scaling(model_fluxes, obs_fluxes, obs_errors)

    # χ² of the comparison of each model to each observation.
    chi2 = np.zeros(model_fluxes.shape[0])
    for i in range(obs_fluxes.size):
        if np.isfinite(obs_fluxes[i]) and obs_errors[i] > 0.:
            chi2 += np.square(
                (obs_fluxes[i] - model_fluxes[:, i] * scaling) / obs_errors[i])

    # Some observations may not have flux values in some filter(s), but
    # they can have upper limit(s).
    if (lim_flag and np.any(obs_errors <= 0.)) == True:
        for imod in range(len(model_fluxes)):
            scaling[imod] = optimize.root(dchi2_over_ds2, scaling[imod],
                                          args=(obs_fluxes, obs_errors,
                                                model_fluxes[imod, :])).x
        mask_lim = (obs_errors <= 0.)
        chi2 += -2. * np.sum(
            np.log(
                np.sqrt(np.pi/2.)*(-obs_errors[mask_lim])*(
                    1.+erf(
                        (obs_fluxes[mask_lim]-model_fluxes[:, mask_lim] *
                         scaling[:, np.newaxis]) /
                        (np.sqrt(2)*(-obs_errors[mask_lim]))))), axis=1)

    return chi2, scaling


def weighted_param(param, weights):
    """Compute the weighted mean and standard deviation of an array of data.

    Parameters
    ----------
    param: array
        Values of the parameters for the entire grid of models
    weights: array
        Weights by which to weight the parameter values

    Returns
    -------
    mean: float
        Weighted mean of the parameter values
    std: float
        Weighted standard deviation of the parameter values

    """

    mean = np.average(param, weights=weights)
    std = np.sqrt(np.average((param-mean)**2, weights=weights))

    return (mean, std)
