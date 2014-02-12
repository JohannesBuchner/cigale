# -*- coding: utf-8 -*-
# Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Authors: Denis Burgarella, Yannick Roehlly

import math
import numpy as np


def ig_transmission_meiksin(wavelength, redshift):
    """Intergalactic transmission (Meiksin, 2006)

    Compute the intergalactic transmission as described in Meiksin, 2006.

    Parameters
    ----------
    wavelength : array like of floats
        The wavelength(s) in nm.
    redshift : float
        The redshift.

    Returns
    -------
    ig_transmission : numpy array of floats
        The intergalactic transmission at each input wavelength.

    """

    redshift = float(redshift)
    wavelength = np.array(wavelength, dtype=float)

    # For redshift below 2, we assume that there is no IGM attenuation.
    if redshift <= 2.:
        ig_tranmission = np.ones(len(wavelength))
        return ig_tranmission

    n_transitions_low = 10
    n_transitions_max = 31
    gamma = 0.2788  # Gamma(0.5,1) i.e., Gamma(2-beta,1) with beta = 1.5
    n0 = 0.25
    lambda_limit = 91.2  # Lyman limit in nm

    lambda_n = np.zeros(n_transitions_low)
    z_n = np.zeros((n_transitions_low, len(wavelength)))
    for n in range(2, n_transitions_low):
        lambda_n[n] = lambda_limit / (1. - 1. / float(n*n))
        z_n[n, :] = wavelength / lambda_n[n] - 1.

    # From Table 1 in Meiksin (2006), only n >=3 are relevant. fact has a
    # length equal to n_transistions_low.
    fact = np.array([1., 1., 1., 0.348, 0.179, 0.109, 0.0722, 0.0508, 0.0373,
                     0.0283])

    # First, tau_alpha is the mean Lyman alpha transmitted flux,
    # Here n = 2 => tau_2 = tau_alpha
    tau_n = np.zeros((n_transitions_max, len(wavelength)))
    if redshift <= 4:
        tau_a = 0.00211 * np.power(1. + redshift, 3.7)
        tau_n[2, :] = 0.00211 * np.power(1. + z_n[2, :], 3.7)
    elif redshift >= 4:
        tau_a = 0.00058 * np.power(1. + redshift, 4.5)
        tau_n[2, :] = 0.00058 * np.power(1. + z_n[2, :], 4.5)

    # Then, tau_n is the mean optical depth value for transitions
    # n = 3 - 9 -> 1
    for n in range(3, n_transitions_low):
        if n <= 5:
            for l, index in enumerate(wavelength):
                if z_n[n, l] < 3:
                    tau_n[n, l] = (tau_a * fact[n] *
                                   np.power(0.25 * (1. + z_n[n, l]),
                                            (1. / 3.)))
                else:
                    tau_n[n, l] = (tau_a * fact[n] *
                                   np.power(0.25 * (1. + z_n[n, l]),
                                            (1. / 6.)))

        elif 5 < n <= 9:
            for l, index in enumerate(wavelength):
                tau_n[n, l] = (tau_a * fact[n] *
                               np.power(0.25 * (1. + z_n[n, l]),
                                        (1. / 3.)))

        else:
            for l, index in enumerate(wavelength):
                tau_n[n, l] = (tau_n[9, l] * 720. /
                               (float(n) * (float(n*n - 1.))))

    for n in range(2, n_transitions_low):
        for l, index in enumerate(wavelength):
            if z_n[n, l] >= redshift:
                tau_n[n, l] = 0.

    z_l = wavelength / lambda_limit - 1.
    tau_l_igm = np.zeros(len(wavelength))
    for l, index in enumerate(wavelength):
        if z_l[l] >= redshift:
            tau_l_igm[l] = 0.
        else:
            tau_l_igm[l] = (0.805 * np.power(1. + z_l[l], 3) *
                            (1. / (1. + z_l[l]) - 1. / (1. + redshift)))

    term1 = gamma - np.exp(-1.)
    term2 = 0.
    for n in range(n_transitions_max):
        term2 = term2 + np.power(-1., n) / (math.factorial(n) * float(2*n - 1))
    term3 = ((1.+redshift) * np.power(wavelength/lambda_limit, 1.5) -
             np.power(wavelength/lambda_limit, 2.5))
    term4 = np.zeros(len(wavelength))
    for n in range(1, n_transitions_max):
        term4 = (term4 + (2.*np.power(-1., n) / (math.factorial(n) *
                                                 float((6*n - 5)*(2*n - 1)))) *
                 (np.power(1.+redshift, 2.5-float(3 * n)) *
                  np.power(wavelength/lambda_limit, 3*n) -
                  np.power(wavelength/lambda_limit, 2.5)))

    tau_l_lls = np.zeros(len(wavelength))
    for l, index in enumerate(wavelength):
        if z_l[l] >= redshift:
            tau_l_lls[l] = 0.
        else:
            tau_l_lls[l] = n0 * ((term1 - term2) * term3[l] - term4[l])

    tau_taun = 0.
    for n in range(n_transitions_max):
        tau_taun = tau_taun + tau_n[n, :]

    weight = np.zeros(len(wavelength))
    for l, index in enumerate(wavelength):
        weight[l] = 0.5*(1.+math.erf(0.01*(wavelength[l]-(1+redshift)*60.)))

    tau = tau_taun + tau_l_igm + tau_l_lls
    ig_tranmission = np.exp(- tau) * weight

    return ig_tranmission
