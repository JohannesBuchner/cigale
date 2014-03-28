# -*- coding: utf-8 -*-
# Copyright (C) 2014 Centre de données Astrophysiques de Marseille
# Copyright (C) 2014 Institute of Astronomy, University of Cambridge
# Copyright (C) 2014 Yannick Roehlly
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Denis Burgarella, Yannick Roehlly, Médéric Boquien

"""
Inter-galactic medium attenuation module
========================================

This module implement the attenuation caused by the inter-galactic medium
depending on the redshift of the galaxy.

For now, only the Meiksin (2006) attenuation is implemented and the module is
parameter-less.  When other method will be implemented, a parameter will allow
to choose the method to use.  This will allow the comparison of the various
attenuations.

"""

import numpy as np
from scipy.misc import factorial
from ..creation_modules import CreationModule


class IgmAttenuation(CreationModule):
    """Add IGM attenuation to the SED

    Add the inter-galactic attenuation to the SED.

    """

    # We cache the IGM transmission computation.
    igm_att = {}

    def process(self, sed):
        """Add IGM attenuation

        Parameters
        ----------
        sed: pcigale.sed.SED object

        """
        redshift = sed.info['redshift']
        if redshift > 0:
            # To memoize the IGM attenuation computation (with depends on the
            #  wavelength grid and the redshift), we have to make the
            # wavelength grid hashable. We do this by converting to a string.
            wave_redshift = (sed.wavelength_grid.tostring(), redshift)

            if wave_redshift not in self.igm_att:
                self.igm_att[wave_redshift] = (
                    igm_transmission_meiksin(sed.wavelength_grid,
                                             redshift) - 1)

            igm_effect = self.igm_att[wave_redshift] * sed .luminosity

            sed.add_module(self.name, self.parameters)
            sed.add_contribution('igm', sed.wavelength_grid, igm_effect)


def igm_transmission_meiksin(wavelength, redshift):
    """Intergalactic transmission (Meiksin, 2006)

    Compute the intergalactic transmission as described in Meiksin, 2006.

    Parameters
    ----------
    wavelength: array like of floats
        The wavelength(s) in nm.
    redshift: float
        The redshift. Must be strictly positive.

    Returns
    -------
    igm_transmission: numpy array of floats
        The intergalactic transmission at each input wavelength.

    """

    redshift = float(redshift)
    wavelength = np.array(wavelength, dtype=float)

    # The redshift must be strictly positive
    if redshift <= 0.:
        raise Exception("The redshift provided must be strictly positive "
                        "<{}>.".format(redshift))

    n_transitions_low = 10
    n_transitions_max = 32
    gamma = 0.2788  # Gamma(0.5,1) i.e., Gamma(2-beta,1) with beta = 1.5
    n0 = 0.25
    lambda_limit = 91.2  # Lyman limit in nm

    lambda_n = np.empty(n_transitions_max)
    z_n = np.empty((n_transitions_max, len(wavelength)))
    for n in range(2, n_transitions_max):
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
    for n in range(3, n_transitions_max):
        if n <= 5:
            w = np.where(z_n[n, :] < 3)
            tau_n[n, w] = (tau_a * fact[n] *
                           np.power(0.25 * (1. + z_n[n, w]), (1. / 3.)))
            w = np.where(z_n[n, :] >= 3)
            tau_n[n, w] = (tau_a * fact[n] *
                           np.power(0.25 * (1. + z_n[n, w]), (1. / 6.)))
        elif 5 < n <= 9:
            tau_n[n, :] = (tau_a * fact[n] *
                           np.power(0.25 * (1. + z_n[n, :]), (1. / 3.)))
        else:
            tau_n[n, :] = (tau_n[9, :] * 720. /
                           (float(n) * (float(n*n - 1.))))

    w = np.where(z_n[2:n_transitions_max, :] > redshift)
    tau_n[w] = 0.
    z_l = wavelength / lambda_limit - 1.
    tau_l_igm = np.zeros_like(wavelength)
    w = np.where(z_l < redshift)
    tau_l_igm[w] = (0.805 * np.power(1. + z_l[w], 3) *
                    (1. / (1. + z_l[w]) - 1. / (1. + redshift)))

    term1 = gamma - np.exp(-1.)

    n = np.arange(n_transitions_low - 1)
    term2 = np.sum(np.power(-1., n) / (factorial(n) * (2*n - 1)))

    w = np.where(z_l < redshift)
    term3 = ((1.+redshift) * np.power(wavelength[w]/lambda_limit, 1.5) -
             np.power(wavelength[w]/lambda_limit, 2.5))

    term4 = np.sum(np.array(
        [((2.*np.power(-1., n) / (factorial(n) * ((6*n - 5)*(2*n - 1)))) *
          ((1.+redshift) ** (2.5-(3 * n)) *
           (wavelength[w]/lambda_limit) ** (3*n) -
           (wavelength[w]/lambda_limit) ** 2.5))
         for n in np.arange(1, n_transitions_low)]), axis=0)

    tau_l_lls = np.zeros_like(wavelength)
    tau_l_lls[w] = n0 * ((term1 - term2) * term3 - term4)

    tau_taun = 0.
    for n in range(n_transitions_max):
        tau_taun = tau_taun + tau_n[n, :]

    lambda_min_igm = (1+redshift)*70.
    weight = np.ones_like(wavelength)
    w = np.where(wavelength < lambda_min_igm)
    weight[w] = np.power(wavelength[w]/lambda_min_igm, 2.)
    # Another weight using erf function can be used.
    # weight[w] = 0.5*(1.+erf(0.05*(wavelength[w]-lambda_min_igm)))

    tau = tau_taun + tau_l_igm + tau_l_lls
    igm_transmission = np.exp(-tau) * weight

    return igm_transmission

# CreationModule to be returned by get_module
Module = IgmAttenuation
