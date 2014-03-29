# -*- coding: utf-8 -*-
# Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Authors: Yannick Roehlly, Médéric Boquien

import numpy as np
from scipy.constants import c, pi


def lambda_to_nu(wavelength):
    """Convert wavelength (nm) to frequency (Hz)

    Parameters
    ----------
    wavelength: float or array of floats
        The wavelength(s) in nm.

    Returns
    -------
    nu: float or array of floats
        The frequency(ies) in Hz.

    """
    return c / (wavelength * 1.e-9)


def nu_to_lambda(frequency):
    """Convert frequency (Hz) to wavelength (nm)

    Parameters
    ----------
    frequency: float or numpy.array of floats
        The frequency(ies) in Hz.

    Returns
    -------
    wavelength: float or numpy.array of floats
        The wavelength(s) in nm.

    """
    return 1.e-9 * c / frequency


def best_grid(wavelengths1, wavelengths2):
    """
    Return the best wavelength grid to regrid to arrays

    Considering the two wavelength grids passed in parameters, this function
    compute the best new grid that will be used to regrid the two spectra
    before combining them. We do not use np.unique as it is much slowe than
    finding the unique elements by hand.

    Parameters
    ----------
    wavelengths1, wavelengths2: array of floats
        The wavelength grids to be 'regridded'.

    Returns
    -------
    new_grid: array of floats
        Array containing all the wavelengths found in the input arrays.

    """
    wl = np.concatenate((wavelengths1, wavelengths2))
    wl.sort()
    flag = np.ones(len(wl), dtype=bool)
    np.not_equal(wl[1:], wl[:-1], out=flag[1:])

    return wl[flag]


def luminosity_to_flux(luminosity, dist):
    """
    Convert a luminosity (or luminosity density) to a flux (or flux density).

    F = L / (4πDl2)

    Parameters
    ----------
    luminosity: float or array of floats
        Luminosity (typically in W) or luminosity density (W/nm or W/Hz).
    dist: float
        Luminosity distance of the object in metres

    Returns
    -------
    flux: float or array of floats
        The flux (typically in W/m²) of flux density (W/m²/nm or W/m²/Hz).

    """

    flux = luminosity / (4. * pi * dist * dist)

    return flux


def lambda_flambda_to_fnu(wavelength, flambda):
    """
    Convert a Fλ vs λ spectrum to Fν vs λ

    Parameters
    ----------
    wavelength: list-like of floats
        The wavelengths in nm.
    flambda: list-like of floats
        Fλ flux density in W/m²/nm (or Lλ luminosity density in W/nm).

    Returns
    -------
    fnu: array of floats
        The Fν flux density in mJy (or the Lν luminosity density in
        1.e-29 W/Hz).

    """
    wavelength = np.array(wavelength, dtype=float)
    flambda = np.array(flambda, dtype=float)

    # Factor 1e+29 is to switch from W/m²/Hz to mJy
    # Factor 1e-9 is to switch from nm to m (only one because the other nm
    # wavelength goes with the Fλ in W/m²/nm).
    fnu = 1e+29 * 1e-9 * flambda * wavelength * wavelength / c

    return fnu


def lambda_fnu_to_flambda(wavelength, fnu):
    """
    Convert a Fν vs λ spectrum to Fλ vs λ

    Parameters
    ----------
    wavelength: list-like of floats
        The wavelengths in nm.
    fnu: list-like of floats
        The Fν flux density in mJy (of the  Lν luminosity density in
        1.e-29 W/Hz).

    Returns
    -------
    flambda: array of floats
        Fλ flux density in W/m²/nm (or Lλ luminosity density in W/nm).

    """
    wavelength = np.array(wavelength, dtype=float)
    fnu = np.array(fnu, dtype=float)

    # Factor 1e-29 is to switch from Jy to W/m²/Hz
    # Factor 1e+9 is to switch from m to nm
    flambda = 1e-29 * 1e+9 * fnu / (wavelength * wavelength) * c

    return flambda


def redshift_spectrum(wavelength, flux, redshift, is_fnu=False):
    """Redshit a spectrum

    Parameters
    ----------
    wavelength: array like of floats
        The wavelength in nm.
    flux: array like of floats
        The flux or luminosity density.
    redshift: float
        The redshift.
    is_fnu: boolean
        If false (default) the flux is a Fλ density in W/m²/nm (or a Lλ
        luminosity density in W/nm). If true, the flux is a Fν density in mJy
        (or a Lν luminosity density in 1.e-29 W/Hz).

    Results
    -------
    wavelength, flux: tuple of numpy arrays of floats
        The redshifted spectrum with the same kind of flux (or luminosity)
        density as the input.

    """
    wavelength = np.array(wavelength, dtype=float)
    flux = np.array(flux, dtype=float)
    redshift = float(redshift)

    if redshift < 0:
        redshift_factor = 1. / (1. - redshift)
    else:
        redshift_factor = 1. + redshift

    if is_fnu:
        # Switch to Fλ
        flux = lambda_fnu_to_flambda(wavelength, flux)

    wavelength *= redshift_factor
    flux /= redshift_factor

    if is_fnu:
        # Switch back to Fλ
        flux = lambda_flambda_to_fnu(wavelength, flux)

    return wavelength, flux
