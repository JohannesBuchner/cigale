# -*- coding: utf-8 -*-
"""
Copyright (C) 2012 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>
@author: Médéric Boquien <mederic.boquien@oamp.fr>

"""


import numpy as np
from scipy import integrate
from scipy.constants import c, pi, parsec


def lambda_to_nu(wavelength):
    """Convert wavelength (nm) to frequency (Hz)

    Parameters
    ----------
    wavelength : float or array of floats
        The wavelength(s) in nm.

    Returns
    -------
    nu : float or array of floats
        The frequency(ies) in Hz.

    """
    return c / (wavelength * 1.e-9)


def nu_to_lambda(frequency):
    """Convert frequency (Hz) to wavelength (nm)

    Parameters
    ----------
    frequency : float or numpy.array of floats
        The frequency(ies) in Hz.

    Returns
    -------
    wavelength : float or numpy.array of floats
        The wavelength(s) in nm.

    """
    return 1.e-9 * c / frequency


def best_grid(wavelengths1, wavelengths2):
    """
    Return the best wavelength grid to regrid to arrays

    Considering the two wavelength grids passed in parameters, this function
    compute the best new grid that will be used to regrid the two spectra
    before combining them.

    Parameters
    ----------
    wavelengths1, wavelengths2 : array of floats
        The wavelength grids to be 'regrided'.

    Returns
    -------
    new_grid : array of floats
        Array containing all the wavelengths found in the input arrays.

    """
    new_grid = np.hstack((wavelengths1, wavelengths2))
    new_grid.sort()
    new_grid = np.unique(new_grid)

    return new_grid


def luminosity_distance(z, h0=71., omega_m=0.27, omega_l=0.73):
    """
    Computes luminosity distance at redshift z in Mpc for given Λ cosmology
    (H_0 in (km/s)/Mpc, Ω_M, and Ω_Λ) Ref.: Hogg (1999) astro-ph/9905116

    Parameters
    ----------
    z : float
        Redshift
    h0 : float
        Hubble's constant
    omega_m : float
        Omega matter.
    omega_l : float
        Omega vacuum

    Returns
    -------
    luminosity_distance : float
        The luminosity distance in Mpc.

    """

    omega_k = 1. - omega_m - omega_l

    if z > 0.:
        dist, edist = integrate.quad(
            lambda x: (omega_m * (1. + x) ** 3
                       + omega_k * (1 + x) ** 2 + omega_l) ** (-.5),
            0.,
            z,
            epsrel=1e-3)
    else:
        # Bad idea as there is something *wrong* going on
        print('LumDist: z <= 0 -> Assume z = 0!')
        z = 0.
        dist = 0.

    if omega_k > 0.:
        dist = np.sinh(dist * np.sqrt(omega_k)) / np.sqrt(omega_k)
    elif omega_k < 0.:
        dist = np.sin(dist * np.sqrt(-omega_k)) / np.sqrt(-omega_k)

    return c / (h0 * 1.e3) * (1. + z) * dist


def luminosity_to_flux(luminosity, redshift=0):
    """
    Convert a luminosity (or luminosity density) to a flux (or flux density).

    F = L / (4πDl2)

    Parameters
    ----------
    luminosity : float or array of floats
        Luminosity (typically in W) or luminosity density (W/nm or W/Hz).
    redshift :
        Redshift. If redshift is 0 (the default) the flux at a luminosity
        distance of 10 pc is returned.

    Returns
    -------
    flux : float or array of floats
        The flux (typically in W/m²) of flux density (W/m²/nm or W/m²/Hz).

    """
    if redshift == 0:
        dist = 10 * parsec
    else:
        dist = luminosity_distance(redshift) * 1.e6 * parsec

    return luminosity / (4 * pi * np.square(dist))


def redshift_wavelength(wavelength, redshift):
    """Redshift a wavelength grid

    Parameters
    ----------
    wavelength : array of floats
        Wavelength vector.
    redshift : float
        Redshift.

    Returns
    -------
    redshifted_wavelength : array of floats
        Redshifted wavelength grid.

    """
    if redshift < 0:
        return wavelength / (1.0 - redshift)
    else:
        return wavelength * (1.0 + redshift)


def lambda_flambda_to_lambda_fnu(spectrum):
    """
    Convert a Fλ vs λ spectrum to Fν vs λ

    Parameters
    ----------
    spectrum : array of floats
        spectrum[0] must contain the wavelength in nm and spectrum[1] must
        contain the Fλ flux in erg/cm^2/s/nm.

    Returns
    -------
    lambda_fnu : array of floats
        lambda_fnu[0] contains the wavelength in nm and lambda_fnu[1] contains
        the Fν flux in Jansky

    """
    wavelength, flambda = spectrum
    # Factor 1e+23 is to switch from erg/s/cm^2/Hz to Jy
    # Factor 1e-9 is to switch from nm to m (only one because the other nm
    # wavelength goes with the Fλ in ergs/s/cm^2/nm).
    fnu = 1e+23 * 1e-9 * flambda * wavelength * wavelength / c

    return np.vstack((wavelength, fnu))


def lambda_fnu_to_lambda_flambda(spectrum):
    """
    Convert a Fν vs λ spectrum to Fλ vs λ

    Parameters
    ----------
    spectrum : array of floats
        spectrum[0] must contain the wavelength in nm and spectrum[1] must
        contain the Fν flux in Jansky

    Returns
    -------
    lambda_flambda : array of floats
        lambda_flambda[0] contains the wavelength in nm and lambda_flambda[1]
        contains the Fλ flux in erg/cm^2/s/nm.

    """
    wavelength, fnu = spectrum
    # Factor 1e-23 is to switch from Jy to erg/s/cm^2/Hz
    # Factor 1e+9 is to switch from m to nm
    flambda = 1e-23 * 1e+9 * fnu / (wavelength * wavelength) * c

    return np.vstack((wavelength, flambda))


def redshift_spectrum(spectrum, redshift, dimming=False, is_fnu=False):
    """
    Redshit a spectrum, optionally adding cosmological dimming

    FIXME: Is this usefull?

    Parameters
    ----------
    spectrum : array of floats
        spectrum[0] must contain the wavelength in nm and spectrum[1] must
        contain the flux. The default is to have Fλ in erg/cm^2/s/nm if is_fnu
        is set to true, the Fν in Jansky is expected (it's only important when
        dimming).

    dimming : boolean
        If set to true, the cosmological dimming is applied to the fluxes.

    is_fnu : boolean
        If set to true, the flux are Fν fluxes, else they are assumed to be Fλ.

    Results
    -------
    spectrum : array of floats
        The redshifted spectrum with the same kind of fluxes as the input.

    """

    wavelength = redshift_wavelength(spectrum[0])
    flux = np.copy(spectrum[1])

    if dimming:
        # If the flux is Fnu, we must switch to Flambda to compute the
        # dimming.
        if is_fnu:
            flux = lambda_fnu_to_lambda_flambda(spectrum)[1]

        # Now flux is Flambda, we can apply cosmological dim.
        if redshift < 0:
            flux = flux * (1.0 - redshift)
        else:
            flux = flux / (1.0 + redshift)

        # If the initial flux was Fnu, convert it back from Flambda
        if is_fnu:
            flux = lambda_flambda_to_lambda_fnu(
                np.vstack((wavelength, flux)))[:, 1]

    return np.vstack((wavelength, flux))
