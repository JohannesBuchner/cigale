# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2014 Laboratoire d'Astrophysique de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly, Denis Burgarella

"""
Charlot and Fall (2000) power law attenuation module
====================================================

This module implements the attenuation based on a power law as defined
in Charlot and Fall (2000) with a UV bump added.

"""

from collections import OrderedDict

import numpy as np

from . import CreationModule


def power_law(wavelength, delta):
    """Compute the power law (λ / λv)^δ

    Parameters
    ----------
    wavelength: array of float
        Wavelength grid in nm.
    delta: float
        Power law slope.

    Returns
    -------
    a numpy array of floats

    """
    wave = np.array(wavelength)
    return (wave / 550) ** delta


def uv_bump(wavelength, central_wave, gamma, ebump):
    """Compute the Lorentzian-like Drude profile.

    Parameters
    ----------
    wavelength: array of floats
        Wavelength grid in nm.
    central_wave: float
        Central wavelength of the bump in nm.
    gamma: float
        Width (FWHM) of the bump in nm.
    ebump: float
        Amplitude of the bump.

    Returns
    -------
    a numpy array of floats

    """
    return (ebump * wavelength ** 2 * gamma ** 2 /
            ((wavelength ** 2 - central_wave ** 2) ** 2 +
             wavelength ** 2 * gamma ** 2))


def alambda_av(wavelength, delta, bump_wave, bump_width, bump_ampl):
    """Compute the complete attenuation curve A(λ)/Av

    The continuum is a power law (λ / λv) ** δ to which is added a UV bump.
    Over the Lyman continuum, there is no attenuation.

    Parameters
    ----------
    wavelength: array of floats
        The wavelength grid (in nm) to compute the attenuation curve on.
    delta: float
        Slope of the power law.
    bump_wave: float
        Central wavelength (in nm) of the UV bump.
    bump_width: float
        Width (FWHM, in nm) of the UV bump.
    bump_ampl: float
        Amplitude of the UV bump.

    Returns
    -------
    attenuation: array of floats
        The A(λ)/Av attenuation at each wavelength of the grid.

    """
    wave = np.array(wavelength)

    attenuation = power_law(wave, delta)
    attenuation += uv_bump(wavelength, bump_wave, bump_width, bump_ampl)

    # Lyman continuum not attenuated.
    attenuation[wavelength <= 91.2] = 0.

    return attenuation


class PowerLawAtt(CreationModule):
    """Power law attenuation module

    This module computes the attenuation using a power law
    as defined in Charlot and Fall (2000).

    The attenuation can be computed on the whole spectrum or on a specific
    contribution and is added to the SED as a negative contribution.

    """

    parameter_list = OrderedDict([
        ("Av_young", (
            "float",
            "V-band attenuation of the young population.",
            1.
        )),
        ("Av_old_factor", (
            "float",
            "Reduction factor for the V-band attenuation of the old "
            "population compared to the young one (<1).",
            0.44
        )),
        ("uv_bump_wavelength", (
            "float",
            "Central wavelength of the UV bump in nm.",
            217.5
        )),
        ("uv_bump_width", (
            "float",
            "Width (FWHM) of the UV bump in nm.",
            35.
        )),
        ("uv_bump_amplitude", (
            "float",
            "Amplitude of the UV bump. For the Milky Way: 3.",
            0.
        )),
        ("powerlaw_slope", (
            "float",
            "Slope delta of the power law continuum.",
            -0.7
        )),
        ("filters", (
            "string",
            "Filters for which the attenuation will be computed and added to "
            "the SED information dictionary. You can give several filter "
            "names separated by a & (don't use commas).",
            "V_B90 & FUV"
        ))
    ])

    def _init_code(self):
        self.filter_list = [item.strip() for item in
                            self.parameters["filters"].split("&")]
        # We cannot compute the attenuation until we know the wavelengths. Yet,
        # we reserve the object.
        self.sel_attenuation = None

    def process(self, sed):
        """Add the dust attenuation to the SED.

        Parameters
        ----------
        sed: pcigale.sed.SED object

        """
        av = {}
        wavelength = sed.wavelength_grid
        av['young'] = float(self.parameters["Av_young"])
        av['old'] = float(self.parameters["Av_old_factor"] * av['young'])
        uv_bump_wavelength = float(self.parameters["uv_bump_wavelength"])
        uv_bump_width = float(self.parameters["uv_bump_width"])
        uv_bump_amplitude = float(self.parameters["uv_bump_amplitude"])
        powerlaw_slope = float(self.parameters["powerlaw_slope"])

        # Fλ fluxes (only from continuum)) in each filter before attenuation.
        flux_noatt = {filt: sed.compute_fnu(filt) for filt in self.filter_list}

        # Compute attenuation curve
        if self.sel_attenuation is None:
            self.sel_attenuation = alambda_av(wavelength, powerlaw_slope,
                                              uv_bump_wavelength,
                                              uv_bump_width, uv_bump_amplitude)

        attenuation_total = 0.
        contribs = [contrib for contrib in sed.contribution_names if
                    'absorption' not in contrib]
        for contrib in contribs:
            age = contrib.split('.')[-1].split('_')[-1]
            luminosity = sed.get_lumin_contribution(contrib)
            attenuated_luminosity = (luminosity * 10 **
                                     (av[age] * self.sel_attenuation / -2.5))
            attenuation_spectrum = attenuated_luminosity - luminosity
            # We integrate the amount of luminosity attenuated (-1 because the
            # spectrum is negative).
            attenuation = -1 * np.trapz(attenuation_spectrum, wavelength)
            attenuation_total += attenuation

            sed.add_module(self.name, self.parameters)
            sed.add_info("attenuation.Av." + contrib, av[age])
            sed.add_info("attenuation." + contrib, attenuation, True)
            sed.add_contribution("attenuation." + contrib, wavelength,
                                 attenuation_spectrum)

        # Bump and slope of the dust attenuation
        sed.add_info("attenuation.uv_bump_amplitude", uv_bump_amplitude)
        sed.add_info("attenuation.powerlaw_slope", powerlaw_slope)

        # Total attenuation
        if 'dust.luminosity' in sed.info:
            sed.add_info("dust.luminosity",
                         sed.info["dust.luminosity"]+attenuation_total, True,
                         True)
        else:
            sed.add_info("dust.luminosity", attenuation_total, True)

        # Fλ fluxes (only in continuum) in each filter after attenuation.
        flux_att = {filt: sed.compute_fnu(filt) for filt in self.filter_list}

        # Attenuation in each filter
        for filt in self.filter_list:
            sed.add_info("attenuation." + filt,
                         -2.5 * np.log10(flux_att[filt] / flux_noatt[filt]))

# CreationModule to be returned by get_module
Module = PowerLawAtt
