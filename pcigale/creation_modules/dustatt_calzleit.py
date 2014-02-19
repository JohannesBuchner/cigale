# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

"""
Calzetti et al. (2000) and Leitherer et al. (2002) attenuation module
=====================================================================

This module implements the Calzetti et al. (2000) and  Leitherer et al. (2002)
attenuation formulae, adding an UV-bump and a power law.

"""

import numpy as np
from collections import OrderedDict
from . import CreationModule
from ..data import Database


def k_calzetti2000(wavelength):
    """Compute the Calzetti et al. (2000) A(λ)/E(B-V)∗

    Given a wavelength grid, this function computes the selective attenuation
    A(λ)/E(B-V)∗ using the formula from Calzetti at al. (2000). This formula
    is given for wavelengths between 120 nm and 2200 nm, but this function
    makes the computation outside.

    Parameters
    ----------
    wavelength : array of floats
        Wavelength grid in nm.

    Returns
    -------
    a numpy array of floats

    """
    wavelength = np.array(wavelength)
    result = np.zeros(len(wavelength))

    # Attenuation between 120 nm and 630 nm
    mask = (wavelength < 630)
    result[mask] = 2.659 * (-2.156 + 1.509e3 / wavelength[mask]
                            - 0.198e6 / wavelength[mask] ** 2
                            + 0.011e9 / wavelength[mask] ** 3) + 4.05

    # Attenuation between 630 nm and 2200 nm
    mask = (wavelength >= 630)
    result[mask] = 2.659 * (-1.857 + 1.040e3 / wavelength[mask]) + 4.05

    return result


def k_leitherer2002(wavelength):
    """Compute the Leitherer et al. (2002) A(λ)/E(B-V)∗

    Given a wavelength grid, this function computes the selective attenuation
    A(λ)/E(B-V)∗ using the formula from Leitherer at al. (2002). This formula
    is given for wavelengths between 91.2 nm and 180 nm, but this function
    makes the computation outside.

    Parameters
    ----------
    wavelength : array of floats
        Wavelength grid in nm.

    Returns
    -------
    a numpy array of floats

    """
    wavelength = np.array(wavelength)
    result = (5.472 + 0.671e3 / wavelength
              - 9.218e3 / wavelength ** 2
              + 2.620e6 / wavelength ** 3)

    return result


def uv_bump(wavelength, central_wave, gamma, ebump):
    """Compute the Lorentzian-like Drude profile.

    Parameters
    ----------
    wavelength : array of floats
        Wavelength grid in nm.
    central_wave : float
        Central wavelength of the bump in nm.
    gamma : float
        Width (FWHM) of the bump in nm.
    ebump : float
        Amplitude of the bump.

    Returns
    -------
    a numpy array of floats

    """
    return (ebump * wavelength ** 2 * gamma ** 2 /
            ((wavelength ** 2 - central_wave ** 2) ** 2
             + wavelength ** 2 * gamma ** 2))


def power_law(wavelength, delta):
    """Power law 'centered' on 550 nm..

    Parameters
    ----------
    wavelength : array of floats
        The wavelength grid in nm.
    delta : float
        The slope of the power law.

    Returns
    -------
    array of floats

    """
    return (wavelength / 550) ** delta


def a_vs_ebv(wavelength, bump_wave, bump_width, bump_ampl, power_slope):
    """Compute the complete attenuation curve A(λ)/E(B-V)*

    The Leitherer et al. (2002) formula is used bellow 150 nm (even if it is
    defined only after 91.2 nm) and the Calzetti et al. (2000) formula is used
    after 150 (we do an extrapolation after 2200 nm). When the attenuation
    becomes negative, it is kept to 0. This continuum is multiplied by the
    power law and then the UV bump is added.

    Parameters
    ----------
    wavelength : array of floats
        The wavelength grid (in nm) to compute the attenuation curve on.
    bump_wave : float
        Central wavelength (in nm) of the UV bump.
    bump_width : float
        Width (FWHM, in nm) of the UV bump.
    bump_ampl : float
        Amplitude of the UV bump.
    power_slope : float
        Slope of the power law.

    Returns
    -------
    attenuation : array of floats
        The A(λ)/E(B-V)* attenuation at each wavelength of the grid.

    """
    attenuation = np.zeros(len(wavelength))

    # Leitherer et al.
    mask = (wavelength < 150)
    attenuation[mask] = k_leitherer2002(wavelength[mask])
    # Calzetti et al.
    mask = (wavelength >= 150)
    attenuation[mask] = k_calzetti2000(wavelength[mask])
    # We set attenuation to 0 where it becomes negative
    mask = (attenuation < 0)
    attenuation[mask] = 0
    # Power law
    attenuation = attenuation * power_law(wavelength, power_slope)
    # UV bump
    attenuation = attenuation + uv_bump(wavelength, bump_wave,
                                        bump_width, bump_ampl)

    return attenuation


class CalzLeit(CreationModule):
    """Calzetti + Leitherer attenuation module

    This module computes the Cardelli, Clayton and Mathis attenuation using the
    formulae from Calzetti et al. (2000) and Leitherer et al. (2002).

    The attenuation can be computed on the whole spectrum or on a specific
    contribution and is added to the SED as a negative contribution.

    """

    parameter_list = OrderedDict([
        ("E_BVs_young", (
            "float",
            "E(B-V)*, the colour excess of the stellar continuum light for "
            "the young population.",
            None
        )),
        ("E_BVs_old_factor", (
            "float",
            "Reduction factor for the E(B-V)* of the old population compared "
            "to the young one (<1).",
            None
        )),
        ("young_contribution_name", (
            "string",
            "Name of the contribution containing the spectrum of the "
            "young population.",
            "ssp_young"
        )),
        ("old_contribution_name", (
            "string",
            "Name of the contribution containing the spectrum of the "
            "old population. If it is set to 'None', only one population "
            "is considered.",
            "ssp_old"
        )),
        ("uv_bump_wavelength", (
            "float",
            "Central wavelength of the UV bump in nm.",
            217.5
        )),
        ("uv_bump_width", (
            "float",
            "Width (FWHM) of the UV bump in nm.",
            None
        )),
        ("uv_bump_amplitude", (
            "float",
            "Amplitude of the UV bump in nm.",
            None
        )),
        ("powerlaw_slope", (
            "float",
            "Slope delta of the power law modifying the attenuation curve.",
            None
        )),
        ("filters", (
            "string",
            "Filters for which the attenuation will be computed and added to "
            "the SED information dictionary. You can give several filter "
            "names separated by a & (don't use commas).",
            "V_B90 & FUV"
        ))
    ])

    out_parameter_list = OrderedDict([
        ("E_BVs_young", "E(B-V)*, the colour excess of the stellar continuum "
                        "light for the young population."),
        ("E_BVs_old", "E(B-V)*, the colour excess of the stellar "
                      "continuum light for the old population."),
        ("attenuation_young", "Amount of luminosity attenuated from the "
                              "young population in W."),
        ("E_BVs_old_factor", "Ratio of the old population E(B-V)* to the "
                             "young one."),
        ("attenuation_old", "Amount of luminosity attenuated from the "
                            "old population in W."),
        ("attenuation", "Total amount of luminosity attenuated in W."),
        ("uv_bump_wavelength", "Central wavelength of UV bump in nm."),
        ("uv_bump_width", "Width of the UV bump in nm."),
        ("uv_bump_amplitude", "Amplitude of the UV bump in nm."),
        ("powerlaw_slope", "Slope of the power law."),
        ("FILTER_attenuation", "Attenuation in the FILTER filter.")
    ])

    def _init_code(self):
        """Get the filters from the database"""
        filter_list = [item.strip() for item in
                       self.parameters["filters"].split("&")]
        self.filters = {}
        with Database() as base:
            for filter_name in filter_list:
                self.filters[filter_name] = base.get_filter(filter_name)

    def process(self, sed):
        """Add the CCM dust attenuation to the SED.

        Parameters
        ----------
        sed : pcigale.sed.SED object

        """

        wavelength = sed.wavelength_grid
        ebvs_young = float(self.parameters["E_BVs_young"])
        ebvs_old = float(self.parameters["E_BVs_old_factor"]) * ebvs_young
        young_contrib = self.parameters["young_contribution_name"]
        old_contrib = self.parameters["old_contribution_name"]
        uv_bump_wavelength = float(self.parameters["uv_bump_wavelength"])
        uv_bump_width = float(self.parameters["uv_bump_wavelength"])
        uv_bump_amplitude = float(self.parameters["uv_bump_amplitude"])
        powerlaw_slope = float(self.parameters["powerlaw_slope"])
        filters = self.filters

        # Fλ fluxes (only from continuum) in each filter before attenuation.
        flux_noatt = {}
        for filter_name, filter_ in filters.items():
            flux_noatt[filter_name] = sed.compute_fnu(
                filter_.trans_table,
                filter_.effective_wavelength,
                add_line_fluxes=False)

        # Compute attenuation curve
        sel_attenuation = a_vs_ebv(wavelength, uv_bump_wavelength,
                                   uv_bump_width, uv_bump_amplitude,
                                   powerlaw_slope)

        # Young population attenuation
        luminosity = sed.get_lumin_contribution(young_contrib)
        attenuated_luminosity = luminosity * 10 ** (ebvs_young
                                                    * sel_attenuation / -2.5)
        attenuated_luminosity[wavelength < 91.2] = 0
        attenuation_spectrum = attenuated_luminosity - luminosity
        # We integrate the amount of luminosity attenuated (-1 because the
        # spectrum is negative).
        attenuation_young = -1 * np.trapz(attenuation_spectrum, wavelength)

        sed.add_module(self.name, self.parameters)
        sed.add_info("E_BVs_young" + self.postfix, ebvs_young)
        sed.add_info("attenuation_young" + self.postfix, attenuation_young)
        sed.add_contribution("attenuation_young" + self.postfix,
                             wavelength, attenuation_spectrum)

        # Old population (if any) attenuation
        if old_contrib:
            luminosity = sed.get_lumin_contribution(old_contrib)
            attenuated_luminosity = (luminosity *
                                     10 ** (ebvs_old
                                            * sel_attenuation / -2.5))
            attenuated_luminosity[wavelength < 91.2] = 0
            attenuation_spectrum = attenuated_luminosity - luminosity
            attenuation_old = -1 * np.trapz(attenuation_spectrum, wavelength)

            sed.add_info("E_BVs_old" + self.postfix, ebvs_old)
            sed.add_info("E_BVs_old_factor" + self.postfix,
                         self.parameters["E_BVs_old_factor"])
            sed.add_info("attenuation_old" + self.postfix, attenuation_old)
            sed.add_contribution("attenuation_old" + self.postfix,
                                 wavelength, attenuation_spectrum)
        else:
            attenuation_old = 0

        # Total attenuation (we don't take into account the energy attenuated
        # in the spectral lines)
        sed.add_info("attenuation" + self.postfix,
                     attenuation_young + attenuation_old)

        # Fλ fluxes (only from continuum) in each filter after attenuation.
        flux_att = {}
        for filter_name, filter_ in filters.items():
            flux_att[filter_name] = sed.compute_fnu(
                filter_.trans_table,
                filter_.effective_wavelength,
                add_line_fluxes=False)

        # Attenuation in each filter
        for filter_name in filters:
            sed.add_info(filter_name + "_attenuation" + self.postfix,
                         -2.5 * np.log10(flux_att[filter_name] /
                                         flux_noatt[filter_name]))

# CreationModule to be returned by get_module
Module = CalzLeit
