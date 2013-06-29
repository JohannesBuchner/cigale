# -*- coding: utf-8 -*-
"""Copyright (C) 2013 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

This module implements the attenuation base on a power law as defined
in Charlot and Fall (2000) with a UV bump added.

"""


import numpy as np
from . import common
from ...data import Database


def power_law(wavelength, delta):
    """Compute the power law (λ / λv)^δ

    Parameters
    ----------
    wavelength : array of float
        Wavelength grid in nm.
    delta : float
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


def alambda_av(wavelength, delta, bump_wave, bump_width, bump_ampl):
    """Compute the complete attenuation curve A(λ)/Av

    The continuum is a power law (λ / λv) ** δ to with is added a UV bump.

    Parameters
    ----------
    wavelength : array of floats
        The wavelength grid (in nm) to compute the attenuation curve on.
    delta : float
        Slope of the power law.
    bump_wave : float
        Central wavelength (in nm) of the UV bump.
    bump_width : float
        Width (FWHM, in nm) of the UV bump.
    bump_ampl : float
        Amplitude of the UV bump.

    Returns
    -------
    attenuation : array of floats
        The A(λ)/Av attenuation at each wavelength of the grid.

    """
    wave = np.array(wavelength)

    attenuation = power_law(wave, delta)
    attenuation = attenuation + uv_bump(wavelength, bump_wave,
                                        bump_width, bump_ampl)

    return attenuation


class Module(common.SEDCreationModule):
    """Add CCM dust attenuation based on Charlot and Fall (2000) power law.
    """

    parameter_list = {
        "Av_young": (
            "float",
            "V-band attenuation of the young population.",
            None
        ),
        "Av_old_factor": (
            "float",
            "Reduction factor for the V-band attenuation of the old "
            "population compared to the young one (<1).",
            None
        ),
        "young_contribution_name": (
            "string",
            "Name of the contribution containing the spectrum of the "
            "young population.",
            "m2005_young"
        ),
        "old_contribution_name": (
            "string",
            "Name of the contribution containing the spectrum of the "
            "old population. If it is set to 'None', only one population "
            "is considered.",
            "m2005_old"
        ),
        "uv_bump_wavelength": (
            "float",
            "Central wavelength of the UV bump in nm.",
            217.5
        ),
        "uv_bump_width": (
            "float",
            "Width (FWHM) of the UV bump in nm.",
            None
        ),
        "uv_bump_amplitude": (
            "float",
            "Amplitude of the UV bump in nm.",
            None
        ),
        "powerlaw_slope": (
            "float",
            "Slope delta of the power law continuum.",
            -0.7
        ),
        "filters": (
            "list of strings",
            "List of the filters for which the attenuation will be computed.",
            ['V_B90', 'FUV']
        )
    }

    out_parameter_list = {
        "NAME_Av_young": "V-band attenuation of the young population.",
        "NAME_Av_old": "V-band attenuation of the old population.",
        "NAME_attenuation_young": "Amount of luminosity attenuated from the "
                                  "young population in W.",
        "NAME_Av_old_factor": "Reduction factor for the V-band attenuation "
                              "of  the old population compared to the young "
                              "one (<1).",
        "NAME_attenuation_old": "Amount of luminosity attenuated from the "
                                "old population in W.",
        "NAME_attenuation": "Total amount of luminosity attenuated in W.",
        "NAME_uv_bump_wavelength": "Central wavelength of UV bump in nm.",
        "NAME_uv_bump_width": "Width of the UV bump in nm.",
        "NAME_uv_bump_amplitude": "Amplitude of the UV bump in nm.",
        "NAME_powerlaw_slope": "Slope of the power law.",
        "NAME_FILTER": "Attenuation in the FILTER filter.",
    }

    def _init_code(self):
        """Get the filters from the database"""
        self.filters = {}
        base = Database()
        for filter_name in self.parameters["filters"]:
            self.filters[filter_name] = base.get_filter(filter_name)
        base.close()

    def process(self, sed):
        """Add the CCM dust attenuation to the SED.

        Parameters
        ----------
        sed : pcigale.sed.SED object

        """

        # Base name for adding information to the SED.
        name = self.name or 'dustatt_powerlaw_'

        wavelength = sed.wavelength_grid
        av_young = float(self.parameters["Av_young"])
        av_old = float(self.parameters["Av_old_factor"] * av_young)
        young_contrib = self.parameters["young_contribution_name"]
        old_contrib = self.parameters["old_contribution_name"]
        uv_bump_wavelength = float(self.parameters["uv_bump_wavelength"])
        uv_bump_width = float(self.parameters["uv_bump_wavelength"])
        uv_bump_amplitude = float(self.parameters["uv_bump_amplitude"])
        powerlaw_slope = float(self.parameters["powerlaw_slope"])
        filters = self.filters

        # Fλ fluxes in each filter before attenuation.
        flux_noatt = {}
        for filter_name, filter_ in filters.items():
            flux_noatt[filter_name] = sed.compute_fnu(
                filter_.trans_table,
                filter_.effective_wavelength)

        # Compute attenuation curve
        sel_attenuation = alambda_av(wavelength, powerlaw_slope,
                                     uv_bump_wavelength, uv_bump_width,
                                     uv_bump_amplitude)

        # Young population attenuation
        luminosity = sed.get_lumin_contribution(young_contrib)
        attenuated_luminosity = luminosity * 10 ** (av_young
                                                    * sel_attenuation / -2.5)
        attenuation_spectrum = attenuated_luminosity - luminosity
        # We integrate the amount of luminosity attenuated (-1 because the
        # spectrum is negative).
        attenuation_young = -1 * np.trapz(attenuation_spectrum, wavelength)

        sed.add_module(name, self.parameters)
        sed.add_info(name + "_Av_young", av_young)
        sed.add_info(name + "_attenuation_young", attenuation_young)
        sed.add_contribution(name + "_young", wavelength, attenuation_spectrum)

        # Old population (if any) attenuation
        if old_contrib:
            luminosity = sed.get_lumin_contribution(old_contrib)
            attenuated_luminosity = (luminosity *
                                     10 ** (av_old
                                            * sel_attenuation / -2.5))
            attenuation_spectrum = attenuated_luminosity - luminosity
            attenuation_old = -1 * np.trapz(attenuation_spectrum, wavelength)

            sed.add_info(name + "_Av_old", av_old)
            sed.add_info(name + "_Av_old_factor",
                         self.parameters["Av_old_factor"])
            sed.add_info(name + "_attenuation_old", attenuation_old)
            sed.add_contribution(name + "_old",
                                 wavelength, attenuation_spectrum)
        else:
            attenuation_old = 0

        # Total attenuation
        sed.add_info(name + "_attenuation",
                     attenuation_young + attenuation_old)

        # Fλ fluxes in each filter after attenuation.
        flux_att = {}
        for filter_name, filter_ in filters.items():
            flux_att[filter_name] = sed.compute_fnu(
                filter_.trans_table,
                filter_.effective_wavelength)

        # Attenuation in each filter
        for filter_name in filters:
            sed.add_info(name + "_" + filter_name,
                         -2.5 * np.log(flux_att[filter_name] /
                                       flux_noatt[filter_name]))
