# -*- coding: utf-8 -*-
# Copyright (C) 2013-2015 Centre de données Astrophysiques de Marseille
# Copyright (C) 2014 Laboratoire d'Astrophysique de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly, Véronique Buat, Denis Burgarella, Barabara Lo Faro

"""
Double power law attenuation module
===================================

This module implements an attenuation law combining the birth cloud (BC)
attenuation and the interstellar medium (ISM) attenuation, each one modelled by
a power law. The young star emission is attenuated by the BC and the ISM
attenuations whereas the old star emission is only affected by the ISM.

Parameters available for analysis
---------------------------------

- attenuation.Av_BC: Av attenuation in the birth clouds
- attenuation.slope_BC: Slope of the power law in the birth clouds
- attenuation.BC_to_ISM_factor: Av in the ISM / Av in birth clouds
- attenuation.slope_ISM: Slope of the power law in the ISM
- attenuation.<NAME>: amount of total attenuation in the luminosity
    contribution <NAME>
- attenuation.<FILTER>: total attenuation in the filter
"""

from collections import OrderedDict

import numpy as np

from . import SedModule


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


def alambda_av(wavelengths, delta, delta_sec=None, factor=None):
    """Compute the complete attenuation curve A(λ)/Av

    The attenuation curve is a power law (λ / λv) ** δ. If a factor and
    a second delta are given, another power law is added multiplied by the
    given factor.  For instance, for the young star population the delta will
    be the slope of the birth cloud attenuation and the delta_sec will be the
    slope of the ISM attenuation.

    The Lyman continuum is not attenuated.

    Parameters
    ----------
    wavelengths: array of floats
        The wavelength grid (in nm) to compute the attenuation curve on.
    delta: float
        Slope of the main power law.
    delta_sec: float
        Slope of the secondary power law.
    factor: float
        Factor by which the secondary power law is multiplied before being
        added to the main one.

    Returns
    -------
    attenuation: array of floats
        The A(λ)/Av attenuation at each wavelength of the grid.

    """
    wave = np.array(wavelengths)

    attenuation = power_law(wave, delta)

    if factor:
        attenuation += factor * power_law(wave, delta_sec)

    # Lyman continuum not attenuated.
    attenuation[wave <= 91.2] = 0.

    return attenuation


class TwoPowerLawAtt(SedModule):
    """Two power laws attenuation module

    Attenuation module combining the birth cloud (BC) attenuation and the
    interstellar medium (ISM) one.

    The attenuation can be computed on the whole spectrum or on a specific
    contribution and is added to the SED as a negative contribution.

    """

    parameter_list = OrderedDict([
        ("Av_BC", (
            "cigale_list(minvalue=0)",
            "V-band attenuation in the birth clouds.",
            1.
        )),
        ("slope_BC", (
            "cigale_list()",
            "Power law slope of the attenuation in the birth clouds.",
            -1.3
        )),
        ("BC_to_ISM_factor", (
            "cigale_list(minvalue=0., maxvalue=1.)",
            "Av ISM / Av BC (<1).",
            0.44
        )),
        ("slope_ISM", (
            "cigale_list()",
            "Power law slope of the attenuation in the ISM.",
            -0.7
        )),
        ("filters", (
            "string()",
            "Filters for which the attenuation will be computed and added to "
            "the SED information dictionary. You can give several filter "
            "names separated by a & (don't use commas).",
            "V_B90 & FUV"
        ))
    ])

    def _init_code(self):
        self.Av_BC = float(self.parameters['Av_BC'])
        self.slope_BC = float(self.parameters['slope_BC'])
        self.BC_to_ISM_factor = float(self.parameters['BC_to_ISM_factor'])
        self.slope_ISM = float(self.parameters['slope_ISM'])
        self.filter_list = [item.strip() for item in
                            self.parameters["filters"].split("&")]

    def process(self, sed):
        """Add the dust attenuation to the SED.

        Parameters
        ----------
        sed: pcigale.sed.SED object

        """

        def ism_lumin_attenuation(wave):
            """Compute the luminosity attenuation factor for the ISM"""
            return 10 ** (self.Av_BC * self.BC_to_ISM_factor *
                          alambda_av(wave, self.slope_ISM) / -2.5)

        def bc_ism_lumin_attenuation(wave):
            """Compute the luminosity attenuation factor for ISM + BC"""
            return 10 ** (self.Av_BC *
                          alambda_av(wave, self.slope_BC, self.slope_ISM,
                                     self.BC_to_ISM_factor) / -2.5)

        wavelength = sed.wavelength_grid

        # Fλ fluxes in each filter before attenuation.
        flux_noatt = {filt: sed.compute_fnu(filt) for filt in self.filter_list}

        attenuation_total = 0.
        contribs = [contrib for contrib in sed.contribution_names if
                    'absorption' not in contrib]

        for contrib in contribs:

            age = contrib.split('.')[-1].split('_')[-1]

            luminosity = sed.get_lumin_contribution(contrib)

            # Use the ISM attenuation unless we are dealing with young
            # populations.
            if age == 'young':
                extinction_factor = bc_ism_lumin_attenuation(wavelength)
            else:
                extinction_factor = ism_lumin_attenuation(wavelength)

            attenuated_luminosity = luminosity * extinction_factor
            attenuation_spectrum = attenuated_luminosity - luminosity
            # We integrate the amount of luminosity attenuated (-1 because the
            # spectrum is negative).
            attenuation = -1 * np.trapz(attenuation_spectrum, wavelength)
            attenuation_total += attenuation

            sed.add_module(self.name, self.parameters)
            sed.add_info("attenuation." + contrib, attenuation, True)
            sed.add_contribution("attenuation." + contrib, wavelength,
                                 attenuation_spectrum)

        sed.add_info('attenuation.Av_BC', self.Av_BC)
        sed.add_info('attenuation.slope_BC', self.slope_BC)
        sed.add_info('attenuation.BC_to_ISM_factor', self.BC_to_ISM_factor)
        sed.add_info('attenuation.slope_ISM', self.slope_ISM)

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
Module = TwoPowerLawAtt
