# -*- coding: utf-8 -*-
"""
Copyright (C) 2013 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""


import numpy as np
from . import common
from ...data import Database


def k_calzetti2000(wavelength):
    """Compute the Calzetti et al. (2000) A(λ)/E(B-V)∗

    Given a wavelength grid, this function computes the selective attenuation
    A(λ)/E(V-V)∗ using the formula from Calzetti at al. (2000). This formula
    is given for wavelengths between 120 nm and 2200 nm. This function returns
    attenuation values only for wavelengths in this range.

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
    mask = (wavelength > 120) & (wavelength < 630)
    result[mask] = 2.659 * (-2.156 + 1.509e3 / wavelength[mask]
                            - 0.198e6 / wavelength[mask] ** 2
                            + 0.011e9 / wavelength[mask] ** 3) + 4.05

    # Attenuation between 630 nm and 2200 nm
    mask = (wavelength >= 630) & (wavelength < 2200)
    result[mask] = 2.659 * (-1.857 + 1.040e3 / wavelength[mask]) + 4.05

    return result


def k_leitherer2002(wavelength):
    """Compute the Leitherer et al. (2002) A(λ)/E(B-V)∗

    Given a wavelength grid, this function computes the selective attenuation
    A(λ)/E(V-V)∗ using the formula from Leitherer at al. (2002). This formula
    is given for wavelengths between 91.2 nm and 180 nm. This function returns
    attenuation values only for wavelengths in this range.

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

    mask = (wavelength > 91.2) & (wavelength < 180)
    result[mask] = (5.472 + 0.671e3 / wavelength[mask]
                    - 9.218e3 / wavelength[mask] ** 2
                    + 2.620e6 / wavelength[mask] ** 3)

    return result


class Module(common.SEDCreationModule):
    """Add CCM dust attenuation based on the Calzetti formula

    If a contribution name is given in the parameter list, the attenuation
    will be applied only to the flux of this contribution; else, it will be
    applied to the whole spectrum.

    """

    parameter_list = {
        "E_BVs_young": (
            "float",
            "E(B-V)*, the colour excess of the stellar continuum light for "
            "the young population.",
            None
        ),
        "E_BVs_old_factor": (
            "float",
            "Reduction factor for the E(B-V)* of the old population compared "
            "to the young one (<1).",
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
        "filters": (
            "list of strings",
            "List of the filters for which the attenuation will be computed.",
            ['V_B90', 'FUV']
        )
    }

    out_parameter_list = {
        "NAME_E_BVs_young": "E(B-V)*, the colour excess of the stellar "
                            "continuum light for the young population.",
        "NAME_attenuation_young": "Amount of luminosity attenuated from the "
                                  "young population in W.",
        "NAME_E_BVs_old_factor": "Ratio of the old population E(B-V)* to the "
                                 "young one.",
        "NAME_attenuation_old": "Amount of luminosity attenuated from the "
                                "old population in W.",
        "NAME_attenuation": "Total amount of luminosity attenuated in W.",
        "NAME_FILTER": "Attenuation in the FILTER filter.",
    }

    def _process(self, sed, parameters):
        """Add the CCM dust attenuation to the SED.

        Parameters
        ----------
        sed : pcigale.sed.SED object
        parameters : dictionary containing the parameters

        """

        # Base name for adding information to the SED.
        name = self.name or 'dustatt'

        wavelength = sed.wavelength_grid
        ebvs_young = parameters["E_BVs_young"]
        ebvs_old = parameters["E_BVs_old_factor"] * ebvs_young
        young_contrib = parameters["young_contribution_name"]
        old_contrib = parameters["old_contribution_name"]
        filter_list = parameters["filters"]

        # Fλ fluxes in each filter before attenuation.
        flux_noatt = {}
        base = Database()
        for filter_name in filter_list:
            filt = base.get_filter(filter_name)
            flux_noatt[filter_name] = sed.compute_fnu(
                filt.trans_table,
                filt.effective_wavelength)
        base.close()

        # Selective attenuation. We use Leitherer et al. (2002) formula
        # between 91.2 and 150 nm and Calzetti et al. (2000) formula between
        # 150 and 2200 nm. Elsewhere, we set no attenuation.
        sel_attenuation = np.zeros(len(wavelength))
        mask = (wavelength > 91.2) & (wavelength <= 150)
        sel_attenuation[mask] = k_leitherer2002(wavelength[mask])
        mask = (wavelength > 150) & (wavelength < 2200)
        sel_attenuation[mask] = k_calzetti2000(wavelength[mask])

        # Young population attenuation
        luminosity = sed.get_lumin_contribution(young_contrib)
        attenuated_luminosity = luminosity * 10**(ebvs_young
                                                  * sel_attenuation / -2.5)
        attenuation_spectrum = attenuated_luminosity - luminosity
        # We integrate the amount of luminosity attenuated (-1 because the
        # spectrum is negative).
        attenuation_young = -1 * np.trapz(attenuation_spectrum, wavelength)

        sed.add_module(name, parameters)
        sed.add_info(name + "_E_BVs_young", ebvs_young)
        sed.add_info(name + "_attenuation_young", attenuation_young)
        sed.add_contribution(name + "_young", wavelength, attenuation_spectrum)

        # Old population (if any) attenuation
        if old_contrib:
            luminosity = sed.get_lumin_contribution(old_contrib)
            attenuated_luminosity = luminosity * 10**(ebvs_old
                                                      * sel_attenuation / -2.5)
            attenuation_spectrum = attenuated_luminosity - luminosity
            attenuation_old = -1 * np.trapz(attenuation_spectrum, wavelength)

            sed.add_info(name + "_E_BVs_old", ebvs_old)
            sed.add_info(name + "_E_BVs_old_factor",
                         parameters["E_BVs_old_factor"])
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
        base = Database()
        for filter_name in filter_list:
            filt = base.get_filter(filter_name)
            flux_att[filter_name] = sed.compute_fnu(filt.trans_table,
                                                    filt.effective_wavelength)
        base.close()

        # Attenuation in each filter
        for filter_name in filter_list:
            sed.add_info(name + "_" + filter_name,
                         -2.5 * np.log(flux_att[filter_name] /
                                       flux_noatt[filter_name]))
