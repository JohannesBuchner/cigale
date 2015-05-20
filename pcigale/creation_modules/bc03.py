# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

"""
Bruzual and Charlot (2003) stellar emission module
==================================================

This module implements the Bruzual and Charlot (2003) Single Stellar
Populations.

"""

import numpy as np
from collections import OrderedDict
from . import CreationModule
from ..data import Database


class BC03(CreationModule):
    """Bruzual and Charlot (2003) stellar emission module

    This SED creation module convolves the SED star formation history with a
    Bruzual and Charlot (2003) single stellar population to add a stellar
    component to the SED.
    """

    parameter_list = OrderedDict([
        ("imf", (
            "int",
            "Initial mass function: 0 (Salpeter) or 1 (Chabrier).",
            0
        )),
        ("metallicity", (
            "float",
            "Metalicity. Possible values are: 0.0001, 0.0004, 0.004, 0.008, "
            "0.02, 0.05.",
            0.02
        )),
        ("separation_age", (
            "int",
            "Age [Myr] of the separation between the young and the old star "
            "populations. The default value in 10^7 years (10 Myr). Set "
            "to 0 not to differentiate ages (only an old population).",
            10
        ))
    ])

    out_parameter_list = OrderedDict([
        ("sfr", "Instantaneous Star Formation Rate in solar mass per year, "
                "at the age of the galaxy."),
        ('sfr10Myrs', 'Average SFR in the last 10 Myr (default) of the '
                        'galaxy history.'),
        ('sfr100Myrs', 'Average SFR in the last 100 Myr (default) of the '
                        'galaxy history.'),
        ("ssp_m_star", "Total mass in stars in Solar mass."),
        ("ssp_m_gas", "Mass returned to the ISM by evolved stars in Solar "
                      "mass."),
        ("ssp_n_ly", "rate of H-ionizing photons in s^-1, per Solar mass "
                     "of galaxy."),
        ("ssp_b_4000", "Amplitude of 4000 Å break (Bruzual 2003)"),
        ("ssp_b4_vn", "Amplitude of 4000 Å narrow break (Balogh et al. 1999)"),
        ("ssp_b4_sdss", "Amplitude of 4000 Å break (Stoughton et al. 2002)"),
        ("ssp_b_912", "Amplitude of Lyman discontinuity")
    ])

    def _init_code(self):
        """Read the SSP from the database."""
        if self.parameters["imf"] == 0:
            imf = 'salp'
        elif self.parameters["imf"] == 1:
            imf = 'chab'
        metallicity = float(self.parameters["metallicity"])
        with Database() as database:
            self.ssp = database.get_bc03(imf, metallicity)

    def process(self, sed):
        """Add the convolution of a Bruzual and Charlot SSP to the SED

        Parameters
        ----------
        sed: pcigale.sed.SED
            SED object.

        """
        imf = self.parameters["imf"]
        metallicity = float(self.parameters["metallicity"])
        separation_age = int(self.parameters["separation_age"])
        sfh_time, sfh_sfr = sed.sfh
        ssp = self.ssp

        # Age of the galaxy at each time of the SFH
        sfh_age = np.max(sfh_time) - sfh_time

        # First, we process the young population (age lower than the
        # separation age.)
        young_sfh = np.copy(sfh_sfr)
        young_sfh[sfh_age > separation_age] = 0
        young_wave, young_lumin, young_info = ssp.convolve(sfh_time, young_sfh)

        # Then, we process the old population. If the SFH is shorter than the
        # separation age then all the arrays will consist only of 0.
        old_sfh = np.copy(sfh_sfr)
        old_sfh[sfh_age <= separation_age] = 0
        old_wave, old_lumin, old_info = ssp.convolve(sfh_time, old_sfh)

        # We compute the Lyman continuum luminosity as it is important to
        # compute the energy absorbed by the dust before ionising gas.
        w = np.where(young_wave <= 91.1)
        lum_ly_young = np.trapz(young_lumin[w], young_wave[w])
        lum_ly_old = np.trapz(old_lumin[w], old_wave[w])

        sed.add_module(self.name, self.parameters)

        sed.add_info("stellar.imf", imf)
        sed.add_info("stellar.metallicity", metallicity)
        sed.add_info("stellar.old_young_separation_age", separation_age)

        sed.add_info("stellar.m_star_young", young_info["m_star"], True)
        sed.add_info("stellar.m_gas_young", young_info["m_gas"], True)
        sed.add_info("stellar.n_ly_young", young_info["n_ly"], True)
        sed.add_info("stellar.lum_ly_young", lum_ly_young, True)
        sed.add_info("stellar.b_400_young", young_info["b_4000"])
        sed.add_info("stellar.b4_vn_young", young_info["b4_vn"])
        sed.add_info("stellar.b4_sdss_young", young_info["b4_sdss"])
        sed.add_info("stellar.b_912_young", young_info["b_912"])

        sed.add_info("stellar.m_star_old", old_info["m_star"], True)
        sed.add_info("stellar.m_gas_old", old_info["m_gas"], True)
        sed.add_info("stellar.n_ly_old", old_info["n_ly"], True)
        sed.add_info("stellar.lum_ly_old", lum_ly_old, True)
        sed.add_info("stellar.b_400_old", old_info["b_4000"])
        sed.add_info("stellar.b4_vn_old", old_info["b4_vn"])
        sed.add_info("stellar.b4_sdss_old", old_info["b4_sdss"])
        sed.add_info("stellar.b_912_old", old_info["b_912"])

        sed.add_info("stellar.m_star",
                     young_info["m_star"] + old_info["m_star"],
                     True)
        sed.add_info("stellar.m_gas",
                     young_info["m_gas"] + old_info["m_gas"],
                     True)

        sed.add_contribution("stellar.old", old_wave, old_lumin)
        sed.add_contribution("stellar.young", young_wave, young_lumin)

# CreationModule to be returned by get_module
Module = BC03
