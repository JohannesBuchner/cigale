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

from collections import OrderedDict

import numpy as np

from . import SedModule
from ..data import Database


class BC03(SedModule):
    """Bruzual and Charlot (2003) stellar emission module

    This SED creation module convolves the SED star formation history with a
    Bruzual and Charlot (2003) single stellar population to add a stellar
    component to the SED.
    """

    parameter_list = OrderedDict([
        ("imf", (
            "cigale_list(dtype=int, options=0. & 1.)",
            "Initial mass function: 0 (Salpeter) or 1 (Chabrier).",
            0
        )),
        ("metallicity", (
            "cigale_list(options=0.0001 & 0.0004 & 0.004 & 0.008 & 0.02 & "
            "0.05)",
            "Metalicity. Possible values are: 0.0001, 0.0004, 0.004, 0.008, "
            "0.02, 0.05.",
            0.02
        )),
        ("separation_age", (
            "cigale_list(dtype=int, minvalue=0)",
            "Age [Myr] of the separation between the young and the old star "
            "populations. The default value in 10^7 years (10 Myr). Set "
            "to 0 not to differentiate ages (only an old population).",
            10
        ))
    ])

    def _init_code(self):
        """Read the SSP from the database."""
        self.imf = int(self.parameters["imf"])
        self.metallicity = float(self.parameters["metallicity"])
        self.separation_age = int(self.parameters["separation_age"])

        with Database() as database:
            if self.imf == 0:
                self.ssp = database.get_bc03('salp', self.metallicity)
            elif self.imf == 1:
                self.ssp = database.get_bc03('chab', self.metallicity)
            else:
                raise Exception("IMF #{} unknown".format(self.imf))

    def process(self, sed):
        """Add the convolution of a Bruzual and Charlot SSP to the SED

        Parameters
        ----------
        sed: pcigale.sed.SED
            SED object.

        """
        # First, we process the young population (age lower than the
        # separation age.)
        young_wave, young_lumin, young_info = self.ssp.convolve(
            sed.sfh[-self.separation_age:])

        # Then, we process the old population. If the SFH is shorter than the
        # separation age then all the arrays will consist only of 0.
        old_sfh = np.copy(sed.sfh)
        old_sfh[-self.separation_age:] = 0.
        old_wave, old_lumin, old_info = self.ssp.convolve(old_sfh)

        # We compute the Lyman continuum luminosity as it is important to
        # compute the energy absorbed by the dust before ionising gas.
        w = np.where(young_wave <= 91.1)
        lum_ly_young = np.trapz(young_lumin[w], young_wave[w])
        lum_ly_old = np.trapz(old_lumin[w], old_wave[w])

        sed.add_module(self.name, self.parameters)

        sed.add_info("stellar.imf", self.imf)
        sed.add_info("stellar.metallicity", self.metallicity)
        sed.add_info("stellar.old_young_separation_age", self.separation_age)

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

# SedModule to be returned by get_module
Module = BC03
