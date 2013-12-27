# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly <yannick.roehlly@oamp.fr>

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
            "string",
            "Initial mass function: salp (Salpeter) or chab (Chabrier)",
            None
        )),
        ("metallicity", (
            "float",
            "Mettalicity, 0.02 for Solar metallicity.",
            None
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
        ("average_sfr", "Average SFR in the last 100 Myr (default) of the "
                        "galaxy history."),
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
        imf = self.parameters["imf"]
        metallicity = float(self.parameters["metallicity"])
        database = Database()
        self.ssp = database.get_ssp_bc03(imf, metallicity)
        database.session.close_all()

    def process(self, sed):
        """Add the convolution of a Bruzual and Charlot SSP to the SED

        Parameters
        ----------
        sed  : pcigale.sed.SED
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

        sed.add_module(self.name, self.parameters)

        sed.add_info("ssp_imf" + self.postfix, imf)
        sed.add_info("ssp_metallicity" + self.postfix, metallicity)
        sed.add_info("ssp_old_young_separation_age" + self.postfix,
                     separation_age)

        sed.add_info("ssp_m_star_young" + self.postfix,
                     young_info["m_star"], True)
        sed.add_info("ssp_m_gas_young" + self.postfix,
                     young_info["m_gas"], True)
        sed.add_info("ssp_n_ly_young" + self.postfix, young_info["n_ly"])
        sed.add_info("ssp_b_400_young" + self.postfix, young_info["b_4000"])
        sed.add_info("ssp_b4_vn_young" + self.postfix, young_info["b4_vn"])
        sed.add_info("ssp_b4_sdss_young" + self.postfix, young_info["b4_sdss"])
        sed.add_info("ssp_b_912_young" + self.postfix, young_info["b_912"])

        sed.add_info("ssp_m_star_old" + self.postfix, old_info["m_star"], True)
        sed.add_info("ssp_m_gas_old" + self.postfix, old_info["m_gas"], True)
        sed.add_info("ssp_n_ly_old" + self.postfix, old_info["n_ly"])
        sed.add_info("ssp_b_400_old" + self.postfix, old_info["b_4000"])
        sed.add_info("ssp_b4_vn_old" + self.postfix, old_info["b4_vn"])
        sed.add_info("ssp_b4_sdss_old" + self.postfix, old_info["b4_sdss"])
        sed.add_info("ssp_b_912_old" + self.postfix, old_info["b_912"])

        sed.add_contribution("ssp_old" + self.postfix,
                             old_wave,
                             old_lumin)
        sed.add_contribution("ssp_young" + self.postfix,
                             young_wave,
                             young_lumin)

# CreationModule to be returned by get_module
Module = BC03
