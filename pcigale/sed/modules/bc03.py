# -*- coding: utf-8 -*-
"""
Copyright (C) 2013 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""


import numpy as np
from . import common
from pcigale.data import Database

# Time lapse used to compute the average star formation rate. We use a
# constant to keep it easily changeable for advanced user while limiting the
# number of parameters. The value is in Myr.
AV_LAPSE = 100


class Module(common.SEDCreationModule):
    """Module computing the Star Formation History contribution bases on the
    Bruzual and Charlot (2003) models.
    """

    parameter_list = {
        "imf": (
            "string",
            "Initial mass function: salp (Salpeter) or chab (Chabrier)",
            None
        ),
        "metallicity": (
            "float",
            "Mettalicity, 0.02 for Solar metallicity.",
            None
        ),
        "separation_age": (
            "integer",
            "Age [Myr] of the separation between the young and the old star "
            "populations. The default value in 10^7 years (10 Myr). Set "
            "to 0 not to differentiate ages (only an old population).",
            10
        )
    }

    out_parameter_list = {
        "sfr": "Instantaneous Star Formation Rate in solar mass per year, "
               "at the age of the galaxy.",
        "average_sfr": "Average SFR in the last 100 Myr (default) of the "
                       "galaxy history.",
        "m_star": "Total mass in stars in Solar mass.",
        "m_gas": "Mass returned to the ISM by evolved stars in Solar mass.",
        "n_ly": "rate of H-ionizing photons in s^-1, per Solar mass "
                "of galaxy.",
        "b_4000": "Amplitude of 4000 Å break (Bruzual 2003)",
        "b4_vn": "Amplitude of 4000 Å narrow break (Balogh et al. 1999)",
        "b4_sdss": "Amplitude of 4000 Å break (Stoughton et al. 2002)",
        "b_912": "Amplitude of Lyman discontinuity"
    }

    def _init_code(self):
        """Read the SSP from the database."""
        imf = self.parameters["imf"]
        metallicity = float(self.parameters["metallicity"])
        database = Database()
        self.ssp = database.get_ssp_bc03(imf, metallicity)
        database.session.close_all()

    def _process(self, sed, parameters):
        """Add the convolution of a Bruzual and Charlot SSP to the SED

        Parameters
        ----------
        sed  : pcigale.sed.SED
            SED object.
        parameters : dictionary
            Dictionary containing the parameters

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

        # SFR of the galaxy
        sfr = sfh_sfr[len(sfh_sfr) - 1]

        # Average SFR on the last AV_LAPSE Myr of its history
        average_sfr = np.mean(sfh_sfr[sfh_age <= AV_LAPSE])

        # Base name for adding information to the SED.
        name = self.name or "bc03"

        sed.add_module(name, parameters)

        sed.add_info(name + "_imf", imf)
        sed.add_info(name + "_metallicity", metallicity)
        sed.add_info(name + '_old_young_separation_age', separation_age)

        sed.add_info(name + '_sfr', sfr, True)
        sed.add_info(name + '_average_sfr', average_sfr, True)

        sed.add_info(name + "_m_star_young", young_info["m_star"], True)
        sed.add_info(name + "_m_gas_young", young_info["m_gas"], True)
        sed.add_info(name + "_n_ly_young", young_info["n_ly"])
        sed.add_info(name + "_b_400_young", young_info["b_4000"])
        sed.add_info(name + "_b4_vn_young", young_info["b4_vn"])
        sed.add_info(name + "_b4_sdss_young", young_info["b4_sdss"])
        sed.add_info(name + "_b_912_young", young_info["b_912"])

        sed.add_info(name + "_m_star_old", old_info["m_star"], True)
        sed.add_info(name + "_m_gas_old", old_info["m_gas"], True)
        sed.add_info(name + "_n_ly_old", old_info["n_ly"])
        sed.add_info(name + "_b_400_old", old_info["b_4000"])
        sed.add_info(name + "_b4_vn_old", old_info["b4_vn"])
        sed.add_info(name + "_b4_sdss_old", old_info["b4_sdss"])
        sed.add_info(name + "_b_912_old", old_info["b_912"])

        sed.add_contribution(name + '_old',
                             old_wave,
                             old_lumin)
        sed.add_contribution(name + '_young',
                             young_wave,
                             young_lumin)
