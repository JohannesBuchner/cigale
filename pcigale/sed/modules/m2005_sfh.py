# -*- coding: utf-8 -*-
"""
Copyright (C) 2012 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""


import numpy as np
from . import common
from pcigale.data import Database


class Module(common.SEDCreationModule):
    """Population synthesis based on Maraston (2005)

    Implements the population syntesis based on the SSP described in Maraston,
    2005, MNRAS, 362, 799.

    """

    parametre_list = {
        'imf': (
            'string',
            None,
            "Initial mass function: ss (Salpeter) or kr (Krupa)",
            None
        ),
        'metallicity': (
            'float',
            None,
            "Z/H abundance of heavy elements with respect to hydrogen, "
            "normalised to solar values.",
            None
        ),
        'sfr': (
            'numpy array of floats',
            "Solar mass per year",
            "sfr is the star formation rate vector (in solar mass per year) "
            "based on the age grid used for the Maraston 2005 SSP, i.e. "
            "np.arange(1e-3,13.701,1e-3).",
            None
        ),
        'oldest_age': (
            'float',
            'Gyr',
            "Age of the oldest stars in the galaxy.",
            None
        ),
        'separation_age': (
            'float',
            'Gyr',
            "Age of the separation between the young and the old star "
            "populations. The default value in 10^7 years (0.01 Gyr). Set "
            "to 0 not to differentiate ages.",
            0.01
        )
    }

    def _process(self, sed, parametres):
        """Add the convolution of a Maraston 2005 SSP to the SED

        Parametres
        ----------
        sed  : pcigale.sed.SED
            SED object.
        parametres : dictionnary
            Dictionnary containing the parametres

        """

        imf = self.parametres['imf']
        metallicity = self.parametres['metallicity']
        sfr = np.copy(self.parametres['sfr'])
        oldest_age = self.parametres['oldest_age']
        separation_age = self.parametres['separation_age']

        # First, we take the SSP out of the database.
        database = Database()
        ssp = database.get_ssp_m2005(imf, metallicity)
        database.session.close_all()
        # We check that the sfr array length corresponds to the one of the
        # SSP age grid.
        if len(ssp.age_grid) != len(sfr):
            raise ValueError("The SFR array must be based on the same age "
                             "grid as the Maraston 2005 SSPs.")

        # We limit the star formation history to the age of the oldest stars
        # and make the distinction between the old a the young populations.

        # Index of the age in sfr[0] nearest to the age of the oldest stars.
        idx_max = np.abs(ssp.age_grid - oldest_age).argmin()
        # Index of the age in sfr[0] separating the old and the young stars.
        idx_sep = np.abs(ssp.age_grid - oldest_age + separation_age).argmin()

        # The age of the galaxy corresponds to the age of its oldest stars
        # and is taken on the age grid.
        galaxy_age = ssp.age_grid[idx_max]

        # We normalise the SFH to have 1 solar mass formed at the age of the
        # galaxy.
        sfr = sfr / np.trapz(sfr[0:idx_max] * 1e9,
                             ssp.age_grid[0:idx_max])

        # First, we process the old population
        old_sfr = np.copy(sfr)
        old_sfr[idx_sep:] = 0
        old_masses, old_spectrum = ssp.convolve(old_sfr, galaxy_age)

        # If there is a separation between young and old galaxies, we process
        # the young population.
        if separation_age:
            sed.add_info('old_young_separation_age', separation_age)
            young_sfr = np.copy(sfr)
            young_sfr[:idx_sep] = 0
            young_sfr[idx_max:] = 0
            young_masses, young_spectrum = ssp.convolve(young_sfr, galaxy_age)

        # The way we enter the information in the SED object is slightly
        # different, depending on if there is the age separation or not.
        if not separation_age:
            sed.add_component(
                'm2005_sfh',
                parametres,
                'm2005_sfh',
                ssp.wavelength_grid,
                old_spectrum,
                {
                    'stellar_mass_total': old_masses[0],
                    'stellar_mass_alive': old_masses[1],
                    'stellar_mass_white_dwarf': old_masses[2],
                    'stellar_mass_neutrino': old_masses[3],
                    'stellar_mass_black_hole': old_masses[4],
                    'stellar_mass_turn_off': old_masses[5]
                }
            )
        else:
            sed.add_info('stellar_mass_total',
                         old_masses[0] + young_masses[0])
            sed.add_info('stellar_mass_alive',
                         old_masses[1] + young_masses[1])
            sed.add_info('stellar_mass_white_dwarf',
                         old_masses[2] + young_masses[2])
            sed.add_info('stellar_mass_neutrino',
                         old_masses[3] + young_masses[3])
            sed.add_info('stellar_mass_black_hole',
                         old_masses[4] + young_masses[4])
            sed.add_info('stellar_mass_turn_off',
                         old_masses[5] + young_masses[5])
            sed.add_component(
                'm2005_sfh',
                parametres,
                'm2005_sfh_old',
                ssp.wavelength_grid,
                old_spectrum,
                {
                    'stellar_mass_total_old': old_masses[0],
                    'stellar_mass_alive_old': old_masses[1],
                    'stellar_mass_white_dwarf_old': old_masses[2],
                    'stellar_mass_neutrino_old': old_masses[3],
                    'stellar_mass_black_hole_old': old_masses[4],
                    'stellar_mass_turn_off_old': old_masses[5]
                }
            )
            sed.add_component(
                'm2005_sfh',
                parametres,
                'm2005_sfh_young',
                ssp.wavelength_grid,
                young_spectrum,
                {
                    'stellar_mass_total_young': young_masses[0],
                    'stellar_mass_alive_young': young_masses[1],
                    'stellar_mass_white_dwarf_young': young_masses[2],
                    'stellar_mass_neutrino_young': young_masses[3],
                    'stellar_mass_black_hole_young': young_masses[4],
                    'stellar_mass_turn_off_young': young_masses[5]
                }
            )
