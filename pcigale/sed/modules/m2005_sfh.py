# -*- coding: utf-8 -*-
"""
Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""


import numpy as np
from . import common
from pcigale.data import Database


# Time lapse used to compute the average star formation rate. We use a
# constant to keep it easily changeable for advanced user while limiting the
# number of parameters. The value is in Gyr.
AV_LAPSE = 0.1


class Module(common.SEDCreationModule):
    """Module computing the Star Formation History contribution bases on the
    Maraston (2005) models.

    Implements the population syntesis based on the SSP described in Maraston,
    2005, MNRAS, 362, 799.

    Information added to the SED:
        - imf, metallicity, galaxy_age
        - sfr: star formation rate normalised to 1 solar mass formed at the
              age of the galaxy.
        - average_sfr: SFR averaged on the last 0.1 Gyr of the galaxy history
        - mass_total, mass_alive, mass_white_dwarf,mass_neutrino,
          mass_black_hole, mass_turn_off : stellar masses in solar mass.
        - age: age of the oldest stars in the galaxy.
        - old_young_separation_age: age (in Gyr) separating the young and the
              old star populations (if 0, there is only one population)
        - mass_total_old, mass_alive_old, mass_white_dwarf_old,
          mass_neutrino_old, mass_black_hole_old, mass_turn_off_old: old
              star population masses.
        - mass_total_young, mass_alive_young, mass_white_dwarf_young,
          mass_neutrino_young, mass_black_hole_young, mass_turn_off_young:
              young star population masses.

    """

    parameter_list = {
        'imf': (
            'string',
            "Initial mass function: ss (Salpeter) or kr (Krupa)",
            None
        ),
        'metallicity': (
            'float',
            "Z/H abundance of heavy elements with respect to hydrogen, "
            "normalised to solar values.",
            None
        ),
        'sfh': (
            'numpy array of floats',
            "Star formation history (in solar mass per year) as an array "
            "based on the age grid used for the Maraston 2005 SSP, i.e. "
            "np.arange(1e-3,13.701,1e-3).",
            None
        ),
        'oldest_age': (
            'float',
            "Age of the oldest stars in the galaxy in Gyr.",
            None
        ),
        'separation_age': (
            'float',
            "Age [Gyr] of the separation between the young and the old star "
            "populations. The default value in 10^7 years (0.01 Gyr). Set "
            "to 0 not to differentiate ages.",
            0.01
        )
    }

    out_parameter_list = {
        'galaxy_age': 'Age (in Gyr) of the oldest stars in the galaxy.',
        'sfr': 'Instantaneous Star Formation Rate in solar mass per year, '
               'at the age of the galaxy.',
        'average_sfr': 'Average SFR in the last 0.1 Gyr (default) of the '
                       'galaxy history.',
        'mass_total': 'Total stellar mass of the galaxy in solar mass.',
        'mass_alive': 'Mass of alive stars in solar mass.',
        'mass_white_dwarf': 'Mass of white dwarf stars in solar mass.',
        'mass_neutrino': 'Mass of neutrino stars in solar mass.',
        'mass_black_hole': 'Mass of black holes in solar mass.',
        'mass_turn_off': 'Mass in the turn-off in solar mass.',
        'age': 'Age of the oldest stars in the galaxy.',
        'old_young_separation_age': 'Age (in Gyr) separating the old and '
                                    'the young star populations (0 if there '
                                    'is only one population).',
        'mass_total_old': 'Total stellar mass of the old population in solar '
                          'mass.',
        'mass_alive_old': 'Mass of alive stars in solar mass (old '
                          'population).',
        'mass_white_dwarf_old': 'Mass of white dwarf stars in solar mass '
                                '(old population).',
        'mass_neutrino_old': 'Mass of neutrino stars in solar mass '
                             '(old population).',
        'mass_black_hole_old': 'Mass of black holes in solar mass '
                               '(old population).',
        'mass_turn_off_old': 'Mass in the turn-off in solar mass '
                             '(old population).',
        'mass_total_young': 'Total stellar mass of the young population '
                            'in solar mass.',
        'mass_alive_young': 'Mass of alive stars in solar mass '
                            '(young population).',
        'mass_white_dwarf_young': 'Mass of white dwarf stars in solar mass '
                                  '(young population).',
        'mass_neutrino_young': 'Mass of neutrino stars in solar mass '
                               '(young population).',
        'mass_black_hole_young': 'Mass of black holes in solar mass '
                                 '(young population).',
        'mass_turn_off_young': 'Mass in the turn-off in solar mass '
                               '(young population).'
    }

    def _process(self, sed, parameters):
        """Add the convolution of a Maraston 2005 SSP to the SED

        Parameters
        ----------
        sed  : pcigale.sed.SED
            SED object.
        parameters : dictionary
            Dictionary containing the parameters

        """

        imf = self.parameters['imf']
        metallicity = self.parameters['metallicity']
        sfh = np.copy(self.parameters['sfh'])
        oldest_age = self.parameters['oldest_age']
        separation_age = self.parameters['separation_age']

        # First, we take the SSP out of the database.
        database = Database()
        ssp = database.get_ssp_m2005(imf, metallicity)
        database.session.close_all()
        # We check that the sfh array length corresponds to the one of the
        # SSP age grid.
        if len(ssp.time_grid) != len(sfh):
            raise ValueError("The SFH array must be based on the same age "
                             "grid as the Maraston 2005 SSPs.")

        # We limit the star formation history to the age of the oldest stars
        # and set SFR=0 for every age after (we need to keep the same age
        # grid for the computations.)
        sfh[ssp.time_grid > oldest_age] = 0

        # The age of the galaxy is taken from the age grid.
        galaxy_age = np.max(ssp.time_grid[ssp.time_grid <= oldest_age])

        # We normalise the SFH to have 1 solar mass formed at the age of the
        # galaxy.
        sfh = sfh / np.trapz(sfh * 1e9, ssp.time_grid)

        # First, we process the old population (age greater than the
        # separation age).
        old_sfh = np.copy(sfh)
        old_sfh[(galaxy_age - ssp.time_grid) < separation_age] = 0
        old_masses, old_spectrum = ssp.convolve(old_sfh, galaxy_age)

        # If there is a separation between young and old galaxies, we process
        # the young population else we set everything to 0 for the young
        # population (always having a young population can be help full e.g.
        # to probe various separation ages including 0).
        if separation_age:
            young_sfh = np.copy(sfh)
            young_sfh[(galaxy_age - ssp.time_grid) >= separation_age] = 0
            young_masses, young_spectrum = ssp.convolve(young_sfh, galaxy_age)
        else:
            young_masses = np.zeros(6)
            young_spectrum = np.zeros(len(ssp.wavelength_grid))

        # SFR of the galaxy
        sfr = sfh[ssp.time_grid == galaxy_age][0]

        # Average SFR on the last AV_LAPSE Gyr of its history
        average_sfr = np.mean(sfh[(ssp.time_grid <= oldest_age) &
                                  (ssp.time_grid >= (oldest_age - AV_LAPSE))])

        # Base name for adding information to the SED.
        name = self.name or 'm2005_sfh'

        sed.add_module(name, parameters)

        sed.add_info('imf', imf)
        sed.add_info('metallicity', metallicity)
        sed.add_info('age', oldest_age)
        sed.add_info('old_young_separation_age', separation_age)

        sed.add_info('sfr', sfr, True)
        sed.add_info('average_sfr', average_sfr, True)

        sed.add_info('mass_total_old', old_masses[0], True)
        sed.add_info('mass_alive_old', old_masses[1], True)
        sed.add_info('mass_white_dwarf_old', old_masses[2], True)
        sed.add_info('mass_neutrino_old', old_masses[3], True)
        sed.add_info('mass_black_hole_old', old_masses[4], True)
        sed.add_info('mass_turn_off_old', old_masses[5], True)

        sed.add_info('mass_total_young', young_masses[0], True)
        sed.add_info('mass_alive_young', young_masses[1], True)
        sed.add_info('mass_white_dwarf_young', young_masses[2], True)
        sed.add_info('mass_neutrino_young', young_masses[3], True)
        sed.add_info('mass_black_hole_young', young_masses[4], True)
        sed.add_info('mass_turn_off_young', young_masses[5], True)

        sed.add_info('mass_total', old_masses[0] + young_masses[0], True)
        sed.add_info('mass_alive', old_masses[1] + young_masses[1], True)
        sed.add_info('mass_white_dwarf', old_masses[2] + young_masses[2], True)
        sed.add_info('mass_neutrino', old_masses[3] + young_masses[3], True)
        sed.add_info('mass_black_hole', old_masses[4] + young_masses[4], True)
        sed.add_info('mass_turn_off', old_masses[5] + young_masses[5], True)

        sed.add_contribution(name + '_old',
                             ssp.wavelength_grid,
                             old_spectrum)
        sed.add_contribution(name + '_young',
                             ssp.wavelength_grid,
                             young_spectrum)
