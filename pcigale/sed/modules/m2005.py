# -*- coding: utf-8 -*-
"""
Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""


import numpy as np
from . import common
from ...data import Database

# Time lapse used to compute the average star formation rate. We use a
# constant to keep it easily changeable for advanced user while limiting the
# number of parameters. The value is in Myr.
AV_LAPSE = 100


class Module(common.SEDCreationModule):
    """Module computing the Star Formation History contribution based on the
    Maraston (2005) models.

    Implements the population synthesis based on the SSP described in Maraston,
    2005, MNRAS, 362, 799.

    Information added to the SED:
        - imf, metallicity, galaxy_age
        - sfr: star formation rate normalised to 1 solar mass formed at the
              age of the galaxy.
        - average_sfr: SFR averaged on the last 100 Myr of the galaxy history
        - mass_total, mass_alive, mass_white_dwarf,mass_neutrino,
          mass_black_hole, mass_turn_off : stellar masses in solar mass.
        - age: age of the oldest stars in the galaxy.
        - old_young_separation_age: age (in Myr) separating the young and the
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
            "Initial mass function: salp (Salpeter) or krou (Krupa)",
            None
        ),
        'metallicity': (
            'float',
            "Metallicity Z.",
            None
        ),
        'separation_age': (
            'integer',
            "Age [Myr] of the separation between the young and the old star "
            "populations. The default value in 10^7 years (10 Myr). Set to "
            "0 not to differentiate ages (only an old population).",
            10
        )
    }

    out_parameter_list = {
        'sfr': 'Instantaneous Star Formation Rate in solar mass per year, '
               'at the age of the galaxy.',
        'average_sfr': 'Average SFR in the last 100 Myr (default) of the '
                       'galaxy history.',
        'mass_total': 'Total stellar mass of the galaxy in solar mass.',
        'mass_alive': 'Mass of alive stars in solar mass.',
        'mass_white_dwarf': 'Mass of white dwarf stars in solar mass.',
        'mass_neutrino': 'Mass of neutrino stars in solar mass.',
        'mass_black_hole': 'Mass of black holes in solar mass.',
        'mass_turn_off': 'Mass in the turn-off in solar mass.',
        'old_young_separation_age': 'Age (in Myr) separating the old and '
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

    def _init_code(self):
        """Read the SSP from the database."""
        imf = self.parameters["imf"]
        metallicity = float(self.parameters["metallicity"])
        database = Database()
        self.ssp = database.get_ssp_m2005(imf, metallicity)
        database.session.close_all()

    def process(self, sed):
        """Add the convolution of a Maraston 2005 SSP to the SED

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
        young_masses, young_spectrum = ssp.convolve(sfh_time, young_sfh)

        # Then, we process the old population. If the SFH is shorter than the
        # separation age then all the arrays will consist only of 0.
        old_sfh = np.copy(sfh_sfr)
        old_sfh[sfh_age <= separation_age] = 0
        old_masses, old_spectrum = ssp.convolve(sfh_time, old_sfh)

        # SFR of the galaxy
        sfr = sfh_sfr[len(sfh_sfr) - 1]

        # Average SFR on the last AV_LAPSE Myr of its history
        average_sfr = np.mean(sfh_sfr[sfh_age <= AV_LAPSE])

        # Base name for adding information to the SED.
        name = self.name or 'm2005_sfh'

        sed.add_module(name, self.parameters)

        sed.add_info('imf', imf)
        sed.add_info('metallicity', metallicity)
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
