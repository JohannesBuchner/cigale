# -*- coding: utf-8 -*-
# Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""
Maraston (2005) stellar emission module
=======================================

This module implements the Maraston (2005) Single Stellar Populations.

"""

import numpy as np
from collections import OrderedDict
from . import CreationModule
from ..data import Database


class M2005(CreationModule):
    """Maraston (2005) stellar emission module

    This SED creation module convolves the SED star formation history with
    a Maraston (2005) single stellar population to add a stellar component to
    the SED.

    Information added to the SED:
        - imf, metallicity, galaxy_age
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

    parameter_list = OrderedDict([
        ('imf', (
            'string',
            "Initial mass function, salp (Salpeter) or krou (Krupa)",
            None
        )),
        ('metallicity', (
            'float',
            "Metallicity Z.",
            None
        )),
        ('separation_age', (
            'integer',
            "Age [Myr] of the separation between the young and the old star "
            "populations. The default value in 10^7 years (10 Myr). Set to "
            "0 not to differentiate ages (only an old population).",
            10
        ))
    ])

    out_parameter_list = OrderedDict([
        ('sfr', 'Instantaneous Star Formation Rate in solar mass per year, '
                'at the age of the galaxy.'),
        ('average_sfr', 'Average SFR in the last 100 Myr (default) of the '
                        'galaxy history.'),
        ('mass_total', 'Total stellar mass of the galaxy in solar mass.'),
        ('mass_alive', 'Mass of alive stars in solar mass.'),
        ('mass_white_dwarf', 'Mass of white dwarf stars in solar mass.'),
        ('mass_neutrino', 'Mass of neutrino stars in solar mass.'),
        ('mass_black_hole', 'Mass of black holes in solar mass.'),
        ('mass_turn_off', 'Mass in the turn-off in solar mass.'),
        ('old_young_separation_age', 'Age (in Myr) separating the old and '
                                     'the young star populations (0 if there '
                                     'is only one population).'),
        ('mass_total_old', 'Total stellar mass of the old population in solar '
                           'mass.'),
        ('mass_alive_old', 'Mass of alive stars in solar mass (old '
                           'population).'),
        ('mass_white_dwarf_old', 'Mass of white dwarf stars in solar mass '
                                 '(old population).'),
        ('mass_neutrino_old', 'Mass of neutrino stars in solar mass '
                              '(old population).'),
        ('mass_black_hole_old', 'Mass of black holes in solar mass '
                                '(old population).'),
        ('mass_turn_off_old', 'Mass in the turn-off in solar mass '
                              '(old population).'),
        ('mass_total_young', 'Total stellar mass of the young population '
                             'in solar mass.'),
        ('mass_alive_young', 'Mass of alive stars in solar mass '
                             '(young population).'),
        ('mass_white_dwarf_young', 'Mass of white dwarf stars in solar mass '
                                   '(young population).'),
        ('mass_neutrino_young', 'Mass of neutrino stars in solar mass '
                                '(young population).'),
        ('mass_black_hole_young', 'Mass of black holes in solar mass '
                                  '(young population).'),
        ('mass_turn_off_young', 'Mass in the turn-off in solar mass '
                                '(young population).')
    ])

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

        sed.add_module(self.name, self.parameters)

        sed.add_info('ssp_imf' + self.postfix, imf)
        sed.add_info('ssp_metallicity' + self.postfix, metallicity)
        sed.add_info('ssp_old_young_separation_age' + self.postfix,
                     separation_age)

        sed.add_info('ssp_mass_total_old' + self.postfix, old_masses[0], True)
        sed.add_info('ssp_mass_alive_old' + self.postfix, old_masses[1], True)
        sed.add_info('ssp_mass_white_dwarf_old' + self.postfix, old_masses[2],
                     True)
        sed.add_info('ssp_mass_neutrino_old' + self.postfix, old_masses[3],
                     True)
        sed.add_info('ssp_mass_black_hole_old' + self.postfix, old_masses[4],
                     True)
        sed.add_info('ssp_mass_turn_off_old' + self.postfix, old_masses[5],
                     True)

        sed.add_info('ssp_mass_total_young' + self.postfix, young_masses[0],
                     True)
        sed.add_info('ssp_mass_alive_young' + self.postfix, young_masses[1],
                     True)
        sed.add_info('ssp_mass_white_dwarf_young' + self.postfix,
                     young_masses[2], True)
        sed.add_info('ssp_mass_neutrino_young' + self.postfix, young_masses[3],
                     True)
        sed.add_info('ssp_mass_black_hole_young' + self.postfix,
                     young_masses[4], True)
        sed.add_info('ssp_mass_turn_off_young' + self.postfix, young_masses[5],
                     True)

        sed.add_info('ssp_mass_total' + self.postfix,
                     old_masses[0] + young_masses[0], True)
        sed.add_info('ssp_mass_alive' + self.postfix,
                     old_masses[1] + young_masses[1], True)
        sed.add_info('ssp_mass_white_dwarf' + self.postfix,
                     old_masses[2] + young_masses[2], True)
        sed.add_info('ssp_mass_neutrino' + self.postfix,
                     old_masses[3] + young_masses[3], True)
        sed.add_info('ssp_mass_black_hole' + self.postfix,
                     old_masses[4] + young_masses[4], True)
        sed.add_info('ssp_mass_turn_off' + self.postfix,
                     old_masses[5] + young_masses[5], True)

        sed.add_contribution("ssp_old" + self.postfix,
                             ssp.wavelength_grid,
                             old_spectrum)
        sed.add_contribution("ssp_young" + self.postfix,
                             ssp.wavelength_grid,
                             young_spectrum)

# CreationModule to be returned by get_module
Module = M2005
