# -*- coding: utf-8 -*-
# Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

"""
Maraston (2005) stellar emission module
=======================================

This module implements the Maraston (2005) Single Stellar Populations.

"""

from collections import OrderedDict

import numpy as np

from . import SedModule
from ..data import Database


class M2005(SedModule):
    """Maraston (2005) stellar emission module

    This SED creation module convolves the SED star formation history with
    a Maraston (2005) single stellar population to add a stellar component to
    the SED.

    Information added to the SED:
        - imf, metallicity, galaxy_age
        - mass_total, mass_alive, mass_white_dwarf,mass_neutron,
          mass_black_hole: stellar masses in solar mass.
        - age: age of the oldest stars in the galaxy.
        - old_young_separation_age: age (in Myr) separating the young and the
              old star populations (if 0, there is only one population)
        - mass_total_old, mass_alive_old, mass_white_dwarf_old,
          mass_neutron_old, mass_black_hole_old, : old star population masses.
        - mass_total_young, mass_alive_young, mass_white_dwarf_young,
          mass_neutron_young, mass_black_hole_young: young star population
              masses.

    """

    parameter_list = OrderedDict([
        ('imf', (
            'cigale_list(dtype=int, options=0. & 1.)',
            "Initial mass function: 0 (Salpeter) or 1 (Kroupa)",
            0
        )),
        ('metallicity', (
            'cigale_list(options=0.001 & 0.01 & 0.02 & 0.04)',
            "Metallicity. Possible values are: 0.001, 0.01, 0.02, 0.04.",
            0.02
        )),
        ('separation_age', (
            'cigale_list(dtype=int, minvalue=0.)',
            "Age [Myr] of the separation between the young and the old star "
            "populations. The default value in 10^7 years (10 Myr). Set to "
            "0 not to differentiate ages (only an old population).",
            10
        ))
    ])

    def _init_code(self):
        """Read the SSP from the database."""
        self.imf = int(self.parameters["imf"])
        self.metallicity = float(self.parameters["metallicity"])
        self.separation_age = int(self.parameters["separation_age"])

        if self.imf == 0:
            with Database() as database:
                self.ssp = database.get_m2005('salp', self.metallicity)
        elif self.imf == 1:
            with Database() as database:
                self.ssp = database.get_m2005('krou', self.metallicity)
        else:
            raise Exception("IMF #{} unknown".format(self.imf))

    def process(self, sed):
        """Add the convolution of a Maraston 2005 SSP to the SED

        Parameters
        ----------
        sed: pcigale.sed.SED
            SED object.

        """
        sfh_time, sfh_sfr = sed.sfh
        ssp = self.ssp

        # Age of the galaxy at each time of the SFH
        sfh_age = np.max(sfh_time) - sfh_time

        # First, we process the young population (age lower than the
        # separation age.)
        young_sfh = np.copy(sfh_sfr)
        young_sfh[sfh_age > self.separation_age] = 0
        young_masses, young_spectrum = ssp.convolve(sfh_time, young_sfh)

        # Then, we process the old population. If the SFH is shorter than the
        # separation age then all the arrays will consist only of 0.
        old_sfh = np.copy(sfh_sfr)
        old_sfh[sfh_age <= self.separation_age] = 0
        old_masses, old_spectrum = ssp.convolve(sfh_time, old_sfh)

        sed.add_module(self.name, self.parameters)

        sed.add_info('stellar.imf', self.imf)
        sed.add_info('stellar.metallicity', self.metallicity)
        sed.add_info('stellar.old_young_separation_age', self.separation_age)

        sed.add_info('stellar.mass_total_old', old_masses[0], True)
        sed.add_info('stellar.mass_alive_old', old_masses[1], True)
        sed.add_info('stellar.mass_white_dwarf_old', old_masses[2], True)
        sed.add_info('stellar.mass_neutron_old', old_masses[3], True)
        sed.add_info('stellar.mass_black_hole_old', old_masses[4], True)

        sed.add_info('stellar.mass_total_young', young_masses[0], True)
        sed.add_info('stellar.mass_alive_young', young_masses[1], True)
        sed.add_info('stellar.mass_white_dwarf_young', young_masses[2], True)
        sed.add_info('stellar.mass_neutron_young', young_masses[3], True)
        sed.add_info('stellar.mass_black_hole_young', young_masses[4], True)

        sed.add_info('stellar.mass_total',
                     old_masses[0] + young_masses[0], True)
        sed.add_info('stellar.mass_alive',
                     old_masses[1] + young_masses[1], True)
        sed.add_info('stellar.mass_white_dwarf',
                     old_masses[2] + young_masses[2], True)
        sed.add_info('stellar.mass_neutron',
                     old_masses[3] + young_masses[3], True)
        sed.add_info('stellar.mass_black_hole',
                     old_masses[4] + young_masses[4], True)

        sed.add_contribution("stellar.old",
                             ssp.wavelength_grid,
                             old_spectrum)
        sed.add_contribution("stellar.young",
                             ssp.wavelength_grid,
                             young_spectrum)

# SedModule to be returned by get_module
Module = M2005
