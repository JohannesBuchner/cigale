# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013 Institute of Astronomy, University of Cambridge
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Laure Ciesla

"""
Schreiber (2016) IR models module
=====================================

This module implements the Schreiber et al. (2016) infra-red models.

"""

from collections import OrderedDict
import numpy as np
from pcigale.data import Database
from . import CreationModule


class Schreiber2016(CreationModule):
    """Schreiber et al. (2016) templates IR re-emission module

    Given an amount of attenuation (e.g. resulting from the action of a dust
    attenuation module) this module normalises the Schreiber et al. (2016)
    template corresponding to a given α to this amount of energy and add it
    to the SED.

    """

    parameter_list = OrderedDict([
        ('tdust', (
            'float',
            "Dust temperature. "
            "Between 15 and 60K, with 1K step.",
            20.
        )),
        ('fpah', (
            'float',
            "Mass fraction of PAH. "
            "Between 0 and 1.",
            0.05
        ))
    ])

    out_parameter_list = OrderedDict([
        ('tdust', 'Dust temperature'),
        ('fpah', 'Mass fraction of PAH')
    ])

    def _init_code(self):
        """Get the model out of the database"""

        self.tdust = float(self.parameters["tdust"])
        self.fpah = float(self.parameters["fpah"])
        with Database() as database:
            self.model_dust = database.get_schreiber2016(0, self.tdust)
            self.model_pah = database.get_schreiber2016(1, self.tdust)

        # The models in memory are in W/nm/kg. At the same time we
        # need to normalize them to 1 W here to easily scale them from the
        # power absorbed in the UV-optical. If we want to retrieve the dust
        # mass at a later point, we have to save their "emissivity" per unit
        # mass in W kg¯¹, The gamma parameter does not affect the fact that it
        # is for 1 kg because it represents a mass fraction of each component.

        self.emissivity = np.trapz((1. - self.fpah) * self.model_dust.lumin +
                                   self.fpah * self.model_pah.lumin,
                                   x=self.model_dust.wave)

        # We want to be able to display the respective contributions of both
        # components, therefore we keep they separately.
        self.model_dust.lumin *= (1. - self.fpah) / self.emissivity
        self.model_pah.lumin *= self.fpah / self.emissivity

    def process(self, sed):
        """Add the IR re-emission contributions

        Parameters
        ----------
        sed: pcigale.sed.SED object
        parameters: dictionary containing the parameters

        """
        if 'dust.luminosity' not in sed.info.keys():
            sed.add_info('dust.luminosity', 1., True)
        luminosity = sed.info['dust.luminosity']

        sed.add_module(self.name, self.parameters)
        sed.add_info('dust.tdust', self.tdust)
        sed.add_info('dust.fpah', self.fpah)
        # To compute the dust mass we simply divide the luminosity by the
        # emissivity and then by the expected MH/Mdust as the emissivity was
        # computed for 1 kg of H. Note that we take 100 here but it should vary
        # with the exact model. Fix that later. Maybe directly in the database.
        sed.add_info('dust.mass', luminosity / self.emissivity, True)

        sed.add_contribution('dust.dust_continuum', self.model_dust.wave,
                             luminosity * self.model_dust.lumin)
        sed.add_contribution('dust.pah', self.model_pah.wave,
                             luminosity * self.model_pah.lumin)

# CreationModule to be returned by get_module
Module = Schreiber2016
