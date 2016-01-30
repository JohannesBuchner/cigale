# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013 Institute of Astronomy, University of Cambridge
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Médéric Boquien

"""
Draine and Li (2007) IR models module
=====================================

This module implements the Draine and Li (2007) infra-red models.

"""

from collections import OrderedDict

import numpy as np

from pcigale.data import Database
from . import CreationModule


class DL2007(CreationModule):
    """Draine and Li (2007) templates IR re-emission module

    Given an amount of attenuation (e.g. resulting from the action of a dust
    attenuation module) this module normalises the Draine and Li (2007)
    template corresponding to a given α to this amount of energy and add it
    to the SED.

    Information added to the SED: NAME_alpha.

    """

    parameter_list = OrderedDict([
        ('qpah', (
            'float',
            "Mass fraction of PAH. Possible values are: 0.47, 1.12, 1.77, "
            "2.50, 3.19, 3.90, 4.58.",
            2.50
        )),
        ('umin', (
            'float',
            "Minimum radiation field. Possible values are: 0.10, 0.15, 0.20, "
            "0.30, 0.40, 0.50, 0.70, 0.80, 1.00, 1.20, 1.50, 2.00, 2.50, "
            "3.00, 4.00, 5.00, 7.00, 8.00, 10.0, 12.0, 15.0, 20.0, 25.0.",
            1.0
        )),
        ('umax', (
            'float',
            "Maximum radiation field. Possible values are: 1e3, 1e4, 1e5, "
            "1e6.",
            1e6
        )),
        ('gamma', (
            'float',
            "Fraction illuminated from Umin to Umax. Possible values between "
            "0 and 1.",
            0.1
        ))
    ])

    def _init_code(self):
        """Get the model out of the database"""

        self.qpah = float(self.parameters["qpah"])
        self.umin = float(self.parameters["umin"])
        self.umax = float(self.parameters["umax"])
        self.gamma = float(self.parameters["gamma"])

        with Database() as database:
            self.model_minmin = database.get_dl2007(self.qpah, self.umin,
                                                    self.umin)
            self.model_minmax = database.get_dl2007(self.qpah, self.umin,
                                                    self.umax)

        # The models in memory are in W/nm for 1 kg of dust. At the same time
        # we need to normalize them to 1 W here to easily scale them from the
        # power absorbed in the UV-optical. If we want to retrieve the dust
        # mass at a later point, we have to save their "emissivity" per unit
        # mass in W (kg of dust)¯¹, The gamma parameter does not affect the
        # fact that it is for 1 kg because it represents a mass fraction of
        # each component.
        self.emissivity = np.trapz((1.-self.gamma) * self.model_minmin.lumin +
                                   self.gamma * self.model_minmax.lumin,
                                   x=self.model_minmin.wave)

        # We want to be able to display the respective contributions of both
        # components, therefore we keep they separately.
        self.model_minmin.lumin *= (1. - self.gamma) / self.emissivity
        self.model_minmax.lumin *= self.gamma / self.emissivity

    def process(self, sed):
        """Add the IR re-emission contributions

        Parameters
        ----------
        sed: pcigale.sed.SED object
        parameters: dictionary containing the parameters

        """
        if 'dust.luminosity' not in sed.info:
            sed.add_info('dust.luminosity', 1., True)
        luminosity = sed.info['dust.luminosity']

        sed.add_module(self.name, self.parameters)
        sed.add_info('dust.qpah', self.qpah)
        sed.add_info('dust.umin', self.umin)
        sed.add_info('dust.umax', self.umax)
        sed.add_info('dust.gamma', self.gamma)
        # To compute the dust mass we simply divide the luminosity in W by the
        # emissivity in W/kg of dust.
        sed.add_info('dust.mass', luminosity / self.emissivity, True)

        sed.add_contribution('dust.Umin_Umin', self.model_minmin.wave,
                             luminosity * self.model_minmin.lumin)
        sed.add_contribution('dust.Umin_Umax', self.model_minmax.wave,
                             luminosity * self.model_minmax.lumin)

# CreationModule to be returned by get_module
Module = DL2007
