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
            "Mass fraction of PAH",
            None
        )),
        ('umin', (
            'float',
            "Minimum radiation field",
            None
        )),
        ('umax', (
            'float',
            "Maximum radiation field",
            None
        )),
        ('gamma', (
            'float',
            "Fraction illuminated from Umin to Umax",
            None
        ))
    ])

    out_parameter_list = OrderedDict([
        ('qpah', 'Mass fraction of PAH'),
        ('umin', 'Minimum radiation field'),
        ('umax', 'Maximum radiation field'),
        ('gamma', 'Fraction illuminated from Umin to Umax')
    ])

    def _init_code(self):
        """Get the model out of the database"""

        qpah = self.parameters["qpah"]
        umin = self.parameters["umin"]
        umax = self.parameters["umax"]
        gamma = self.parameters["gamma"]

        with Database() as database:
            self.model_minmin = database.get_dl2007(qpah, umin, umin)
            self.model_minmax = database.get_dl2007(qpah, umin, umax)

        # The models in memory are in W/nm for 1 kg of dust. At the same time
        # we need to normalize them to 1 W here to easily scale them from the
        # power absorbed in the UV-optical. If we want to retrieve the dust
        # mass at a later point, we have to save their "emissivity" per unit
        # mass in W kg¯¹, The gamma parameter does not affect the fact that it
        # is for 1 kg because it represents a mass fraction of each component.
        self.emissivity = np.trapz((1. - gamma) * self.model_minmin.lumin +
                                   gamma * self.model_minmax.lumin,
                                   x=self.model_minmin.wave)

        # We want to be able to display the respective contributions of both
        # components, therefore we keep they separately.
        self.model_minmin.lumin *= (1. - gamma) / self.emissivity
        self.model_minmax.lumin *= gamma / self.emissivity

    def process(self, sed):
        """Add the IR re-emission contributions

        Parameters
        ----------
        sed  : pcigale.sed.SED object
        parameters : dictionary containing the parameters

        """
        if 'dust.luminosity' not in sed.info.keys():
            sed.add_info('dust.luminosity', 1., True)
        luminosity = 1.

        sed.add_module(self.name, self.parameters)
        sed.add_info('dust.qpah', self.parameters["qpah"])
        sed.add_info('dust.umin', self.parameters["umin"])
        sed.add_info('dust.umax', self.parameters["umax"])
        sed.add_info('dust.gamma', self.parameters["gamma"])

        sed.add_contribution('dust.Umin_Umin', self.model_minmin.wave,
                             luminosity * self.model_minmin.lumin)
        sed.add_contribution('dust.Umin_Umax', self.model_minmax.wave,
                             luminosity * self.model_minmax.lumin)

# CreationModule to be returned by get_module
Module = DL2007
