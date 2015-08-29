# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013 Institute of Astronomy, University of Cambridge
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Médéric Boquien

"""
Updated Draine and Li (2007) IR models module
=====================================

This module implements the updated Draine and Li (2007) infrared models.

"""

from collections import OrderedDict
import numpy as np
from pcigale.data import Database
from . import CreationModule


class DL2014(CreationModule):
    """Updated Draine and Li (2007) templates IR re-emission module

    Given an amount of attenuation (e.g. resulting from the action of a dust
    attenuation module) this module normalises the updated Draine and Li (2007)
    model corresponding to a given set of parameters to this amount of energy
    and add it to the SED.

    Information added to the SED: NAME_alpha.

    """

    parameter_list = OrderedDict([
        ('qpah', (
            'float',
            "Mass fraction of PAH. Possible values are: 0.47, 1.12, 1.77, "
            "2.50, 3.19, 3.90, 4.58, 5.26, 5.95, 6.63, 7.32.",
            None
        )),
        ('umin', (
            'float',
            "Minimum radiation field. Possible values are: 0.100, 0.120, "
            "0.150, 0.170, 0.200, 0.250, 0.300, 0.350, 0.400, 0.500, 0.600, "
            "0.700, 0.800, 1.000, 1.200, 1.500, 1.700, 2.000, 2.500, 3.000, "
            "3.500, 4.000, 5.000, 6.000, 7.000, 8.000, 10.00, 12.00, 15.00, "
            "17.00, 20.00, 25.00, 30.00, 35.00, 40.00, 50.00.",
            None
        )),
        ('alpha', (
            'float',
            "Powerlaw slope dU/dM propto U^alpha. Possible values are: 1.0, "
            "1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, "
            "2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0.",
            None
        )),
        ('gamma', (
            'float',
            "Fraction illuminated from Umin to Umax. Possible values between "
            "0 and 1.",
            None
        ))
    ])

    out_parameter_list = dict([
        ('qpah', 'Mass fraction of PAH'),
        ('umin', 'Minimum radiation field'),
        ('alpha', 'Power law slope dU/dM∝U¯ᵅ'),
        ('gamma', 'Fraction illuminated from Umin to Umax')
    ])

    def _init_code(self):
        """Get the model out of the database"""

        self.qpah = self.parameters["qpah"]
        self.umin = self.parameters["umin"]
        self.alpha = self.parameters["alpha"]
        self.gamma = self.parameters["gamma"]
        self.umax = 1e7

        with Database() as database:
            self.model_minmin = database.get_dl2014(self.qpah, self.umin,
                                                    self.umin, 1.)
            self.model_minmax = database.get_dl2014(self.qpah, self.umin,
                                                    self.umax, self.alpha)

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
        sed.add_info('dust.alpha', self.alpha)
        sed.add_info('dust.gamma', self.gamma)
        # To compute the dust mass we simply divide the luminosity in W by the
        # emissivity in W/kg of dust.
        sed.add_info('dust.mass', luminosity / self.emissivity, True)

        sed.add_contribution('dust.Umin_Umin', self.model_minmin.wave,
                             luminosity * self.model_minmin.lumin)
        sed.add_contribution('dust.Umin_Umax', self.model_minmax.wave,
                             luminosity * self.model_minmax.lumin)

# CreationModule to be returned by get_module
Module = DL2014
