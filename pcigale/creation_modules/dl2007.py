# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Médéric Boquien <mederic.boquien@oamp.fr>

from collections import OrderedDict
import numpy as np
from pcigale.data import Database
from . import SEDCreationModule


class Module(SEDCreationModule):
    """
    Module computing the infra-red re-emission corresponding to an amount of
    attenuation using the Draine and Li (2007) models.

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
        )),
        ('attenuation_value_keys', (
            'string',
            "Keys of the SED information dictionary where the module will "
            "look for the attenuation (in W) to re-emit. You can give several "
            "keys separated with a & (don't use commas), a re-emission "
            "contribution will be added for each key.",
            "attenuation"
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

        database = Database()
        self.model_minmin = database.get_dl2007(qpah, umin, umin)
        self.model_minmax = database.get_dl2007(qpah, umin, umax)
        database.session.close_all()

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
        attenuation_value_keys = [
            item.strip() for item in
            self.parameters["attenuation_value_keys"].split("&")]

        # Base name for adding information to the SED.
        name = self.name or 'dl2007'

        sed.add_module(name, self.parameters)
        sed.add_info('qpah' + self.postfix, self.parameters["qpah"])
        sed.add_info('umin' + self.postfix, self.parameters["umin"])
        sed.add_info('umax' + self.postfix, self.parameters["umax"])
        sed.add_info('gamma' + self.postfix, self.parameters["gamma"])

        for attenuation in attenuation_value_keys:
            sed.add_contribution(
                name + '_Umin_Umin_' + attenuation,
                self.model_minmin.wave,
                sed.info[attenuation] * self.model_minmin.lumin
            )
            sed.add_contribution(
                name + '_Umin_Umax_' + attenuation,
                self.model_minmax.wave,
                sed.info[attenuation] * self.model_minmax.lumin
            )
