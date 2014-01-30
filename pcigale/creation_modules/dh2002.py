# -*- coding: utf-8 -*-
# Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

"""
Dale and Helou (2002) IR models module
======================================

This module implements the Dale and Helou (2002) infra-red models.

"""

from collections import OrderedDict
from . import CreationModule
from ..data import Database


class DH2002(CreationModule):
    """Dale and Helou (2002) templates IR re-emission module

    Given an amount of attenuation (e.g. resulting from the action of a dust
    attenuation module) this module normalises the Dale and Helou (2002)
    template corresponding to a given α to this amount of energy and add it
    to the SED.

    Information added to the SED: NAME_alpha.

    """

    parameter_list = OrderedDict([
        ('alpha', (
            'float',
            "Alpha slope.",
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

    out_parameter_list = OrderedDict([('alpha', 'Alpha slope.')])

    def _init_code(self):
        """Get the template set out of the database"""
        with Database() as database:
            self.dh2002 = database.get_dh2002_infrared_templates()

    def process(self, sed):
        """Add the IR re-emission contributions

        Parameters
        ----------
        sed  : pcigale.sed.SED object

        """
        alpha = float(self.parameters["alpha"])
        attenuation_value_keys = [
            item.strip() for item in
            self.parameters["attenuation_value_keys"].split("&")]

        ir_template = self.dh2002.get_template(alpha)

        sed.add_module(self.name, self.parameters)
        sed.add_info("alpha" + self.postfix, alpha)

        for attenuation in attenuation_value_keys:
            sed.add_contribution(
                self.name + "_" + attenuation,
                self.dh2002.wavelength_grid,
                sed.info[attenuation] * ir_template
            )

# CreationModule to be returned by get_module
Module = DH2002
