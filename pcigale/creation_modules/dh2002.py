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
            "Alpha slope. Possible values between 0.0625 and 4.000.",
            2.
        ))
    ])

    out_parameter_list = OrderedDict([('alpha', 'Alpha slope.')])

    def _init_code(self):
        """Get the template set out of the database"""
        with Database() as database:
            self.dh2002 = database.get_dh2002()

    def process(self, sed):
        """Add the IR re-emission contributions

        Parameters
        ----------
        sed: pcigale.sed.SED object

        """
        alpha = float(self.parameters["alpha"])

        if 'dust.luminosity' not in sed.info.keys():
            sed.add_info('dust.luminosity', 1., True)
        luminosity = sed.info['dust.luminosity']

        ir_template = self.dh2002.get_template(alpha)

        sed.add_module(self.name, self.parameters)
        sed.add_info("dust.alpha", alpha)

        sed.add_contribution('dust', self.dh2002.wavelength_grid,
                             luminosity * ir_template)

# CreationModule to be returned by get_module
Module = DH2002
