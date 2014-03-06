# -*- coding: utf-8 -*-
# Copyright (C) 2013 Department of Physics, University of Crete
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Laure Ciesla

"""
Dale et al. (2014) IR models module
===================================

This module implements the Dale (2014) infra-red models.

"""

from collections import OrderedDict
from pcigale.data import Database
from . import CreationModule


class Dale2014(CreationModule):
    """Dale et al. (2014) templates IR re-emission

    Given an amount of attenuation (e.g. resulting from the action of a dust
    attenuation module) this module normalises the Dale et al (2014)
    template corresponding to a given α to this amount of energy and add it
    to the SED.

    Information added to the SED: NAME_fracAGN, NAME_alpha.

    """

    parameter_list = OrderedDict([
        ('fracAGN', (
            'float',
            "Contribution of the AGN",
            None
        )),
        ('alpha', (
            'float',
            "Alpha slope",
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
        ('fracAGN', 'Contribution of the AGN'),
        ('alpha', 'Alpha slope'),
        ('lir', 'Total IR luminosity between 8 and 1000 microns (AGN + SB)')
    ])

    def _init_code(self):
        """
        Get the models out of the database
        model_sb corresponds to an AGN fraction of 0%: only the dust heated by
        star formation
        model_quasar corresponds to an AGN fraction of 100%: only the SED of
        the quasar
        The energy attenuated is re-injected in model_sb only.
        """
        alpha = self.parameters["alpha"]

        with Database() as database:
            self.model_sb = database.get_dale2014(0.00, alpha)
            self.model_quasar = database.get_dale2014(1.00, alpha)

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

        frac_agn = self.parameters["fracAGN"]

        sed.add_module(self.name, self.parameters)
        sed.add_info("agn.fracAGN", self.parameters["fracAGN"])
        sed.add_info("dust.alpha", self.parameters["alpha"])

        sed.add_contribution('dust', self.model_sb.wave,
                             luminosity * self.model_sb.lumin)
        sed.add_contribution('agn', self.model_quasar.wave,
                             frac_agn * luminosity * self.model_quasar.lumin)

# CreationModule to be returned by get_module
Module = Dale2014
