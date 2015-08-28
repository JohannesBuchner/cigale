# -*- coding: utf-8 -*-
# Copyright (C) 2013 Department of Physics, University of Crete
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Laure Ciesla

"""
Dale et al. (2014) IR models module
===================================

This module implements the Dale (2014) infra-red models.

"""

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

    parameter_list = dict([
        ('fracAGN', (
            'float',
            "AGN fraction. It is not recommended to combine this AGN emission "
            "with the of Fritz et al. (2006) models.",
            0.0
        )),
        ('alpha', (
            'float',
            "Alpha slope. Possible values are: 0.0625, 0.1250, 0.1875, "
            "0.2500, 0.3125, 0.3750, 0.4375, 0.5000, 0.5625, 0.6250, 0.6875, "
            "0.7500, 0.8125, 0.8750, 0.9375, 1.0000, 1.0625, 1.1250, 1.1875, "
            "1.2500, 1.3125, 1.3750, 1.4375, 1.5000, 1.5625, 1.6250, 1.6875, "
            "1.7500, 1.8125, 1.8750, 1.9375, 2.0000, 2.0625, 2.1250, 2.1875, "
            "2.2500, 2.3125, 2.3750, 2.4375, 2.5000, 2.5625, 2.6250, 2.6875, "
            "2.7500, 2.8125, 2.8750, 2.9375, 3.0000, 3.0625, 3.1250, 3.1875, "
            "3.2500, 3.3125, 3.3750, 3.4375, 3.5000, 3.5625, 3.6250, 3.6875, "
            "3.7500, 3.8125, 3.8750, 3.9375, 4.0000",
            2.
        ))
    ])

    out_parameter_list = dict([
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
            self.model_quasar = database.get_dale2014(1.00, 0.0)

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

        frac_agn = self.parameters["fracAGN"]

        if frac_agn < 1.:
            L_AGN = luminosity * (1./(1.-frac_agn) - 1.)
        else:
            raise Exception("AGN fraction is exactly 1. Behaviour "
                            "undefined.")

        sed.add_module(self.name, self.parameters)
        sed.add_info("agn.fracAGN", self.parameters["fracAGN"])
        sed.add_info("dust.alpha", self.parameters["alpha"])

        sed.add_contribution('dust', self.model_sb.wave,
                             luminosity * self.model_sb.lumin)
        if frac_agn != 0.:
            sed.add_contribution('agn', self.model_quasar.wave,
                                 L_AGN * self.model_quasar.lumin)

# CreationModule to be returned by get_module
Module = Dale2014
