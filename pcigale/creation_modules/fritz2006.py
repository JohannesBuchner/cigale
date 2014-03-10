# -*- coding: utf-8 -*-
# Copyright (C) 2013, 2014 Department of Physics, University of Crete
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Laure Ciesla

"""
Fritz et al. (2006) AGN dust torus emission module
==================================================

TODO: Describe the module

"""

from collections import OrderedDict
from pcigale.data import Database
from . import CreationModule


class Fritz2006(CreationModule):
    """Fritz et al. (2006) AGN dust torus emission

    TODO: Describe de module actions.

    Information added to the SED: fracAGN, L_AGN.

    """

    parameter_list = OrderedDict([
        ('r_ratio', (
            'float',
            "Ratio of the maximum and minimum radii of the dust torus. "
            "Possible values are: 10, 30, 60, 100, 150",
            None
        )),
        ('tau', (
            'float',
            "Tau at 9.7microns. "
            "Possible values are: 0.1, 0.3, 0.6, 1.0, 2.0, 3.0, 6.0, 10.0",
            None
        )),
        ('beta', (
            'float',
            "Beta. Possible values are:-1.00, -0.75, -0.50, -0.25, 0.00",
            None
        )),
        ('gamma', (
            'float',
            "Gamma. Possible values are: 0.0, 2.0, 4.0, 6.0",
            None
        )),
        ('opening_angle', (
            'float',
            "Opening angle of the dust torus. Possible values are: 20, 40, 60",
            None
        )),
        ('psy', (
            'float',
            "Angle between AGN axis and line of sight. Possible values are: "
            "0.001, 10.100, 20.100, 30.100, 40.100, 50.100, 60.100, 70.100,"
            "80.100, 89.990",
            None
        )),
        ('fracAGN', (
            'float',
            "Contribution of the AGN"
            "",
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
        ('L_AGN', 'Luminosity of the AGN contribution')
    ])


    def _init_code(self):
        """Get the template set out of the database"""
        r_ratio = self.parameters["r_ratio"]
        tau = self.parameters["tau"]
        beta = self.parameters["beta"]
        gamma = self.parameters["gamma"]
        opening_angle = self.parameters["opening_angle"]
        psy = self.parameters["psy"]

        with Database() as base:
            self.fritz2006 = base.get_fritz2006(r_ratio, tau, beta, gamma,
                                                    opening_angle, psy)



    def process(self, sed):
        """Add the IR re-emission contributions

        Parameters
        ----------
        sed  : pcigale.sed.SED object
        parameters : dictionary containing the parameters

        """

        if 'dust.luminosity' not in sed.info.keys():
            sed.add_info('dust.luminosity',1.,True)
        luminosity = sed.info['dust.luminosity']


        fracAGN = self.parameters["fracAGN"]

        sed.add_module(self.name, self.parameters)
        sed.add_info('r_ratio', self.parameters["r_ratio"])
        sed.add_info('tau', self.parameters["tau"])
        sed.add_info('beta', self.parameters["beta"])
        sed.add_info('gamma', self.parameters["gamma"])
        sed.add_info('opening_angle', self.parameters["opening_angle"])
        sed.add_info('psy', self.parameters["psy"])
        sed.add_info('fracAGN', self.parameters["fracAGN"])


        # Compute the AGN luminosity
        L_AGN = fracAGN * (luminosity + 1)
        #sed.add_info("L_AGN" + self.postfix, self.parameters["L_AGN"])
        sed.add_contribution(
            'agn_fritz2006_therm',
                self.fritz2006.wave,
                L_AGN * self.fritz2006.lumin_therm
            )

        sed.add_contribution(
            'agn_fritz2006_scatt',
                self.fritz2006.wave,
                L_AGN * self.fritz2006.lumin_scatt
            )

        sed.add_contribution(
            'agn_fritz2006_agn',
                self.fritz2006.wave,
                L_AGN * self.fritz2006.lumin_agn
            )

# CreationModule to be returned by get_module
Module = Fritz2006
