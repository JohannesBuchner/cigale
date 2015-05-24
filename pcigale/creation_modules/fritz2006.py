# -*- coding: utf-8 -*-
# Copyright (C) 2013, 2014 Department of Physics, University of Crete
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Laure Ciesla

"""
Fritz et al. (2006) AGN dust torus emission module
==================================================

This module implements the Fritz et al. (2006) models.

"""
import numpy as np
from collections import OrderedDict
from pcigale.data import Database
from . import CreationModule
from pcigale.sed.cosmology import cosmology


class Fritz2006(CreationModule):
    """Fritz et al. (2006) AGN dust torus emission

    The AGN emission is computed from the library of Fritz et al. (2006) from
    which all of the models are available. They take into account two emission
    components linked to the AGN. The first one is the isotropic emission of
    the central source, which is assumed to be point-like. This emission is a
    composition of power laws with variable indices, in the wavelength range of
    0.001-20 microns. The second one is the thermal and scattering dust torus
    emission. The conservation of the energy is always verified within 1% for
    typical solutions, and up to 10% in the case of very high optical depth and
    non-constant dust density. We refer the reader to Fritz et al. (2006) for
    more information on the library.

    The relative normalization of these components is handled through a
    parameter which is the fraction of the total IR luminosity due to the AGN
    so that: L_AGN = fracAGN * L_IRTOT, where L_AGN is the AGN luminosity,
    fracAGN is the contribution of the AGN to the total IR luminosity
    (L_IRTOT), i.e. L_Starburst+L_AGN.

    """

    parameter_list = OrderedDict([
        ('r_ratio', (
            'float',
            "Ratio of the maximum to minimum radii of the dust torus. "
            "Possible values are: 10, 30, 60, 100, 150.",
            60.
        )),
        ('tau', (
            'float',
            "Optical depth at 9.7 microns. "
            "Possible values are: 0.1, 0.3, 0.6, 1.0, 2.0, 3.0, 6.0, 10.0.",
            1.0
        )),
        ('beta', (
            'float',
            "Beta. Possible values are:-1.00, -0.75, -0.50, -0.25, 0.00.",
            -0.50
        )),
        ('gamma', (
            'float',
            "Gamma. Possible values are: 0.0, 2.0, 4.0, 6.0.",
            4.0
        )),
        ('opening_angle', (
            'float',
            "Full opening angle of the dust torus (Fig 1 of Fritz 2006). "
            "Possible values are: 60., 100., 140.",
            100.
        )),
        ('psy', (
            'float',
            "Angle between equatorial axis and line of sight. "
            "Psy = 90◦ for type 1 and Psy = 0° for type 2. Possible values "
            "are: 0.001, 10.100, 20.100, 30.100, 40.100, 50.100, 60.100, "
            "70.100, 80.100, 89.990.",
            50.100
        )),
        ('fracAGN', (
            'float',
            "AGN fraction.",
            0.1
        ))
    ])

    out_parameter_list = OrderedDict([
        ('fracAGN', 'Contribution of the AGN'),
        ('agn.therm_luminosity', 'Luminosity of the AGN contribution due to '
                                 'the dust torus'),
        ('agn.scatt_luminosity', 'Luminosity of the AGN contribution due to '
                                 'the photon scattering'),
        ('agn.agn_luminosity', 'Luminosity of the AGN contribution due to the '
                               'central source'),
        ('agn.luminosity', 'Total luminosity of the AGN contribution')
    ])

    def _init_code(self):
        """Get the template set out of the database"""
        r_ratio = self.parameters["r_ratio"]
        tau = self.parameters["tau"]
        beta = self.parameters["beta"]
        gamma = self.parameters["gamma"]
        opening_angle = (180. - self.parameters["opening_angle"]) / 2.
        psy = self.parameters["psy"]

        with Database() as base:
            self.fritz2006 = base.get_fritz2006(r_ratio, tau, beta, gamma,
                                                opening_angle, psy)

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
        if fracAGN < 1.:
            agn_power = luminosity * (1./(1.-fracAGN) - 1.)
            l_agn_therm = agn_power
            l_agn_scatt = np.trapz(agn_power * self.fritz2006.lumin_scatt,
                                   x=self.fritz2006.wave)
            l_agn_agn = np.trapz(agn_power * self.fritz2006.lumin_agn,
                                 x=self.fritz2006.wave)
            l_agn_total = l_agn_therm + l_agn_scatt + l_agn_agn

        else:
            raise Exception("AGN fraction is exactly 1. Behaviour "
                            "undefined.")

        sed.add_info('agn.therm_luminosity', l_agn_therm, True)
        sed.add_info('agn.scatt_luminosity', l_agn_scatt, True)
        sed.add_info('agn.agn_luminosity', l_agn_agn, True)
        sed.add_info('agn.luminosity', l_agn_total, True)

        sed.add_contribution('agn.fritz2006_therm', self.fritz2006.wave,
                             agn_power * self.fritz2006.lumin_therm)
        sed.add_contribution('agn.fritz2006_scatt', self.fritz2006.wave,
                             agn_power * self.fritz2006.lumin_scatt)
        sed.add_contribution('agn.fritz2006_agn', self.fritz2006.wave,
                             agn_power * self.fritz2006.lumin_agn)

# CreationModule to be returned by get_module
Module = Fritz2006
