# -*- coding: utf-8 -*-
# Copyright (C) 2013, 2014 Department of Physics, University of Crete
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Laure Ciesla

"""
Fritz et al. (2006) AGN dust torus emission module
==================================================

This module implements the Fritz et al. (2006) models.

"""
from collections import OrderedDict

import numpy as np

from pcigale.data import Database
from . import CreationModule


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
            "cigale_list(options=10. & 30. & 60. & 100. & 150.)",
            "Ratio of the maximum to minimum radii of the dust torus. "
            "Possible values are: 10, 30, 60, 100, 150.",
            60.
        )),
        ('tau', (
            "cigale_list(options=0.1 & 0.3 & 0.6 & 1.0 & 2.0 & 3.0 & 6.0 & "
            "10.0)",
            "Optical depth at 9.7 microns. "
            "Possible values are: 0.1, 0.3, 0.6, 1.0, 2.0, 3.0, 6.0, 10.0.",
            1.0
        )),
        ('beta', (
            "cigale_list(options=-1.00 & -0.75 & -0.50 & -0.25 & 0.00)",
            "Beta. Possible values are: -1.00, -0.75, -0.50, -0.25, 0.00.",
            -0.50
        )),
        ('gamma', (
            'cigale_list(options=0.0 & 2.0 & 4.0 & 6.0)',
            "Gamma. Possible values are: 0.0, 2.0, 4.0, 6.0.",
            4.0
        )),
        ('opening_angle', (
            'cigale_list(options=60. & 100. & 140.)',
            "Full opening angle of the dust torus (Fig 1 of Fritz 2006). "
            "Possible values are: 60., 100., 140.",
            100.
        )),
        ('psy', (
            'cigale_list(options=0.001 & 10.100 & 20.100 & 30.100 & 40.100 & '
            '50.100 & 60.100 & 70.100 & 80.100 & 89.990)',
            "Angle between equatorial axis and line of sight. "
            "Psy = 90◦ for type 1 and Psy = 0° for type 2. Possible values "
            "are: 0.001, 10.100, 20.100, 30.100, 40.100, 50.100, 60.100, "
            "70.100, 80.100, 89.990.",
            50.100
        )),
        ('fracAGN', (
            'cigale_list(minvalue=0., maxvalue=1.)',
            "AGN fraction.",
            0.1
        ))
    ])

    def _init_code(self):
        """Get the template set out of the database"""
        self.r_ratio = float(self.parameters["r_ratio"])
        self.tau = float(self.parameters["tau"])
        self.beta = float(self.parameters["beta"])
        self.gamma = float(self.parameters["gamma"])
        self.opening_angle = (180. - self.parameters["opening_angle"]) / 2.
        self.psy = float(self.parameters["psy"])
        self.fracAGN = float(self.parameters["fracAGN"])

        with Database() as base:
            self.fritz2006 = base.get_fritz2006(self.r_ratio, self.tau,
                                                self.beta, self.gamma,
                                                self.opening_angle, self.psy)

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
        sed.add_info('agn.r_ratio', self.r_ratio)
        sed.add_info('agn.tau', self.tau)
        sed.add_info('agn.beta', self.beta)
        sed.add_info('agn.gamma', self.gamma)
        sed.add_info('agn.opening_angle', self.parameters["opening_angle"])
        sed.add_info('agn.psy', self.psy)
        sed.add_info('agn.fracAGN', self.fracAGN)

        # Compute the AGN luminosity
        if self.fracAGN < 1.:
            agn_power = luminosity * (1./(1.-self.fracAGN) - 1.)
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
