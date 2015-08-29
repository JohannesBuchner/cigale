# -*- coding: utf-8 -*-
# Copyright (C) 2015 Centre de données Astrophysiques de Marseille
# Copyright (C) 2015 Institute of Astronomy, University of Cambridge
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Denis Burgarella

"""
Modified Black Body module
======================================

This module implements a modified black body (MBB). This MMB can be 
a stand-alone IR emission or it can come as an additional component. 
The energy balance can be set to check the presence of a component 
or not to account for the emission of a region completely embedded 
in dust and not visible in the wavelength range.

"""

from collections import OrderedDict
import numpy as np
import scipy.constants as cst
from . import CreationModule


class MBB(CreationModule):
    """One modified black body IR re-emission

    Given an amount of attenuation (e.g. resulting from the action of a dust
    attenuation module) this module normalises MBB plus any previous IR
    contribution to this amount of energy. The final SED allows to keep the
    energy balance or not..

    """

    parameter_list = OrderedDict([
        ("epsilon_mbb", (
            "float",
            "Fraction [>= Ø] of L_dust(energy balance) in the MBB",
            0.5
        )),
        ("t_mbb", (
            "float",
            "Temperature of black body in K.",
            50.
        )),
        ("beta_mbb", (
            "float",
            "Emissivity index of modified black body.",
            1.5
        )),
        ("energy_balance", (
            "boolean",
            "Energy balance checked?"
            "If False, Lum[MBB] not taken into account in energy balance",
            False
        )),
    ])

    def _init_code(self):
        """Build the model for a given set of parameters."""

        epsilon = float(self.parameters["epsilon_mbb"])
        T = float(self.parameters["t_mbb"])
        beta = float(self.parameters["beta_mbb"])
        
        if epsilon < 0:
           epsilon = 0.
           print("epsilon_mbb must >= 0, we set epsilon_mbb = 0.0")
        
        # We define various constants necessary to compute the model. For
        # consistency, we define the speed of light in nm s¯¹ rather than in
        # m s¯¹.
        c = cst.c * 1e9
        lambda_0 = 200e3

        self.wave = np.logspace(3., 6., 1000.)
        conv = c / (self.wave * self.wave)

        self.lumin_mbb = (conv * (1. - np.exp(-(lambda_0 / self.wave)
                                ** beta)) * (c / self.wave) ** 3. / (np.exp(
                                cst.h * c / (self.wave * cst.k * T)) - 1.))

        # TODO, save the right normalisation factor to retrieve the dust mass
        norm = np.trapz(self.lumin_mbb, x=self.wave)
        self.lumin_mbb /= norm

    def process(self, sed):
        """Add the IR re-emission contributions.

        Parameters
        ----------
        sed: pcigale.sed.SED object

        """
        if 'dust.luminosity' not in sed.info:
            sed.add_info('dust.luminosity', 1., True)
        luminosity = sed.info['dust.luminosity']

        sed.add_module(self.name, self.parameters)
        sed.add_info("dust.t_mbb", self.parameters["t_mbb"])
        sed.add_info("dust.beta_mbb", self.parameters["beta_mbb"])
        sed.add_info("dust.epsilon_mbb", self.parameters["epsilon_mbb"])
        epsilon = float(self.parameters["epsilon_mbb"])
        energy_balance = (self.parameters["energy_balance"].lower() == "true")

        if energy_balance:
           # Since we can have another contribution to L_dust and the modified black body
           # enters into the energy budget, energy balance, we have to save a new negative 
           # component for each one, previously existing as:
           # luminosity_other_balance = -luminosity_other * (1-epsilon); 
           # epsilon being the contribution of the present MBB to L_dust.             
           other_dust_contributions = [contrib for contrib in sed.contribution_names 
                                       if "dust" in contrib]
           for item in other_dust_contributions:
              item_balance = item + '_balance'
              lumin = sed.get_lumin_contribution(item)
              wavelength = sed.wavelength_grid
              sed.add_info(item_balance, 1., True)
              sed.add_contribution(item_balance, wavelength, -lumin * epsilon)

        # If the modified black body does not enter into the energy budget, 
        # we do not change the luminosity of other dust contributions.
        # The luminosity of the modified black body L_MBB = epsilon * L_dust. 
        # The total dust luminosity will be : L_dust(inside energy budget) + L_MBB.
         
        # We add the contribution of the MBB to L_dust.
        sed.add_contribution('dust.mbb', self.wave,
                             luminosity * epsilon * self.lumin_mbb)
#
# CreationModule to be returned by get_module
Module = MBB
