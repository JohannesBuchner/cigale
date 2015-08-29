# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013 Institute of Astronomy, University of Cambridge
# Copyright (C) 2014 University of Crete
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Laure Ciesla & Médéric Boquien

"""
Radio module
=============================

This module implements the radio emission of galaxies, taking into account only
 the non-thermal emission. The thermal emission is handled by the nebular
 module. The parameters that this module takes as input are:
- the value of the coefficient of the FIR/radio correlation
- the value of the spectral index of the power law emission from synchrotron.

"""

from collections import OrderedDict
import numpy as np
import scipy.constants as cst
from . import CreationModule


class Radio(CreationModule):
    """Radio emission

    Given the number of Lyman photons, the module computes the free-free
    (thermal) emission of galaxies. Based on the SN collapse rate, the module
    computes the synchrotron (non-thermal) emission of galaxies.

    """

    parameter_list = OrderedDict([
        ("qir", (
            "float",
            "The value of the FIR/radio correlation coefficient.",
            2.58
        )),
        ("alpha", (
            "float",
            "The slope of the power-law synchrotron emission.",
            0.8
        ))
    ])

    out_parameter_list = dict([
        ("qir", "The value of the FIR/radio correlation coefficient."),
        ("alpha", "The slope of the power-law synchrotron emission.")
    ])

    def _init_code(self):
        """Build the model for a given set of parameters."""

        qir = float(self.parameters["qir"])
        alpha = float(self.parameters["alpha"])

        # We define various constants necessary to compute the model. For
        # consistency, we define the speed of light in nm s¯¹ rather than in
        # m s¯¹.
        c = cst.c * 1e9
        # We define the wavelength range for the non thermal emission
        self.wave = np.logspace(5., 9., 1000.)
        # We compute the synchrotron emission normalised at 21cm
        self.lumin_nonthermal = ((1./self.wave)**(-alpha + 2) /
                                 (1./2.1e8)**(-alpha + 2))
        # Normalisation factor from the FIR/radio correlation to apply to the
        # IR luminosity
        S21cm = (1. / (10**qir*3.75e12)) * (c/(2.1e8)**2)
        self.lumin_nonthermal *= S21cm

    def process(self, sed):
        """Add the radio contribution.

        Parameters
        ----------
        sed: pcigale.sed.SED object

        """
        if 'dust.luminosity' not in sed.info:
            sed.add_info('dust.luminosity', 1., True)
        luminosity = sed.info['dust.luminosity']

        sed.add_module(self.name, self.parameters)
        sed.add_info("radio.qir", self.parameters["qir"])
        sed.add_info("radio.alpha", self.parameters["alpha"])
        sed.add_contribution('radio_nonthermal', self.wave,
                             self.lumin_nonthermal * luminosity)

# CreationModule to be returned by get_module
Module = Radio
