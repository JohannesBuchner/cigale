# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Manuela Raimbault

import numpy as np
from collections import OrderedDict
from scipy.constants import c
from . import common


class Module(common.SEDCreationModule):

    parameter_list = OrderedDict([
        ("Nly_key", (
            "string",
            "Key in the SED info dictionary where the number of ionising "
            "photons is found.",
            "bc03_n_ly_young"
        )),
        ("metallicity_key", (
            "string",
            "Key in the SED info dictionary where the metallicity is found.",
            "bc03_metallicity"
        )),
        ("set_name", (
            "string",
            "Name of the lines component.",
            None
        ))
    ])

    out_parameter_list = OrderedDict()

    def process(self, sed):

        n_ly = sed.info[self.parameters["Nly_key"]]
        metallicity = sed.info[self.parameters["metallicity_key"]]

        f_esc = 0.  # fraction of n_ly which escapes from the galaxy
        f = 1.  # n_ly fraction which contributes to the ionisation, (1-f) goes
                # in the dust

        flow_H = 4.757e-13 * n_ly * (1-f_esc) * f  # erg.s^-1
        flow_H *= 1e-7  # W.s^-1

        wav = np.array([1215.67, 1335.00, 1663.00, 1909.00, 2141.00, 2326.00,
                        2798.00, 3727.00, 3869.00, 3889.00, 3970.00, 4026.00,
                        4068.60, 4076.35, 4101.73, 4363.00, 4471.00, 4711.00,
                        4861.32, 4958.91, 5006.84, 5199.00, 5755.00, 5876.00,
                        6300.00, 6312.00, 6548.05, 6562.8, 6583.45, 6678.00,
                        6716.00, 6730.00, 7065.00, 7135.79, 7319.99, 7330.73,
                        7751.11, 9068.60, 9530.85, 10286.73, 10320.49,
                        10336.41])

        # ratio for metallicity = 0.0004 et 0.0001 for now
        if (metallicity == 0.0004 or metallicity == 0.0001):
            ratio = np.array([22.20, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                              0.489, 0.295, 0.203, 0.270, 0.015, 0.005, 0.002,
                              0.256, 0.109, 0.036, 0.010, 1.000, 1.097, 3.159,
                              0.003, 0.000, 0.096, 0.008, 0.009, 0.005, 2.870,
                              0.015, 0.026, 0.037, 0.029, 0.028, 0.027, 0.012,
                              0.007, 0.067, 0.000, 0.000, 0.000, 0.000, 0.000])
        # ratio for metallicity = 0.004
        if (metallicity == 0.004):
            ratio = np.array([22.20, 0.000, 0.058, 0.000, 0.000, 0.000, 0.310,
                              1.791, 0.416, 0.192, 0.283, 0.015, 0.017, 0.007,
                              0.256, 0.066, 0.036, 0.014, 1.000, 1.617, 4.752,
                              0.010, 0.000, 0.108, 0.041, 0.017, 0.059, 2.870,
                              0.175, 0.030, 0.188, 0.138, 0.023, 0.071, 0.027,
                              0.014, 0.176, 0.510, 0.000, 0.000, 0.000, 0.000])
        # ratio for metallicity = 0.008 or 0.02 or 0.05
        if (metallicity == 0.008 or metallicity == 0.02 or
                metallicity == 0.05):
            ratio = np.array([22.20, 0.110, 0.010, 0.180, 0.010, 0.290, 0.070,
                              3.010, 0.300, 0.107, 0.159, 0.015, 0.029, 0.011,
                              0.256, 0.010, 0.050, 0.000, 1.000, 1.399, 4.081,
                              0.030, 0.010, 0.140, 0.130, 0.030, 0.136, 2.870,
                              0.404, 0.030, 0.300, 0.210, 0.040, 0.035, 0.026,
                              0.014, 0.086, 0.945, 0.365, 0.048, 0.058, 0.054])


            v = 300000  # m.s^-1
            fwhm = (wav * v) / c  # line width depending on z

            #x=n.arange(1100,10500) # range of wavelength

            sigma = fwhm / (2*np.sqrt(2*np.log(2)))

            # attenuation linked to the distance

            intensity = ratio*flow_H

            sed.add_lines(self.parameters["set_name"], wav, intensity, sigma)
