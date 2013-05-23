# -*- coding: utf-8 -*-
"""
Copyright (C) 2013 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""

import numpy as np


class AgnFritz2006(object):
    """Fritz et al. (2006) AGN dust torus emision model.

    This class holds the UV-optical data associated with a Fritz et al. (2006)
    AGN model.

    """

    def __init__(self, model_nb, agn_type, r_ratio, tau, beta, gamma, theta,
                 psy, wave, luminosity):
        """Create a new AGN model

        Parameters
        ----------
        model_nb : integer
            Number identifying the AGN model.
        agn_type : integer
            Type of AGN.
        r_ratio : float
            Ratio of the maximum and minimum radii of the dust torus.
        tau : float
            Tau at 9.7µm
        beta : float
            Beta
        gamma : float
            Gamma
        theta : float
            Opening angle of the dust torus.
        psy : float
            Angle between AGN axis and line of sight.
        wave : array of float
            Wavelength grid in nm.
        luminosity : array of float
            Luminosity density at each wavelength in W/nm.

        """
        self.model_nb = model_nb
        self.agn_type = agn_type
        self.r_ratio = r_ratio
        self.tau = tau
        self.beta = beta
        self.gamma = gamma
        self.theta = theta
        self.psy = psy
        self.wave = np.array(wave)
        self.luminosity = np.array(luminosity)
