# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly


class Fritz2006(object):
    """Fritz et al. (2006) AGN dust torus emission model.

    This class holds the UV-optical data associated with a Fritz et al. (2006)
    AGN model.

    """

    def __init__(self, r_ratio, tau, beta, gamma, opening_angle, psy, wave,
                 lumin_therm, lumin_scatt, lumin_agn):
        """Create a new AGN model

        Parameters
        ----------
        r_ratio: float
            Ratio of the maximum and minimum radii of the dust torus.
        tau: float
            Tau at 9.7µm
        beta: float
            Beta
        gamma: float
            Gamma
        opening_angle: float
            Opening angle of the dust torus.
        psy: float
            Angle between AGN axis and line of sight.
        """

        self.r_ratio = r_ratio
        self.tau = tau
        self.beta = beta
        self.gamma = gamma
        self.opening_angle = opening_angle
        self.psy = psy
        self.wave = wave
        self.lumin_therm = lumin_therm
        self.lumin_scatt = lumin_scatt
        self.lumin_agn = lumin_agn
