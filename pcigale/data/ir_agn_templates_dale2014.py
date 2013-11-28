## -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly <yannick.roehlly@oamp.fr>, Laure Ciesla <ciesla@physics.uoc.gr>


class DALE2014(object):
    """Dale et al (2014) IR templates containing an AGN component.

    This class holds the data associated with the Dale et al (2014)
    dust templates.

    """

    def __init__(self, fracAGN, alpha, wave, lumin):
        """Create a new IR model

        Parameters
        ----------
        fracAGN: float
            Contribution of the AGN
        alpha: float
            Dale & Helou (2002) alpha slope.
        lir: float
            Total IR luminosity between 8 and 1000 microns (AGN + SB)
        wave: array
            Vector of the λ grid used in the templates [nm]
        lumin: array
            Model data in an array containing the luminosity density
            versus the wavelength λ

        """

        self.fracAGN = fracAGN
        self.alpha = alpha
        self.wave = wave
        self.lumin = lumin
        