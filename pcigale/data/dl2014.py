## -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013 Institute of Astronomy, University of Cambridge
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Médéric Boquien


class DL2014(object):
    """2014 update of the Draine and Li (2007) dust models, 

    This class holds the data associated with the updated Draine and Li (2007)
    dust models.

    """

    def __init__(self, qpah, umin, umax, alpha, wave, lumin):
        """Create a new IR model

        Parameters
        ----------
        qpah: float
            Mass fraction of PAH
        umin: float
            Minimum radiation field illuminating the dust
        umax: float
            Maximum radiation field illuminating the dust
        alpha: float
            Powerlaw slope dU/dM∝U¯ᵅ
        wave: array
            Vector of the λ grid used in the templates [nm]
        lumin: array
            Model data in an array containing the luminosity density
            versus the wavelength λ

        """

        self.qpah = qpah
        self.umin = umin
        self.umax = umax
        self.alpha = alpha
        self.wave = wave
        self.lumin = lumin
