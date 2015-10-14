# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2013 Institute of Astronomy, University of Cambridge
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Laure Ciesla


class Schreiber2016(object):
    """Schreiber et al. (2016) dust models

    This class holds the data associated with the Schreiber et al. (2016)
    dust models.

    """

    def __init__(self, type, tdust, wave, lumin):
        """Create a new IR model

        Parameters
        ----------
        type: integer
            Type of the SED: PAH or dust continuum
        tdust: float
            Dust temperature
        wave: array
            Vector of the λ grid used in the templates [nm]
        lumin: array
            Model data in an array containing the luminosity density
            versus the wavelength λ

        """

        self.type = type
        self.tdust = tdust
        self.wave = wave
        self.lumin = lumin
