# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2014 Institute of Astronomy, University of Cambridge
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Médéric Boquien


class NebularContinuum(object):
    """Nebular lines templates.

    This class holds the data associated with the line templates

    """

    def __init__(self, metallicity, logU, wave, lumin):
        """Create a new nebular lines template

        Parameters
        ----------
        metallicity: float
            Gas phase metallicity
        logU: float
            Ionisation parameter
        wave: array
            Vector of the λ grid used in the templates [nm]
        lumin: array
            Luminosity density of the nebular continuum in Fλ

        """

        self.metallicity = metallicity
        self.logU = logU
        self.wave = wave
        self.lumin = lumin
