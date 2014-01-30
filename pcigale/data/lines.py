## -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2014 Institute of Astronomy, University of Cambridge
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Médéric Boquien


class Lines(object):
    """Nebular lines templates.

    This class holds the data associated with the line templates

    """

    def __init__(self, metallicity, logU, wave, ratio):
        """Create a new nebular lines template

        Parameters
        ----------
        metallicity: float
            Gas phase metallicity
        logU: float
            Radiation field intensity
        wave: array
            Vector of the λ grid used in the templates [nm]
        ratio: array
            Line intensities relative to Hβ

        """

        self.metallicity = metallicity
        self.logU = logU
        self.wave = wave
        self.ratio = ratio
