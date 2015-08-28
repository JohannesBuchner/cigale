# -*- coding: utf-8 -*-
# Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

"""
Read spectrum from file module
==============================

This module reads a SED spectrum from a file.

"""

from astropy.table import Table
from ..utils import read_table
from . import CreationModule


class LoadSpecFile(CreationModule):
    """Module reading a spectrum from a file and adding it to the SED.

    """

    parameter_list = dict([
        ("filename", (
            'str',
            "Name of the file to load and to add to the SED table. This "
            "file must be loadable with astropy",
            None
        )),
        ("lambda_column", (
            'str',
            "Name of the column containing the wavelength in nm.",
            None
        )),
        ("l_lambda_column", (
            'str',
            "Name of the column containing the Lλ luminosity in W/nm.",
            None
        ))
    ])

    def process(self, sed):
        """Add the spectrum from the file to the SED object

        Parameters
        ----------
        sed: pcigale.sed.SED object

        """
        filename = self.parameters['filename']
        table = read_table(filename)

        sed.add_module(self.name, self.parameters)

        sed.add_contribution(
            filename,
            table[self.parameters['lambda_column']],
            table[self.parameters['l_lambda_column']]
        )

# CreationModule to be returned by get_module
Module = LoadSpecFile
