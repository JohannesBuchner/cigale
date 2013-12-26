# -*- coding: utf-8 -*-
# Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly <yannick.roehlly@oamp.fr>

import atpy
from collections import OrderedDict
from . import common


class Module(common.SEDCreationModule):
    """Module reading a spectrum from a file and adding it to the SED.

    Note that this module uses the atpy module, which is not automatically
    installed when one installs pcigale.

    """

    parameter_list = OrderedDict([
        ("filename", (
            'str',
            "Name of the file to load and to add to the SED table. This "
            "file must be loadable with atpy (that depends on other modules "
            "being installed).",
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
        sed  : pcigale.sed.SED object

        """
        filename = self.parameters['filename']
        table = atpy.Table(filename, verbose=False)

        sed.add_module(self.name, self.parameters)

        sed.add_contribution(
            filename,
            table[self.parameters['lambda_column']],
            table[self.parameters['l_lambda_column']]
        )
