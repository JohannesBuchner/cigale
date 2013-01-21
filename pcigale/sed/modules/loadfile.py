# -*- coding: utf-8 -*-
"""
Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""


import atpy
from . import common


class Module(common.SEDCreationModule):
    """Module reading a spectrum from a file and adding it to the SED.

    Note that this module uses the atpy module, which is not automatically
    installed when one installs pcigale.

    """

    parametre_list = {
        "filename": (
            'str',
            "Name of the file to load and to add to the SED table. This "
            "file must be loadable with atpy (that depends on other modules "
            "being installed).",
            None
        ),
        "lambda_column": (
            'str',
            "Name of the column containing the wavelength in nm.",
            None
        ),
        "l_lambda_column": (
            'str',
            "Name of the column containing the Lλ luminosity in W/nm.",
            None
        )
    }

    def _process(self, sed, parametres):
        """Add the spectrum from the file to the SED object

        Parametres
        ----------
        sed  : pcigale.sed.SED object
        parametres : dictionnary containing the parametres

        """
        filename = parametres['filename']
        table = atpy.Table(filename)

        # Base name for adding information to the SED.
        name = self.name or 'loadfile'

        sed.add_module(name, parametres)

        sed.add_contribution(
            name + '_' + filename,
            table[parametres['lambda_column']],
            table[parametres['l_lambda_column']]
        )
