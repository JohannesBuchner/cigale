# -*- coding: utf-8 -*-
"""
Copyright (C) 2012 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""


import numpy as np
from . import common
from pcigale.data import Database


class Module(common.SEDCreationModule):
    """
    Compute the infra-red re-emission of absorbed energy using Dale and Helou
    (2002) templates.

    This module integrates the energy absorbed by the dust and normalises the
    Dale and Helou (2002) template corresponding to the given α to this amount
    of energy. This module can only be used on a SED that has gone through a
    dust extinction module that has produced at least one extinction
    (negative) contribution.

    """

    parametre_list = {
        'alpha': (
            'float',
            None,
            "Alpha slope.",
            None
        ),
        'extinction_contrib_names': (
            'array of strings',
            None,
            "List of the extinction contributions to process. This module "
            "will a new re-emission contribution for each.",
            None
        )
    }

    def _process(self, sed, parametres):
        """Add the IR re-emission contributions

        Parametres
        ----------
        sed  : pcigale.sed.SED object
        parametres : dictionnary containing the parametres

        """
        alpha = parametres["alpha"]
        extinction_contrib_names = parametres["extinction_contrib_names"]

        # Get the template set out of the database
        database = Database()
        dh2002 = database.get_dh2002_infrared_templates()
        database.session.close_all()

        ir_template = dh2002.get_template(alpha)

        sed.add_info('ir_alpha_slope', alpha)

        for contrib_name in extinction_contrib_names:
            # The extinction contribution is negative
            extinction = -1. * np.trapz(
                sed.get_lumin_contribution(contrib_name),
                sed.wavelength_grid
            )
            sed.add_component(
                'dh2002_ir',
                parametres,
                'dh2002_ir_' + contrib_name,
                dh2002.wavelength_grid,
                extinction * ir_template,
                {}
            )
