# -*- coding: utf-8 -*-
"""
Copyright (C) 2012 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""


from . import common
from pcigale.data import Database


class Module(common.SEDCreationModule):
    """
    Module computing the infra-red re-emission corresponding to an amount of
    attenuation using the Dale and Helou (2002) templates.

    Given an amount of attenuation (e.g. resulting from the action of a dust
    extinction module) this module normalises the Dale and Helou (2002)
    template corresponding to a given α to this amount of energy and add it
    to the SED.

    Information added to the SED: NAME_alpha.

    """

    parametre_list = {
        'alpha': (
            'float',
            None,
            "Alpha slope.",
            None
        ),
        'extinction_value_names': (
            'array of strings',
            None,
            "List of extinction value names (in the SED's info dictionary). "
            "A new re-emission contribution will be added for each one.",
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
        extinction_value_names = parametres["extinction_value_names"]

        # Get the template set out of the database
        database = Database()
        dh2002 = database.get_dh2002_infrared_templates()
        database.session.close_all()

        ir_template = dh2002.get_template(alpha)

        # Base name for adding information to the SED.
        name = self.name or 'dh2002_ir'

        sed.add_module(name, parametres)
        sed.add_info(name + '_alpha', alpha)

        # We multiply the extinction * ir_template by -1 because the value
        # of the former is positive in the SED's info and its effects are
        # negative.
        for extinction in extinction_value_names:
            sed.add_contribution(
                name + '_' + extinction,
                dh2002.wavelength_grid,
                -1 * sed.info[extinction] * ir_template
            )
