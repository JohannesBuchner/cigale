# -*- coding: utf-8 -*-
"""
Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
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
    attenuation module) this module normalises the Dale and Helou (2002)
    template corresponding to a given α to this amount of energy and add it
    to the SED.

    Information added to the SED: NAME_alpha.

    """

    parameter_list = {
        'alpha': (
            'float',
            "Alpha slope.",
            None
        ),
        'attenuation_value_names': (
            'array of strings',
            "List of attenuation value names (in the SED's info dictionary). "
            "A new re-emission contribution will be added for each one.",
            None
        )
    }

    out_parameter_list = {'alpha': 'Alpha slope.'}

    def _init_code(self):
        """Get the template set out of the database"""
        database = Database()
        self.dh2002 = database.get_dh2002_infrared_templates()
        database.session.close_all()

    def _process(self, sed, parameters):
        """Add the IR re-emission contributions

        Parameters
        ----------
        sed  : pcigale.sed.SED object
        parameters : dictionary containing the parameters

        """
        alpha = float(parameters["alpha"])
        attenuation_value_names = parameters["attenuation_value_names"]

        dh2002 = self.dh2002
        ir_template = dh2002.get_template(alpha)

        # Base name for adding information to the SED.
        name = self.name or 'dh2002'

        sed.add_module(name, parameters)
        sed.add_info(name + '_alpha', alpha)

        for attenuation in attenuation_value_names:
            sed.add_contribution(
                name + '_' + attenuation,
                dh2002.wavelength_grid,
                sed.info[attenuation] * ir_template
            )
