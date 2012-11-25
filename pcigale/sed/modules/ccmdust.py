# -*- coding: utf-8 -*-
"""
Copyright (C) 2012 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""


from . import common
from ...extern.lsst import Sed as lsst


class Module(common.SEDCreationModule):
    """Add CCM dust model extinction to the SED

    If a contribution name is given in the parametre list, the extinction will
    be applied only to the flux of this contribution; else, it will be applied
    to the whole spectrum.

    This module is based on the code from the Large Synoptic Survey Telescope.
    http://dev.lsstcorp.org/trac/
    """

    parametre_list = {
        'A_v': (
            'float',
            'kesako',
            'kesako',
            'None'
        ),
        'ebv': (
            'float',
            'kesako',
            'kesako',
            'None'
        ),
        'R_v': (
            'float',
            'kesako',
            'kesako',
            3.1
        ),
        'contribution_name': (
            'string',
            None,
            'Name of the contribution the dust extinction will be applied '
            'to. If None, it will be applied to the whole spectrum.',
            'None'
        )

    }

    def _process(self, sed, parametres):
        """Add CCM dust model extinction to the SED

        The computation is done by LSST functions. Two on the three parametres
        (A_v, ebv and R_v) are expected. If that's not the case, the LSST
        function will raise an exception.

        Parametres
        ----------
        sed : pcigale.sed.SED object
        parametres : dictionnary containing the parametres

        """

        # First, we set the parametres that have 'None' (string) as value to
        # None.
        for key in parametres:
            if parametres[key] == 'None':
                parametres[key] = None

        wavelen = sed.wavelength_grid

        # We get either to contribution flux or the whole spectrum.
        if parametres['contribution_name']:
            l_lambda = sed.get_lumin_contribution(
                parametres['contribution_name']
            )
        else:
            l_lambda = sed.luminosity

        # We need a lsst.Sed object to use its methods.
        lsstSed = lsst.Sed()

        a_x, b_x = lsstSed.setupCCMab(wavelen=wavelen)
        wavelen, extended_l_lambda = lsstSed.addCCMDust(
            a_x,
            b_x,
            A_v=parametres['A_v'],
            ebv=parametres['ebv'],
            R_v=parametres['R_v'],
            wavelen=wavelen,
            flambda=l_lambda
        )

        # We only want to add the extinction as a SED component. We compute
        # the difference because extended_flambda and flambda have the same
        # wavelength grit (see addCCMDust definition).
        extinction = extended_l_lambda - l_lambda

        # If the extinction was applied to a specific contribution flux, we
        # suffix the ccmdust contribution name with its own.
        ccmdust_contrib_name = 'ccmdust'
        if parametres['contribution_name']:
            ccmdust_contrib_name += '_' + parametres['contribution_name']

        sed.add_component(
            'ccmdust',
            parametres,
            ccmdust_contrib_name,
            wavelen,
            extinction,
            {}
        )
