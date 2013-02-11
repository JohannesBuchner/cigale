# -*- coding: utf-8 -*-
"""
Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""


import numpy as np
from . import common
from ...extern.lsst import Sed as lsst


class Module(common.SEDCreationModule):
    """Add CCM dust model extinction to the SED

    If a contribution name is given in the parameter list, the extinction will
    be applied only to the flux of this contribution; else, it will be applied
    to the whole spectrum.

    This module is based on the code from the Large Synoptic Survey Telescope.
    http://dev.lsstcorp.org/trac/

    The parameters added and available for the statistical analysis are:
    ccmdust_A_v, ccmdust_ebv and ccmdust_R_v.

    """

    parameter_list = {
        'A_v': (
            'float',
            'kesako',
            'None'
        ),
        'ebv': (
            'float',
            'kesako',
            'None'
        ),
        'R_v': (
            'float',
            'kesako',
            3.1
        ),
        'contribution_name': (
            'string',
            'Name of the contribution the dust extinction will be applied '
            'to. If None, it will be applied to the whole spectrum.',
            'None'
        )

    }
    comments = ("One must indicate only two of the three parameters 'A_v', "
                "'ebv' and 'R_v'. The other must be set to 'None'.")

    def _process(self, sed, parameters):
        """Add CCM dust model extinction to the SED

        The computation is done by LSST functions. Two on the three parameters
        (A_v, ebv and R_v) are expected. If that's not the case, the LSST
        function will raise an exception.

        Parameters
        ----------
        sed : pcigale.sed.SED object
        parameters : dictionary containing the parameters

        """

        # First, we set the parameters that have 'None' (string) as value to
        # None.
        for key in parameters:
            if parameters[key] == 'None':
                parameters[key] = None

        wavelen = sed.wavelength_grid

        # We get either to contribution flux or the whole spectrum.
        if parameters['contribution_name']:
            l_lambda = sed.get_lumin_contribution(
                parameters['contribution_name']
            )
        else:
            l_lambda = sed.luminosity

        # We need a lsst.Sed object to use its methods.
        lsstSed = lsst.Sed()

        a_x, b_x = lsstSed.setupCCMab(wavelen=wavelen)
        wavelen, extended_l_lambda = lsstSed.addCCMDust(
            a_x,
            b_x,
            A_v=parameters['A_v'],
            ebv=parameters['ebv'],
            R_v=parameters['R_v'],
            wavelen=wavelen,
            flambda=l_lambda
        )

        # We only want to add the extinction as a SED component. We compute
        # the difference because extended_flambda and flambda have the same
        # wavelength grid (see addCCMDust definition).
        extinction = extended_l_lambda - l_lambda

        # Integrate the value of the extinction (-1 is because the extinction
        # spectrum is negative.
        extinction_value = -1 * np.trapz(extinction, wavelen)

        # Base name for adding information to the SED.
        name = self.name or 'ccmdust'

        sed.add_module(name, parameters)

        # Add the parameters values to the SED information.
        sed.add_info(name + '_A_v', parameters['A_v'])
        sed.add_info(name + '_ebv', parameters['ebv'])
        sed.add_info(name + '_R_v', parameters['R_v'])

        # Add the extinction value to the SED information
        sed.add_info(name + '_extinction', extinction_value)

        sed.add_contribution(
            name,
            wavelen,
            extinction
        )
