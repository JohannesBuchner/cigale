# -*- coding: utf-8 -*-
# Copyright (C) 2014 Yannick Roehlly
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

"""
Redshifting module
==================

This module implements the redshifting of a SED. The SED must be rest-frame
or the module will raise en exception when processing it.

Note that this module, contrary to the other SED creation modules, actually
changes the individual luminosity contributions as it redshifts everyone.
Also note that doing this, this module does not rely on the SED object
interface but on its inner implementations. That means that if the SED object
is changed, this module may need to be adapted.

"""

from collections import OrderedDict
from ..creation_modules import CreationModule


class Redshifting(CreationModule):
    """Redshift a SED

    This module redshift a rest-frame SED. If the SED is already redshifted, an
    exception is raised.

    """

    parameter_list = OrderedDict([
        ("redshift", (
            "float",
            "Redshift to apply to the galaxy.",
            None
        ))
    ])

    def process(self, sed):
        """Redshift the SED

        Parameters
        ----------
        sed : pcigale.sed.SED object

        """
        redshift = float(self.parameters["redshift"])

        # If the SED is already redshifted, raise an error.
        if sed.redshift > 0:
            raise Exception("The SED is already redshifted <z={}>."
                                .format(sed.redshift))

        # Raise an error when applying a negative redshift. This module is
        # not for blue-shifting.
        if redshift < 0:
            raise Exception("The redshift provided is negative <{}>."
                                .format(redshift))

        # We redshift directly the SED wavelength grid
        sed.wavelength_grid *= 1. + redshift

        # We modify each luminosity contribution to keep energy constant
        sed.luminosities /= 1. + redshift

        sed.redshift = redshift
        sed.add_module(self.name, self.parameters)

# CreationModule to be returned by get_module
Module = Redshifting
