# -*- coding: utf-8 -*-
"""
Copyright (C) 2012 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""
from . import SED
from .modules import common as sed_modules


def create_sed(module_list):
    """Create a new SED using the given modules

    Parametres
    ----------
    module_list :array of tuples
        Array of tuples (module name, parametre dictionnary).

    Returns
    -------
    sed : pcigale.sed
        The SED made from the given modules with the given parametres.

    """

    # We start from an empty SED.
    sed = SED()

    for (module, parametres) in module_list:
        mod = sed_modules.get_module(module)
        mod.parametres = parametres
        mod.process(sed)

    return sed
