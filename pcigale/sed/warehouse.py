# -*- coding: utf-8 -*-
"""
Copyright (C) 2012 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""
from . import SED
from .modules import common as sed_modules


def create_sed(module_list, parametre_list):
    """Create a new SED using the given modules and parametres

    Parametres
    ----------
    module_list : list
        List of module names in the order they have to be used to create
        the SED.
    parametre_list : list of dictionaries
        List of the parametre dictionaries corresponding to each module of
        the module_list list.

    Returns
    -------
    sed : pcigale.sed
        The SED made from the given modules with the given parametres.

    """

    # We start from an empty SED.
    sed = SED()

    for (module, parametres) in zip(module_list, parametre_list):
        mod = sed_modules.get_module(module)
        mod.parametres = parametres
        mod.process(sed)

    return sed
