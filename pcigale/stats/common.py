# -*- coding: utf-8 -*-
"""
Copyright (C) 2012 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""


class AnalysisModule(object):
    """Abstract class, the pCigale analysis modules are based on.
    """

    # parametre_list is a dictionnary containing all the parametres used by
    # the module. Each parametre name is associate to a tuple (variable type,
    # unit [string], description [string], default value). Each module must
    # define its parametre list, unless it does not use any parametre. Using
    # None means that there is no description, unit or default value. If None
    # should be the default value, use the 'None' string instead.
    parametre_list = {}

    def __init__(self, **kwargs):
        """Instantiate a analysis module

        The module parametres values can be passed as keyworded paramatres.
        """
        # parametres is a dictionnary containing the actual values for each
        # module parametre.
        self.parametres = kwargs


def get_module(module_name):
    """Return the main class of the module provided

    Parametres
    ----------
    module_name : string
        The name of the module we want to get the class.

    Returns
    -------
    module_class : class
    """

    try:
        # TODO Find a better way to do dynamic import
        import_string = 'from . import ' + module_name + ' as module'
        exec import_string
        return module.Module()
    except ImportError:
        print('Module ' + module_name + ' does not exists!')
        raise
