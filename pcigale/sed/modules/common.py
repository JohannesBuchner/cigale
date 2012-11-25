# -*- coding: utf-8 -*-
"""
Copyright (C) 2012 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""


class SEDCreationModule(object):
    """Abstract class, the pCigale SED creation modules are based on.
    """

    # parametre_list is a dictionnary containing all the parametres used by the
    # module. Each parametre name is associate to a tuple (variable type, unit
    # [string], description [string], default value). Each module must define
    # its parametre list, unless it does not use any parametre. Using None
    # means that there is non description, unit or default value. If None
    # should be the default value, use the 'None' string instead.
    parametre_list = {}

    def __init__(self, **kwargs):
        """Instantiate a SED creation module

        The module parametres values can be passed as keyworded paramatres.
        """
        # parametres is a dictionnary containing the actual values for each
        # module parametre.
        self.parametres = kwargs

    def printParametreList(self):
        """Pretty print the list of parametres for the module
        """
        for parametre in self.parametre_list:
            pType, unit, description, defValue = \
                self.parametre_list[parametre]
            outStr = str(parametre) + ': ' + str(description)
            if pType:
                outStr += ' (type=' + str(pType) + ')'
            if unit:
                outStr += ' (unit=' + str(unit) + ')'
            if defValue:
                outStr += ' (default=' + str(defValue) + ')'
            print outStr

    def _process(self, sed, parametres):
        """Do the actual processing of the module on a SED object

        This method is called with an object and a complete parametres
        dictionnary. It is not meant to be called directly but through the
        process method.

        Parametres
        ----------
        sed  : pcigale.sed.SED object
        parametres : dictionnary of the module parametres

        """
        raise NotImplementedError()

    def process(self, sed, parametres=None):
        """Process a SED object with the module

        This method is responsible for checking the module parametres (whether
        they are given in the method call or are taken from parametres class
        attribute) before doing the actual processing (_process method). If a
        parametre is not given but exists in the parametre_list with a default
        value, this value is used.

        The SED object is updated during the process, one must take care of
        copying it before, if needed.

        Parametres
        ----------
        sed  : pcigale.sed.SED object
        parametres : dictionnary
            Dictionnary containing the module parametre values, if it is not
            given, the module parametre values are used

        Raises
        ------
        KeyError : when not all the needed parametres are given.

        """

        # If the parametre dictionnary is not passed, use the module one
        if not parametres:
            parametres = self.parametres

        # For parametres that are present on the parametre_list with a default
        # value and that are not in the parametres dictionnary, we add them
        # with their default value.
        for key in self.parametre_list:
            if (not key in parametres) and (
                    self.parametre_list[key][3] is not None):
                parametres[key] = self.parametre_list[key][3]

        # If the keys of the parametres dictionnary are different from the one
        # of the parametre_list dictionnary, we raises a KeyError. That means
        # that a parametre is missing (and has no default value) or that an
        # unexpected one was given.
        if not parametres.keys() == self.parametre_list.keys():
            raise KeyError("The parametres passed are different from the "
                           "expected one.")

        # TODO: We should also check that all parametres is from the right
        # type.

        # We do the actual processing.
        self._process(sed, parametres)


def getModuleClass(moduleName):
    """Return the main class of the module provided

    Parametres
    ----------
    moduleName : string
        The name of the module we want to get the class.

    Returns
    -------
    moduleClass : class
    """

    try:
        # TODO Find a better way to do dynamic import
        importString = 'from . import ' + moduleName + ' as module'
        exec importString
        return module.Module()
    except ImportError:
        print('Module ' + moduleName + ' does not exists!')
        raise
