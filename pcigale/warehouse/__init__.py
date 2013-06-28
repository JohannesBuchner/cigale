# -*- coding: utf-8 -*-
"""
Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""
from json import JSONEncoder
from ..sed import SED
from ..sed.modules import common as sed_modules


class SedWarehouse(object):
    """Create, cache and store SED

    This object is responsible for creating SED and storing them in a memory
    cache or a database.
    """

    def __init__(self, cache_type="memory"):
        """Instantiate a SED warehouse

        Parameters
        ----------
        cache_type : string
            Type of cache used. For now, only in memory caching.
        """
        if cache_type == "memory":
            from .store.memory import SedStore
        elif cache_type == "shelf":
            from .store.shelf import SedStore

        self.storage = SedStore()

    def get_sed(self, module_list, parameter_list):
        """Get the SED corresponding to the module and parameter lists

        If the SED was cached, get it from the cache. If it is not, create it
        and add it the the cache. The method is recursive to permit caching
        partial SED.

        Parameters
        ----------
        module_list : iterable
            List of module names in the order they have to be used to
            create the SED.
        parameter_list : iterable
            List of the parameter dictionaries corresponding to each
            module of the module_list list.

        Returns
        -------
        sed : pcigale.sed
            The SED made from the given modules with the given parameters.

        """
        module_list = list(module_list)
        parameter_list = list(parameter_list)

        # JSon representation of the tuple (module_list, parameter_list)
        # used as a key for storing the SED in the cache.
        encoder = JSONEncoder()
        sed_key = encoder.encode((module_list, parameter_list))

        sed = self.storage.get(sed_key)

        if not sed:
            mod = sed_modules.get_module(module_list.pop(),
                                         **parameter_list.pop())

            if (len(module_list) == 0):
                sed = SED()
            else:
                sed = self.get_sed(module_list, parameter_list)

            mod.process(sed)
            self.storage.add(sed_key, sed)

        return sed

    def close(self):
        """ Close the underlying storage if needed """
        self.storage.close()
