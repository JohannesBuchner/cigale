# -*- coding: utf-8 -*-
"""
Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""
from copy import deepcopy
from json import JSONEncoder
from . import SED
from .modules import common as sed_modules


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
            dict_cache = {}

            def get_from_cache(key):
                """Return the value corresponding to the key in the cache.

                If the key is not in the cache, returns None.

                Parameters
                ----------
                key : any immutable

                Returns
                -------
                object

                """
                # We return a copy not to modify the stored object.
                return deepcopy(dict_cache.get(key))

            def add_to_cache(key, value):
                """Add a new key, value pair to the cache.

                Parameters
                ----------
                key : any immutable
                value : object

                """
                # We store a copy not to modify the stored object.
                dict_cache[key] = deepcopy(value)

        self.get_from_cache = get_from_cache
        self.add_to_cache = add_to_cache

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

        sed = self.get_from_cache(sed_key)

        if not sed:
            mod = sed_modules.get_module(module_list.pop())
            mod.parameters = parameter_list.pop()

            if (len(module_list) == 0):
                sed = SED()
            else:
                sed = self.get_sed(module_list, parameter_list)

            mod.process(sed)
            self.add_to_cache(sed_key, sed)

        return sed
