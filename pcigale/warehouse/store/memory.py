# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

from copy import deepcopy


class SedStore(object):
    """In memory storage"""

    def __init__(self):
        self.dictionary = {}

    def get(self, key):
        """Get a value from the cache dictionary

        Parameters
        ----------
        key: any immutable

        Returns
        -------
        object

        """
        # We return a copy not to modify the stored object.
        return deepcopy(self.dictionary.get(key))

    def add(self, key, value):
        """Add a new key, value pair to the cache.

        Parameters
        ----------
        key: any immutable
        value: object

        """
        # We store a copy not to modify the stored object.
        self.dictionary[key] = deepcopy(value)

    def delete(self, key):
        """Delete a key, value pair from the cache

        Parameters
        ----------
        key: key of the element to be deleted

        """
        del self.dictionary[key]

    def close(self):
        """Do nothing"""
        pass
