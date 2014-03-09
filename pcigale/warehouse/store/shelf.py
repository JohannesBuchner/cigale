# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

import shelve
from copy import deepcopy


class SedStore(object):
    """Shelf (dictionary into a file) storage"""

    def __init__(self, filename="pcigale.cache"):
        self.shelf_storage = shelve.open(filename)

    def get(self, key):
        """Get a value from the shelf

        Parameters
        ----------
        key : any immutable

        Returns
        -------
        object

        """
        # We return a copy not to modify the stored object.
        return deepcopy(self.shelf_storage.get(key))

    def add(self, key, value):
        """Add a new key, value pair to the shelf.

        Parameters
        ----------
        key : any immutable
        value : object

        """
        # We store a copy not to modify the stored object.
        self.shelf_storage[key] = deepcopy(value)

    def delete(self, key):
        """Delete a key, value pair from the cache

        Parameters
        ----------
        key: key of the element to be deleted

        """
        del self.dictionary[key]

    def close(self):
        """Close the shelf file"""
        self.shelf_storage.close()
