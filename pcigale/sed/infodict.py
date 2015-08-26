# -*- coding: utf-8 -*-
# Copyright (C) 2015 Institute of Astronomy
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Médéric Boquien

"""
This class is the implementation of a light weight ordered-dictionary-like
object. Its main advantage is that it is very fast to copy, unlike
OrderedDict(), which can take a significant fraction of the runtime.

It makes some very strong assumptions on the data and how the object is
used. Because of this, it should not be used for anything else other than the
info member of the SED class.

The implementation is based on the emulation of a container type and therefore
gives some standard function to behave like a dictionary, making it a drop-in
replacement. The class members are self-explanatory.

"""


class InfoDict(object):
    """Light weight ordered-dictionary-like object that is fast to copy.
    """

    def __init__(self, keys=[], values=[]):
        self.k = keys[:]
        self.v = values[:]

    def __len__(self):
        return len(self.k)

    def __getitem__(self, key):
        try:
            idx = self.k.index(key)
            return self.v[idx]
        except ValueError:
            raise KeyError

    def __setitem__(self, key, value):
        if key not in self.k:
            self.k.append(key)
            self.v.append(value)
        else:
            self.v[self.k.index(key)] = value

    def __contains__(self, key):
        if key in self.k:
            return True
        else:
            return False

    def keys(self):
        return self.k[:]

    def values(self):
        return self.v[:]

    def items(self):
        return ((k, v) for k, v in zip(self.k, self.v))

    def copy(self):
        return InfoDict(self.k, self.v)
