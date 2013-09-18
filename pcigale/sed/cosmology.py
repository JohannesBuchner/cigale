# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly <yannick.roehlly@oamp.fr>
"""Cosmology configuration

We use astropy.cosmology module. For now, we import the desired cosmology here
so that it can be changed easily for people needing a specific cosmology.

"""

from astropy.cosmology import WMAP7 as cosmology
