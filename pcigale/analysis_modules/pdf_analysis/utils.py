# -*- coding: utf-8 -*-
# Copyright (C) 2013 Centre de données Astrophysiques de Marseille
# Copyright (C) 2014 Yannick Roehlly <yannick@iaora.eu>
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

import numpy as np
from copy import deepcopy
from ...sed.cosmology import cosmology
from ...creation_modules import get_module as get_creation_module


def gen_compute_fluxes_at_redshift(sed, filters, redshifting_module,
                                   igm_module):
    """"Generate function to compute the fluxes of a SED at a given redshift

    Given a SED, a list of filters and a redshift module name, this generator
    returns a function computing the fluxes of the SED in all the filters at
    a given redshift.  If the SED is older than the Universe at the given
    redshift, the fluxes returned by this function will be all -99.

    If the redshift module name is None, the SED is not redshifted by the
    returned function. This means that it will always return the same fluxes,
    whatever its redshift parameter. This can be used for computing
    photometric redshifts.

    Parameters
    ----------
    sed : pcigale.sed
        The pcigale SED object.
    filters : list of pcigale.data.filters
        List of pcigale filters objects.
    redshifting_module : picgale.creation_modules.Module
        A pcigale SED creation module to redshift the SED. It must accept a
        redshift parameter (or be None).
    igm_module : picgale.creation_modules.Module
        A pcigale SED creation module to add the IGM attenuation to the SED.
        Is not used if the redshifting module is None.

    Return
    ------
    gen_fluxes : function
        Function computing the fluxes of the SED in all filters at a given
        redshift.

    """
    # The returned function is memoized.
    cache = {}

    def gen_fluxes(redshift):
        """Compute the flux of the SED in various filters.

        Parameters
        ----------
        redshift : float

        Returns
        -------
        array fo floats

        """

        # If the function is generated without a redshift module, it always
        # computes the fluxes at the SED redshift (the SED may be already
        # redshifted).
        if not redshifting_module:
            redshift = sed.redshift

        if redshift not in cache:
            # Age of the Universe at redshift.
            # astropy 0.3 cosmology functions return quantities
            try:
                age_at_redshift = cosmology.age(redshift).value * 1000
            except AttributeError:
                age_at_redshift = cosmology.age(redshift) * 1000

            if sed.info["age"] > age_at_redshift:
                cache[redshift] = -99 * np.ones(len(filters))
            else:

                if redshifting_module:
                    red_sed = deepcopy(sed)
                    redshifting_module.parameters["redshift"] = redshift
                    redshifting_module.process(red_sed)
                    if igm_module:
                        igm_module.process(red_sed)
                else:
                    red_sed = sed

                cache[redshift] = np.array(
                    [red_sed.compute_fnu(filt_.trans_table,
                                         filt_.effective_wavelength)
                     for filt_ in filters])

        return cache[redshift]

    return gen_fluxes
