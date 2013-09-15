# -*- coding: utf-8 -*-
# Copyright (C) 2012 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly <yannick.roehlly@oamp.fr>

import numpy as np
from scipy import interpolate


class IrTemplatesDH2002(object):
    """Dale an Helou 2002 infra-red templates

    This class holds the data associated with the Dale and Helou (2002)
    infra red templates.

    """

    def __init__(self, alpha_grid, wavelength_grid, templates):
        """Create a new IR template set

        Parameters
        ----------
        alpha_grid : array
            Vector of the various values for the α slope in the templates.
        wavelength_grid : array
            Vector of the λ grid used in the templates [nm]
        templates : array
            Template data in a 2D array containing the luminosity density
            (normalised over the full spectrum) with α in the first axis and
            λ in the second.

        """

        self.alpha_grid = alpha_grid
        self.wavelength_grid = wavelength_grid
        self.templates = templates

    def get_template(self, alpha):
        """Get the IR template corresponding to the given alpha

        The template for the given α is interpolated from the Dale and Helou
        (2002) template set. The new template is normalised in case of errors
        in the interpolation.

        Parameters
        ----------
        alpha : float
            α slope of the IR.

        Returns
        -------
        luminosity : array
            The luminosity density vector base on the template set λ grid.

        """
        result = interpolate.interp1d(self.alpha_grid, self.templates,
                                      axis=0)(alpha)
        # Return the normalised result
        return result / np.trapz(result, x=self.wavelength_grid)
