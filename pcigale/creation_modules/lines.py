# -*- coding: utf-8 -*-
# Copyright (C) 2014 University of Cambridge
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Médéric Boquien <mboquien@ast.cam.ac.uk>

from collections import OrderedDict
import numpy as np
from pcigale.data import Database
import scipy.constants as cst
from . import CreationModule


class Module(CreationModule):
    """
    Module computing the nebular emission lines from the ultraviolet
    to the near-infrared

    Given the number of Lyman continuum photons, we compute the Hβ line
    luminosity. We then compute the other lines using the
    metallicity-dependent templates that provide the ratio between individual
    lines and Hβ.

    Information added to the SED: NAME_escape_fraction.

    """

    parameter_list = OrderedDict([
        ('logU', (
            'float',
            "Radiation field intensity",
            -2.
        )),
        ('escape_fraction', (
            'float',
            "Fraction of Lyman continuum photons escaping the galaxy",
            0.
        )),
        ('f_dust', (
            'float',
            "Fraction of Lyman continuum photons absorbed by dust",
            0.
        )),
        ('lines_width', (
            'float',
            "Line width in km s¯¹",
            300.
        ))
    ])

    out_parameter_list = OrderedDict([
        ('logU', "Radiation field intensity"),
        ('escape_fraction', "Fraction of Lyman continuum photons escaping "
         "the galaxy"),
        ('f_dust', "Fraction of Lyman continuum photons absorbed by dust"),
        ('lines_width', "Width of the lines")
    ])

    def _init_code(self):
        """Get the nebular emission lines out of the database and resample
           them to see the line profile.
        """

        with Database() as database:
            self.lines_template = {m:database.
                                   get_lines(m,self.parameters['logU'])
                                   for m in database.get_lines_metallicities()
                                  }
        
        lines_width = self.parameters['lines_width'] * 1e3
        for lines in self.lines_template.values():
            new_wave = np.array([])
            for line_wave in lines.wave:
                width = line_wave * lines_width / cst.c
                new_wave = np.concatenate((new_wave, np.linspace(line_wave -
                                                                 3. * width,
                                                                 line_wave
                                                                 + 3. * 
                                                                 width,
                                                       19)))
            new_wave = np.sort(new_wave)
            new_flux = np.zeros_like(new_wave)
            for line_flux, line_wave in zip(lines.ratio, lines.wave):
                width = line_wave * lines_width / cst.c
                new_flux += (line_flux * np.exp(- 4. * np.log(2.) *
                            (new_wave - line_wave) ** 2. /  (width * width)) /
                            (width * np.sqrt(np.pi / np.log(2.)) / 2.))
            lines.wave = new_wave
            lines.ratio = new_flux

        # We compute the conversion coefficient to determine the fluxes using
        # the formula of Krüger+95: h×ν(Hβ)*α(Hβ,10000K)/αB(10000K)
        # Osterbrock 1989 gives:
        # α(Hβ,10000K)=3.03×10¯¹⁴
        # αB(10000K)2.59×10¯¹³
        self.conv = (13.598 * cst.electron_volt * (1. / 4. - 1. / 16.) *
                    3.03e-14 / 2.59e-13 *
                    (1. - self.parameters['escape_fraction'] -
                     self.parameters['f_dust'])
                    )


    def process(self, sed):
        """Add the nebular emission lines

        Parameters
        ----------
        sed  : pcigale.sed.SED object
        parameters : dictionary containing the parameters

        """

        # Base name for adding information to the SED.
        name = self.name or 'lines'

        escape_fraction = self.parameters['escape_fraction']
        f_dust = self.parameters['f_dust']
        NLy_old = sed.info['ssp_n_ly_old']
        NLy_young = sed.info['ssp_n_ly_young']
        lines = self.lines_template[sed.info['ssp_metallicity']]

        sed.add_module(name, self.parameters)
        sed.add_info(name + '_escape_fraction', escape_fraction)
        sed.add_info(name + '_f_dust', f_dust)

        sed.add_contribution('lines_old', lines.wave, NLy_old * lines.ratio *
                             self.conv)
        sed.add_contribution('lines_young', lines.wave, NLy_young * lines.ratio *
                             self.conv)
