# -*- coding: utf-8 -*-
# Copyright (C) 2014 University of Cambridge
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Médéric Boquien <mboquien@ast.cam.ac.uk>

import numpy as np
from pcigale.data import Database
import scipy.constants as cst
from . import CreationModule


class NebularEmission(CreationModule):
    """
    Module computing the nebular emission from the ultraviolet to the
    near-infrared. It includes both the nebular lines and the nubular
    continuum. It takes into account the escape fraction and the absorption by
    dust.

    Given the number of Lyman continuum photons, we compute the Hβ line
    luminosity. We then compute the other lines using the
    metallicity-dependent templates that provide the ratio between individual
    lines and Hβ. The nebular continuum is scaled directly from the number of
    ionizing photons.

    """

    parameter_list = dict([
        ('logU', (
            'float',
            "Ionisation parameter",
            -2.
        )),
        ('f_esc', (
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
            "Line width in km/s",
            300.
        ))
    ])

    out_parameter_list = dict([
        ('logU', "Ionisation parameter"),
        ('f_esc', "Fraction of Lyman continuum photons escaping "
         "the galaxy"),
        ('f_dust', "Fraction of Lyman continuum photons absorbed by dust"),
        ('lines_width', "Width of the nebular lines")
    ])

    def _init_code(self):
        """Get the nebular emission lines out of the database and resample
           them to see the line profile. Compute scaling coefficients.
        """

        fesc = self.parameters['f_esc']
        fdust = self.parameters['f_dust']

        if fesc < 0. or fesc > 1:
            raise Exception("Escape fraction must be between 0 and 1")

        if fdust < 0 or fdust > 1:
            raise Exception("Fraction of lyman photons absorbed by dust must "
                            "be between 0 and 1")

        if fesc + fdust > 1:
            raise Exception("Escape fraction+f_dust>1")

        with Database() as db:
            self.lines_template = {m:
                                   db.
                                   get_nebular_lines(m,
                                                     self.parameters['logU'])
                                   for m in db.get_nebular_lines_parameters()
                                   ['metallicity']
                                   }
            self.cont_template = {m:
                                  db.
                                  get_nebular_continuum(m,
                                                        self.parameters['logU'])
                                  for m in db.get_nebular_continuum_parameters()
                                  ['metallicity']
                                  }

        lines_width = self.parameters['lines_width'] * 1e3
        for lines in self.lines_template.values():
            new_wave = np.array([])
            for line_wave in lines.wave:
                width = line_wave * lines_width / cst.c
                new_wave = np.concatenate((new_wave,
                                           np.linspace(line_wave - 3. * width,
                                                       line_wave + 3. * width,
                                                       9)))
            new_wave.sort()
            new_flux = np.zeros_like(new_wave)
            for line_flux, line_wave in zip(lines.ratio, lines.wave):
                width = line_wave * lines_width / cst.c
                new_flux += (line_flux * np.exp(- 4. * np.log(2.) *
                             (new_wave - line_wave) ** 2. / (width * width)) /
                             (width * np.sqrt(np.pi / np.log(2.)) / 2.))
            lines.wave = new_wave
            lines.ratio = new_flux

        # We compute the conversion coefficient to determine the fluxes using
        # the formula of Inoue 2011: LHβ=Q(H)*γHβ(10000K)/αβ(10000K)
        gamma_Hbeta = 1.23e-38  # Inoue 2011, W m³
        alpha_B = 2.58e-19  # Ferland 1980, m³ s¯¹

        # To take into acount the escape fraction and the fraction of Lyman
        # continuum photons absorbed by dust we correct by a factor
        # k=(1-fesc-fdust)/(1+(α1/αβ)*(fesc+fdust))
        alpha_1 = 1.54e-19  # αA-αB, Ferland 1980, m³ s¯¹
        k = (1. - fesc - fdust) / (1. + alpha_1 / alpha_B * (fesc + fdust))

        self.conv_line = gamma_Hbeta / alpha_B * k
        self.conv_cont = k

    def process(self, sed):
        """Add the nebular emission lines

        Parameters
        ----------
        sed: pcigale.sed.SED object
        parameters: dictionary containing the parameters

        """
        NLy_old = sed.info['stellar.n_ly_old']
        NLy_young = sed.info['stellar.n_ly_young']
        lines = self.lines_template[sed.info['stellar.metallicity']]
        cont = self.cont_template[sed.info['stellar.metallicity']]

        sed.add_module(self.name, self.parameters)
        sed.add_info('nebular.logU', self.parameters['logU'])
        sed.add_info('nebular.f_esc', self.parameters['f_esc'])
        sed.add_info('nebular.f_dust', self.parameters['f_dust'])
        sed.add_info('nebular.lines_width', self.parameters['lines_width'])
        sed.add_info('dust.luminosity', (sed.info['stellar.lum_ly_young'] +
                     sed.info['stellar.lum_ly_old']) *
                     self.parameters['f_dust'], True)

        sed.add_contribution('nebular.lines_old', lines.wave,
                             lines.ratio * NLy_old * self.conv_line)
        sed.add_contribution('nebular.lines_young', lines.wave,
                             lines.ratio * NLy_young * self.conv_line)

        sed.add_contribution('nebular.continuum_old', cont.wave,
                             cont.lumin * NLy_old * self.conv_cont)
        sed.add_contribution('nebular.continuum_young', cont.wave,
                             cont.lumin * NLy_young * self.conv_cont)

# CreationModule to be returned by get_module
Module = NebularEmission
