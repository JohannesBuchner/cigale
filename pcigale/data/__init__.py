# -*- coding: utf-8 -*-
# Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

"""
This is the database where we store some data used by pcigale:
 - the information relative to the filters
 - the single stellar populations as defined in Marason (2005)
 - the infra-red templates from Dale and Helou (2002)

The classes for these various objects are described in pcigale.data
sub-packages. The corresponding underscored classes here are used by the
SqlAlchemy ORM to store the data in a unique SQLite3 database.

"""

import pkg_resources
from sqlalchemy import (create_engine, exc, Column, String, Text,
                        Float, PickleType)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import class_mapper, sessionmaker
import numpy as np
from .filters import Filter
from .m2005 import M2005
from .bc03 import BC03
from .dale2014 import Dale2014
from .dl2007 import DL2007
from .dl2014 import DL2014
from .fritz2006 import Fritz2006
from .nebular_continuum import NebularContinuum
from .nebular_lines import NebularLines


DATABASE_FILE = pkg_resources.resource_filename(__name__, 'data.db')

ENGINE = create_engine('sqlite:///' + DATABASE_FILE, echo=False)
BASE = declarative_base()
SESSION = sessionmaker(bind=ENGINE)


class DatabaseLookupError(Exception):
    """
    A custom exception raised when a search in the database does not find a
    result.
    """


class DatabaseInsertError(Exception):
    """
    A custom exception raised when one tries to insert in the database
    something that is already in it.
    """


class _Filter(BASE):
    """ Storage for filters
    """

    __tablename__ = 'filters'

    name = Column(String, primary_key=True)
    description = Column(String)
    trans_type = Column(String)
    trans_table = Column(PickleType)
    effective_wavelength = Column(Float)

    def __init__(self, f):
        self.name = f.name
        self.description = f.description
        self.trans_type = f.trans_type
        self.trans_table = f.trans_table
        self.effective_wavelength = f.effective_wavelength


class _M2005(BASE):
    """Storage for Maraston 2005 SSP
    """

    __tablename__ = 'maraston2005'

    imf = Column(String, primary_key=True)
    metallicity = Column(Float, primary_key=True)
    time_grid = Column(PickleType)
    wavelength_grid = Column(PickleType)
    mass_table = Column(PickleType)
    spec_table = Column(PickleType)

    def __init__(self, ssp):
        self.imf = ssp.imf
        self.metallicity = ssp.metallicity
        self.time_grid = ssp.time_grid
        self.wavelength_grid = ssp.wavelength_grid
        self.mass_table = ssp.mass_table
        self.spec_table = ssp.spec_table


class _BC03(BASE):
    """Storage for Bruzual and Charlot 2003 SSP
    """

    __tablename__ = "bc03"

    imf = Column(String, primary_key=True)
    metallicity = Column(Float, primary_key=True)
    time_grid = Column(PickleType)
    wavelength_grid = Column(PickleType)
    color_table = Column(PickleType)
    lumin_table = Column(PickleType)

    def __init__(self, ssp):
        self.imf = ssp.imf
        self.metallicity = ssp.metallicity
        self.time_grid = ssp.time_grid
        self.wavelength_grid = ssp.wavelength_grid
        self.color_table = ssp.color_table
        self.lumin_table = ssp.lumin_table


class _Dale2014(BASE):
    """Storage for Dale et al (2014) infra-red templates
    """

    __tablename__ = 'dale2014_templates'
    fracAGN = Column(Float, primary_key=True)
    alpha = Column(String, primary_key=True)
    wave = Column(PickleType)
    lumin = Column(PickleType)

    def __init__(self, iragn):
        self.fracAGN = iragn.fracAGN
        self.alpha = iragn.alpha
        self.wave = iragn.wave
        self.lumin = iragn.lumin


class _DL2007(BASE):
    """Storage for Draine and Li (2007) IR models
    """

    __tablename__ = 'DL2007_models'
    qpah = Column(Float, primary_key=True)
    umin = Column(Float, primary_key=True)
    umax = Column(Float, primary_key=True)
    wave = Column(PickleType)
    lumin = Column(PickleType)

    def __init__(self, model):
        self.qpah = model.qpah
        self.umin = model.umin
        self.umax = model.umax
        self.wave = model.wave
        self.lumin = model.lumin


class _DL2014(BASE):
    """Storage for the updated Draine and Li (2007) IR models
    """

    __tablename__ = 'DL2014_models'
    qpah = Column(Float, primary_key=True)
    umin = Column(Float, primary_key=True)
    umax = Column(Float, primary_key=True)
    alpha = Column(Float, primary_key=True)
    wave = Column(PickleType)
    lumin = Column(PickleType)

    def __init__(self, model):
        self.qpah = model.qpah
        self.umin = model.umin
        self.umax = model.umax
        self.alpha = model.alpha
        self.wave = model.wave
        self.lumin = model.lumin


class _Fritz2006(BASE):
    """Storage for Fritz et al. (2006) models
    """

    __tablename__ = 'fritz2006'
    r_ratio = Column(Float, primary_key=True)
    tau = Column(Float, primary_key=True)
    beta = Column(Float, primary_key=True)
    gamma = Column(Float, primary_key=True)
    opening_angle = Column(Float, primary_key=True)
    psy = Column(Float, primary_key=True)
    wave = Column(PickleType)
    lumin_therm = Column(PickleType)
    lumin_scatt = Column(PickleType)
    lumin_agn = Column(PickleType)

    def __init__(self, agn):
        self.r_ratio = agn.r_ratio
        self.tau = agn.tau
        self.beta = agn.beta
        self.gamma = agn.gamma
        self.opening_angle = agn.opening_angle
        self.psy = agn.psy
        self.wave = agn.wave
        self.lumin_therm = agn.lumin_therm
        self.lumin_scatt = agn.lumin_scatt
        self.lumin_agn = agn.lumin_agn


class _NebularLines(BASE):
    """Storage for line templates
    """

    __tablename__ = 'nebular_lines'
    metallicity = Column(Float, primary_key=True)
    logU = Column(Float, primary_key=True)
    wave = Column(PickleType)
    ratio = Column(PickleType)

    def __init__(self, nebular_lines):
        self.metallicity = nebular_lines.metallicity
        self.logU = nebular_lines.logU
        self.wave = nebular_lines.wave
        self.ratio = nebular_lines.ratio


class _NebularContinuum(BASE):
    """Storage for nebular continuum templates
    """

    __tablename__ = 'nebular_continuum'
    metallicity = Column(Float, primary_key=True)
    logU = Column(Float, primary_key=True)
    wave = Column(PickleType)
    lumin = Column(PickleType)

    def __init__(self, nebular_continuum):
        self.metallicity = nebular_continuum.metallicity
        self.logU = nebular_continuum.logU
        self.wave = nebular_continuum.wave
        self.lumin = nebular_continuum.lumin


class Database(object):
    """Object giving access to pcigale database."""

    def __init__(self, writable=False):
        """
        Create a collection giving access to access the pcigale database.

        Parameters
        ----------
        writable: boolean
            If True the user will be able to write new data in the database
            (but he/she must have a writable access to the sqlite file). By
            default, False.
        """
        self.session = SESSION()
        self.is_writable = writable

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def upgrade_base(self):
        """ Upgrade the table schemas in the database
        """
        if self.is_writable:
            BASE.metadata.create_all(ENGINE)
        else:
            raise Exception('The database is not writable.')

    def close(self):
        """ Close the connection to the database

        TODO: It would be better to wrap the database use inside a context
        manager.
        """
        self.session.close_all()

    def add_m2005(self, ssp_m2005):
        """
        Add a Maraston 2005 SSP to pcigale database

        Parameters
        ----------
        ssp: pcigale.base.M2005

        """
        if self.is_writable:
            ssp = _M2005(ssp_m2005)
            self.session.add(ssp)
            try:
                self.session.commit()
            except exc.IntegrityError:
                self.session.rollback()
                raise DatabaseInsertError('The SSP is already in the base.')
        else:
            raise Exception('The database is not writable.')

    def get_m2005(self, imf, metallicity):
        """
        Query the database for a Maraston 2005 SSP corresponding to the given
        initial mass function and metallicity.

        Parameters
        ----------
        imf: string
            Initial mass function (ss for Salpeter, kr for Kroupa)
        metallicity: float
            [Z/H] = Log10(Z/Zsun) - Log10(H/Hsun)

        Returns
        -------
        ssp: pcigale.base.M2005
            The M2005 object.

        Raises
        ------
        DatabaseLookupError: if the requested SSP is not in the database.

        """
        result = self.session.query(_M2005)\
            .filter(_M2005.imf == imf)\
            .filter(_M2005.metallicity == metallicity)\
            .first()
        if result:
            return M2005(result.imf, result.metallicity, result.time_grid,
                         result.wavelength_grid, result.mass_table,
                         result.spec_table)
        else:
            raise DatabaseLookupError(
                "The M2005 SSP for imf <{0}> and metallicity <{1}> is not in "
                "the database.".format(imf, metallicity))

    def get_m2005_parameters(self):
        """Get parameters for the Maraston 2005 stellar models.

        Returns
        -------
        paramaters: dictionary
            dictionary of parameters and their values
        """
        return self._get_parameters(_M2005)

    def add_bc03(self, ssp_bc03):
        """
        Add a Bruzual and Charlot 2003 SSP to pcigale database

        Parameters
        ----------
        ssp: pcigale.data.SspBC03

        """
        if self.is_writable:
            ssp = _BC03(ssp_bc03)
            self.session.add(ssp)
            try:
                self.session.commit()
            except exc.IntegrityError:
                self.session.rollback()
                raise DatabaseInsertError('The SSP is already in the base.')
        else:
            raise Exception('The database is not writable.')

    def get_bc03(self, imf, metallicity):
        """
        Query the database for the Bruzual and Charlot 2003 SSP corresponding
        to the given initial mass function and metallicity.

        Parameters
        ----------
        imf: string
            Initial mass function (salp for Salpeter, chab for Chabrier)
        metallicity: float
            0.02 for Solar metallicity
        Returns
        -------
        ssp: pcigale.data.BC03
            The BC03 object.

        Raises
        ------
        DatabaseLookupError: if the requested SSP is not in the database.

        """
        result = self.session.query(_BC03)\
            .filter(_BC03.imf == imf)\
            .filter(_BC03.metallicity == metallicity)\
            .first()
        if result:
            return BC03(result.imf, result.metallicity, result.time_grid,
                        result.wavelength_grid, result.color_table,
                        result.lumin_table)
        else:
            raise DatabaseLookupError(
                "The BC03 SSP for imf <{0}> and metallicity <{1}> is not in "
                "the database.".format(imf, metallicity))

    def get_bc03_parameters(self):
        """Get parameters for the Bruzual & Charlot 2003 stellar models.

        Returns
        -------
        paramaters: dictionary
            dictionary of parameters and their values
        """
        return self._get_parameters(_BC03)

    def add_dl2007(self, model):
        """
        Add a Draine and Li (2007) model to the database.

        Parameters
        ----------
        model: pcigale.data.DL2007

        """
        if self.is_writable:
            self.session.add(_DL2007(model))
            try:
                self.session.commit()
            except exc.IntegrityError:
                self.session.rollback()
                raise DatabaseInsertError(
                    'The DL07 model is already in the base.')
        else:
            raise Exception('The database is not writable.')

    def get_dl2007(self, qpah, umin, umax):
        """
        Get the Draine and Li (2007) model corresponding to the given set of
        parameters.

        Parameters
        ----------
        qpah: float
            Mass fraction of PAH
        umin: float
            Minimum radiation field
        umax: float
            Maximum radiation field

        Returns
        -------
        model: pcigale.data.DL2007
            The Draine and Li (2007) model.

        Raises
        ------
        DatabaseLookupError: if the requested model is not in the database.

        """
        result = (self.session.query(_DL2007).
                  filter(_DL2007.qpah == qpah).
                  filter(_DL2007.umin == umin).
                  filter(_DL2007.umax == umax).
                  first())
        if result:
            return DL2007(result.qpah, result.umin, result.umax, result.wave,
                          result.lumin)
        else:
            raise DatabaseLookupError(
                "The DL2007 model for qpah <{0}>, umin <{1}>, and umax <{2}> "
                "is not in the database.".format(qpah, umin, umax))

    def get_dl2007_parameters(self):
        """Get parameters for the DL2007 models.

        Returns
        -------
        paramaters: dictionary
            dictionary of parameters and their values
        """
        return self._get_parameters(_DL2007)

    def add_dl2014(self, model):
        """
        Add an updated Draine and Li (2007) model to the database.

        Parameters
        ----------
        model: pcigale.data.DL2014

        """
        if self.is_writable:
            self.session.add(_DL2014(model))
            try:
                self.session.commit()
            except exc.IntegrityError:
                self.session.rollback()
                raise DatabaseInsertError(
                    'The updated DL07 model is already in the base.')
        else:
            raise Exception('The database is not writable.')

    def get_dl2014(self, qpah, umin, umax, alpha):
        """
        Get the Draine and Li (2007) model corresponding to the given set of
        parameters.

        Parameters
        ----------
        qpah: float
            Mass fraction of PAH
        umin: float
            Minimum radiation field
        umin: float
            Maximum radiation field
        alpha: float
            Powerlaw slope dU/dM∝U¯ᵅ

        Returns
        -------
        model: pcigale.data.DL2014
            The updated Draine and Li (2007) model.

        Raises
        ------
        DatabaseLookupError: if the requested model is not in the database.

        """
        result = (self.session.query(_DL2014).
                  filter(_DL2014.qpah == qpah).
                  filter(_DL2014.umin == umin).
                  filter(_DL2014.umax == umax).
                  filter(_DL2014.alpha == alpha).
                  first())
        if result:
            return DL2014(result.qpah, result.umin, result.umax, result.alpha,
                          result.wave, result.lumin)
        else:
            raise DatabaseLookupError(
                "The DL2014 model for qpah <{0}>, umin <{1}>, umax <{2}>, and "
                "alpha <{3}> is not in the database."
                .format(qpah, umin, umax, alpha))

    def get_dl2014_parameters(self):
        """Get parameters for the DL2014 models.

        Returns
        -------
        paramaters: dictionary
            dictionary of parameters and their values
        """
        return self._get_parameters(_DL2014)

    def add_dale2014(self, iragn):
        """
        Add Dale et al (2014) templates the collection.

        Parameters
        ----------
        iragn: pcigale.data.Dale2014

        """

        if self.is_writable:
            template = _Dale2014(iragn)
            self.session.add(template)
            try:
                self.session.commit()
            except exc.IntegrityError:
                self.session.rollback()
                raise DatabaseInsertError(
                    'The Dale2014 template is already in the base.')
        else:
            raise Exception('The database is not writable.')

    def get_dale2014(self, frac_agn, alpha):
        """
        Get the Dale et al (2014) template corresponding to the given set of
        parameters.

        Parameters
        ----------
        frac_agn: float
            contribution of the AGN to the IR luminosity
        alpha: float
            alpha corresponding to the updated Dale & Helou (2002) star
            forming template.

        Returns
        -------
        template: pcigale.data.Dale2014
            The Dale et al. (2014) IR template.

        Raises
        ------
        DatabaseLookupError: if the requested template is not in the database.

        """
        result = (self.session.query(_Dale2014).
                  filter(_Dale2014.fracAGN == frac_agn).
                  filter(_Dale2014.alpha == alpha).
                  first())
        if result:
            return Dale2014(result.fracAGN, result.alpha, result.wave,
                            result.lumin)
        else:
            raise DatabaseLookupError(
                "The Dale2014 template for frac_agn <{0}> and alpha <{1}> "
                "is not in the database.".format(frac_agn, alpha))

    def get_dale2014_parameters(self):
        """Get parameters for the Dale 2014 models.

        Returns
        -------
        paramaters: dictionary
            dictionary of parameters and their values
        """
        return self._get_parameters(_Dale2014)

    def add_fritz2006(self, agn):
        """
        Add a Fritz et al. (2006) AGN model to the database.

        Parameters
        ----------
        agn: pcigale.data.Fritz2006

        """
        if self.is_writable:
            self.session.add(_Fritz2006(agn))
            try:
                self.session.commit()
            except exc.IntegrityError:
                self.session.rollback()
                raise DatabaseInsertError(
                    'The agn model is already in the base.')
        else:
            raise Exception('The database is not writable.')

    def get_fritz2006(self, r_ratio, tau, beta, gamma, opening_angle, psy):
        """
        Get the Fritz et al. (2006) AGN model corresponding to the number.

        Parameters
        ----------
        r_ratio: float
            Ratio of the maximum and minimum radii of the dust torus.
        tau: float
            Tau at 9.7µm
        beta: float
            Beta
        gamma: float
            Gamma
        opening_angle: float
            Opening angle of the dust torus.
        psy: float
            Angle between AGN axis and line of sight.
        wave: array of float
            Wavelength grid in nm.
        lumin_therm: array of float
            Luminosity density of the dust torus at each wavelength in W/nm.
        lumin_scatt: array of float
            Luminosity density of the scattered emission at each wavelength
            in W/nm.
        lumin_agn: array of float
            Luminosity density of the central AGN at each wavelength in W/nm.


        Returns
        -------
        agn: pcigale.data.Fritz2006
            The AGN model.

        Raises
        ------
        DatabaseLookupError: if the requested template is not in the database.

        """
        result = (self.session.query(_Fritz2006).
                  filter(_Fritz2006.r_ratio == r_ratio).
                  filter(_Fritz2006.tau == tau).
                  filter(_Fritz2006.beta == beta).
                  filter(_Fritz2006.gamma == gamma).
                  filter(_Fritz2006.opening_angle == opening_angle).
                  filter(_Fritz2006.psy == psy).
                  first())
        if result:
            return Fritz2006(result.r_ratio, result.tau, result.beta,
                             result.gamma, result.opening_angle, result.psy,
                             result.wave, result.lumin_therm,
                             result.lumin_scatt, result.lumin_agn)
        else:
            raise DatabaseLookupError(
                "The Fritz2006 model is not in the database.")

    def get_fritz2006_parameters(self):
        """Get parameters for the Fritz 2006 AGN models.

        Returns
        -------
        paramaters: dictionary
            dictionary of parameters and their values
        """
        return self._get_parameters(_Fritz2006)

    def add_nebular_lines(self, nebular_lines):
        """
        Add ultraviolet and optical line templates to the database.
        """
        if self.is_writable:
            self.session.add(_NebularLines(nebular_lines))
            try:
                self.session.commit()
            except exc.IntegrityError:
                self.session.rollback()
                raise Exception('The line is already in the base')
        else:
            raise Exception('The database is not writable')

    def get_nebular_lines(self, metallicity, logU):
        """
        Get the line ratios corresponding to the given set of parameters.

        Parameters
        ----------
        metallicity: float
            Gas phase metallicity
        logU: float
            Ionisation parameter
        """
        result = (self.session.query(_NebularLines).
                  filter(_NebularLines.metallicity == metallicity).
                  filter(_NebularLines.logU == logU).
                  first())
        if result:
            return NebularLines(result.metallicity, result.logU, result.wave,
                                result.ratio)
        else:
            return None

    def get_nebular_lines_parameters(self):
        """Get parameters for the nebular lines.

        Returns
        -------
        paramaters: dictionary
            dictionary of parameters and their values
        """
        return self._get_parameters(_NebularLines)

    def add_nebular_continuum(self, nebular_continuum):
        """
        Add nebular continuum templates to the database.
        """
        if self.is_writable:
            self.session.add(_NebularContinuum(nebular_continuum))
            try:
                self.session.commit()
            except exc.IntegrityError:
                self.session.rollback()
                raise Exception('The continuum template is already in the '
                                'base')
        else:
            raise Exception('The database is not writable')

    def get_nebular_continuum(self, metallicity, logU):
        """
        Get the nebular continuum corresponding to the given set of parameters.

        Parameters
        ----------
        metallicity: float
            Gas phase metallicity
        logU: float
            Ionisation parameter
        """
        result = (self.session.query(_NebularContinuum).
                  filter(_NebularContinuum.metallicity == metallicity).
                  filter(_NebularContinuum.logU == logU).
                  first())
        if result:
            return NebularContinuum(result.metallicity, result.logU,
                                    result.wave, result.lumin)
        else:
            return None

    def get_nebular_continuum_parameters(self):
        """Get parameters for the nebular continuum.

        Returns
        -------
        paramaters: dictionary
            dictionary of parameters and their values
        """
        return self._get_parameters(_NebularContinuum)

    def _get_parameters(self, schema):
        """Generic function to get parameters from an arbitrary schema.

        Returns
        -------
        parameters: dictionary
            Dictionary of parameters and their values
        """

        return {k.name: np.sort(
                [v[0] for v in set(self.session.query(schema).values(k))])
                for k in class_mapper(schema).primary_key}

    def add_filter(self, pcigale_filter):
        """
        Add a filter to pcigale database.

        Parameters
        ----------
        pcigale_filter: pcigale.data.Filter
        """
        if self.is_writable:
            self.session.add(_Filter(pcigale_filter))
            try:
                self.session.commit()
            except exc.IntegrityError:
                self.session.rollback()
                raise DatabaseInsertError('The filter is already in the base.')
        else:
            raise Exception('The database is not writable.')

    def get_filter(self, name):
        """
        Get a specific filter from the collection

        Parameters
        ----------
        name: string
            Name of the filter

        Returns
        -------
        filter: pcigale.base.Filter
            The Filter object.

        Raises
        ------
        DatabaseLookupError: if the requested filter is not in the database.

        """
        result = (self.session.query(_Filter).
                  filter(_Filter.name == name).
                  first())
        if result:
            return Filter(result.name, result.description,
                          result.trans_type, result.trans_table,
                          result.effective_wavelength)
        else:
            raise DatabaseLookupError(
                "The filter <{0}> is not in the database".format(name))

    def get_filter_list(self):
        """Get the list of the filters in the database.

        Returns
        -------
        names, lambda_eff: array, dictionary
            names is the list of the filter names and lambda_eff is a
            dictionary associating the effective wavelength (in nm) to the
            filter name
        """
        result = self.session.query(_Filter.name,
                                    _Filter.effective_wavelength).all()
        result = dict(result)
        return result.keys(), result

    def parse_filters(self):
        """Generator to parse the filter database."""
        for filt in self.session.query(_Filter):
            yield Filter(filt.name, filt.description, filt.trans_type,
                         filt.trans_table, filt.effective_wavelength)

    def parse_m2005(self):
        """Generator to parse the Maraston 2005 SSP database."""
        for ssp in self.session.query(_M2005):
            yield M2005(ssp.imf, ssp.metallicity, ssp.time_grid,
                        ssp.wavelength_grid, ssp.mass_table, ssp.spec_table)
