# -*- coding: utf-8 -*-
# Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly <yannick.roehlly@oamp.fr>

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
                        Float, Integer, PickleType)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from .filters import Filter
from .ssp_m2005 import SspM2005
from .ssp_bc03 import SspBC03
from .ir_templates_dh2002 import IrTemplatesDH2002
from .ir_agn_templates_dale2014 import DALE2014
from .ir_models_dl2007 import DL2007
from .agn_fritz2006 import AgnFritz2006


DATABASE_FILE = pkg_resources.resource_filename(__name__, 'data.db')

ENGINE = create_engine('sqlite:///' + DATABASE_FILE, echo=False)
BASE = declarative_base()
SESSION = sessionmaker(bind=ENGINE)


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


class _SspM2005(BASE):
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


class _SspBC03(BASE):
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


class _DH2002InfraredTemplates(BASE):
    """Storage for Dale and Helou (2002) infra-red templates

    The Dale and Helou (2002) template are gathered in a unique numpy array,
    nevertheless, they are stored in their own table with a unique row.

    """

    __tablename__ = 'dh2002_templates'

    name = Column(String, primary_key=True)
    description = Column(Text)
    data = Column(PickleType)

    def __init__(self, name, description, data):
        self.name = name
        self.description = description
        self.data = data

class _DALE2014(BASE):
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


class _Fritz2006AGN(BASE):
    """Storage for Fritz et al. (2006) models
    """

    __tablename__ = 'fritz2006_agn'
    model_nb = Column(Integer, primary_key=True)
    agn_type = Column(Integer)
    r_ratio = Column(Float)
    tau = Column(Float)
    beta = Column(Float)
    gamma = Column(Float)
    theta = Column(Float)
    psy = Column(Float)
    wave = Column(PickleType)
    luminosity = Column(PickleType)

    def __init__(self, agn):
        self.model_nb = agn.model_nb
        self.agn_type = agn.agn_type
        self.r_ratio = agn.r_ratio
        self.tau = agn.tau
        self.beta = agn.beta
        self.gamma = agn.gamma
        self.theta = agn.theta
        self.psy = agn.psy
        self.wave = agn.wave
        self.luminosity = agn.luminosity


class Database(object):
    """Object giving access to pcigale database."""

    def __init__(self, writable=False):
        """
        Create a collection giving access to access the pcigale database.

        Parameters
        ----------
        writable : boolean
            If True the user will be able to write new data in the database
            (but he/she must have a writable access to the sqlite file). By
            default, False.
        """
        self.session = SESSION()
        self.is_writable = writable

    def upgrade_base(self):
        """ Upgrade the table schemas in the database
        """
        if self.is_writable:
            BASE.metadata.create_all(ENGINE)
        else:
            raise StandardError('The database is not writable.')

    def close(self):
        """ Close the connection to the database

        TODO: It would be better to wrap the database use inside a context
        manager.
        """
        self.session.close_all()

    def add_filter(self, pcigale_filter):
        """
        Add a filter to pcigale database.

        Parameters
        ----------
        pcigale_filter : pcigale.data.Filter
        """
        if self.is_writable:
            self.session.add(_Filter(pcigale_filter))
            try:
                self.session.commit()
            except exc.IntegrityError:
                self.session.rollback()
                raise StandardError('The filter is already in the base.')
        else:
            raise StandardError('The database is not writable.')

    def add_ssp_m2005(self, ssp_m2005):
        """
        Add a Maraston 2005 SSP to pcigale database

        Parameters
        ----------
        ssp : pcigale.base.SspM2005

        """
        if self.is_writable:
            ssp = _SspM2005(ssp_m2005)
            self.session.add(ssp)
            try:
                self.session.commit()
            except exc.IntegrityError:
                self.session.rollback()
                raise StandardError('The SSP is already in the base.')
        else:
            raise StandardError('The database is not writable.')

    def add_ssp_bc03(self, ssp_bc03):
        """
        Add a Bruzual and Charlot 2003 SSP to pcigale database

        Parameters
        ----------
        ssp : pcigale.data.SspBC03

        """
        if self.is_writable:
            ssp = _SspBC03(ssp_bc03)
            self.session.add(ssp)
            try:
                self.session.commit()
            except exc.IntegrityError:
                self.session.rollback()
                raise StandardError('The SSP is already in the base.')
        else:
            raise StandardError('The database is not writable.')

    def add_dh2002_infrared_templates(self, data):
        """
        Add Dale and Helou (2002) templates the collection.

        Parameters
        ----------
        data : array
            Array containing the templates data.

        """
        name = 'dh2002'

        description = ("These are the Dale & Helou (2002) infra-red "
                       "templates to which the stellar emission was "
                       "subtracted (Nohl et al., 2009). The data is a "
                       "tuple composed of the alpha grid, the lambda grid "
                       "and a 2D array of luminosity density (normalised "
                       "over the full spectrum) with the alpha in the first "
                       "axis and the lambda in the second.")

        if self.is_writable:
            template = _DH2002InfraredTemplates(name, description, data)
            self.session.add(template)
            try:
                self.session.commit()
            except exc.IntegrityError:
                self.session.rollback()
                raise StandardError('The template is already in the base.')
        else:
            raise StandardError('The database is not writable.')
            
        
    def add_dale2014(self, iragn):
        """
        Add Dale et al (2014) templates the collection.

        Parameters
        ----------
        iragn : pcigale.data.DALE2014

        """

        if self.is_writable:
            template = _DALE2014(iragn)
            self.session.add(template)
            try:
                self.session.commit()
            except exc.IntegrityError:
                self.session.rollback()
                raise StandardError('The DALE2014 template is already in the base.')
        else:
            raise StandardError('The database is not writable.')

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
                raise StandardError('The DL07 model is already in the base.')
        else:
            raise StandardError('The database is not writable.')

    def add_fritz2006_agn(self, agn):
        """
        Add a Fritz et al. (2006) AGN model to the database.

        Parameters
        ----------
        agn : pcigale.data.AgnFritz2006

        """
        if self.is_writable:
            self.session.add(_Fritz2006AGN(agn))
            try:
                self.session.commit()
            except exc.IntegrityError:
                self.session.rollback()
                raise StandardError('The agn model is already in the base.')
        else:
            raise StandardError('The database is not writable.')

    def get_filter(self, name):
        """
        Get a specific filter from the collection

        Parameters
        ----------
        name : string
            Name of the filter

        Returns
        -------
        filter : pcigale.base.Filter
            The Filter object. If the filter is not in the database,
            returns None.
        """
        result = (self.session.query(_Filter).
                  filter(_Filter.name == name).
                  first())
        if result:
            return Filter(result.name, result.description,
                          result.trans_type, result.trans_table,
                          result.effective_wavelength)
        else:
            return None

    def get_ssp_bc03(self, imf, metallicity):
        """
        Query the database for the Bruzual and Charlot 2003 SSP corresponding
        to the given initial mass function and metallicity.

        Parameters
        ----------
        imf : string
            Initial mass function (salp for Salpeter, chab for Chabrier)
        metallicity : float
            0.02 for Solar metallicity
        Returns
        -------
        ssp : pcigale.data.SspBC03
            The SspBC03 object. If no SSP corresponds to the given imf and
            metallicity, returns None.

        """
        result = self.session.query(_SspBC03)\
            .filter(_SspBC03.imf == imf)\
            .filter(_SspBC03.metallicity == metallicity)\
            .first()
        if result:
            return SspBC03(result.imf, result.metallicity, result.time_grid,
                           result.wavelength_grid, result.color_table,
                           result.lumin_table)
        else:
            return None

    def get_ssp_m2005(self, imf, metallicity):
        """
        Query the database for a Maraston 2005 SSP corresponding to the given
        initial mass function and metallicity.

        Parameters
        ----------
        imf : string
            Initial mass function (ss for Salpeter, kr for Kroupa)
        metallicity : float
            [Z/H] = Log10(Z/Zsun) - Log10(H/Hsun)

        Returns
        -------
        ssp : pcigale.base.SspM2005
            The SspM2005 object. If no SSP corresponds to the given imf and
            metallicity, returns None.

        """
        result = self.session.query(_SspM2005)\
            .filter(_SspM2005.imf == imf)\
            .filter(_SspM2005.metallicity == metallicity)\
            .first()
        if result:
            return SspM2005(result.imf, result.metallicity, result.time_grid,
                            result.wavelength_grid, result.mass_table,
                            result.spec_table)
        else:
            return None

    def get_dh2002_infrared_templates(self):
        """
        Get the Dale and Helou infrared templates from the database

        Returns
        -------
        template : pcigale.base.IrTemplatesDH2002
            The Dale and Helou (2002) infrared templates. If they are not in
            the database, return None.

        """
        result = (self.session.query(_DH2002InfraredTemplates).
                  filter(_DH2002InfraredTemplates.name == 'dh2002').
                  first())
        if result:
            return IrTemplatesDH2002(result.data[0], result.data[1],
                                     result.data[2])
        else:
            return None

    def get_dale2014(self, fracAGN, alpha):
        """
        Get the Dale et al (2014) template corresponding to the given set of
        parameters.

        Parameters
        ----------
        fracAGN: float
            contribution of the AGN to the IR luminosity
        alpha: float
            alpha corresponding to the updated Dale & Helou (2002) star forming template.

        """
        result = (self.session.query(_DALE2014).
                  filter(_DALE2014.fracAGN == fracAGN).
                  filter(_DALE2014.alpha == alpha).
                  first())
        if result:
            return DALE2014(result.fracAGN, result.alpha, result.wave,
                          result.lumin)
        else:
            return None

    def get_agn_fritz2006(self, model_nb):
        """
        Get the Fritz et al. (2006) AGN model corresponding to the number.

        Parameters
        ----------
        model_nb : integer
            Model number.

        Returns
        -------
        agn : pcigale.data.AgnFritz2006
            The AGN model or not if no model exists for the given number.

        """
        result = (self.session.query(_Fritz2006AGN).
                  filter(_Fritz2006AGN.model_nb == model_nb).
                  first())
        if result:
            return AgnFritz2006(result.model_nb, result.agn_type,
                                result.r_ratio, result.tau, result.beta,
                                result.gamma, result.theta, result.psy,
                                result.wave, result.luminosity)
        else:
            return None

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
        gamma: float
            Fraction of the dust exposed from Umin to Umax

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
            return None

    def get_filter_list(self):
        """Get the list of the filters in the database.

        Returns
        -------
        names, lambda_eff : array, dictionary
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

    def parse_ssp_m2005(self):
        """Generator to parse the Maraston 2005 SSP database."""
        for ssp in self.session.query(_SspM2005):
            yield SspM2005(ssp.imf, ssp.metallicity, ssp.time_grid,
                           ssp.wavelength_grid, ssp.mass_table,
                           ssp.spec_table)
