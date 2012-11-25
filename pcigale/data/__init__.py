# -*- coding: utf-8 -*-
"""
Copyright (C) 2012 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

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
from sqlalchemy.orm import sessionmaker
from .filters import Filter
from .ssp_m2005 import SspM2005
from .ir_templates_dh2002 import IrTemplatesDH2002


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
    age_grid = Column(PickleType)
    wavelength_grid = Column(PickleType)
    mass_table = Column(PickleType)
    spec_table = Column(PickleType)

    def __init__(self, ssp):
        self.imf = ssp.imf
        self.metallicity = ssp.metallicity
        self.age_grid = ssp.age_grid
        self.wavelength_grid = ssp.wavelength_grid
        self.mass_table = ssp.mass_table
        self.spec_table = ssp.spec_table


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


class Database(object):
    """Object giving access to pcigale database."""

    def __init__(self, writable=False):
        """
        Create a collection giving access to access the pcigale database.

        Parametres
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

        Parametres
        ----------
        pcigale_filter : pcigale.data.Filter
        """
        if self.is_writable:
            self.session.add(_Filter(pcigale_filter))
            try:
                self.session.commit()
            except exc.IntegrityError:
                self.session.rollback()
                raise StandardError('The filter is yet in the base.')
        else:
            raise StandardError('The database is not writable.')

    def add_ssp_m2005(self, ssp_m2005):
        """
        Add a Maraston 2005 SSP to pcigale database

        Parametres
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
                raise StandardError('The SSP is yet in the base.')
        else:
            raise StandardError('The database is not writable.')

    def add_dh2002_infrared_templates(self, data):
        """
        Add Dale and Helou (2002) templates the collection.

        Parametres
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
                raise StandardError('The template is yet in the base.')
        else:
            raise StandardError('The database is not writable.')

    def get_filter(self, name):
        """
        Get a specific filter from the collection

        Parametres
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

    def get_ssp_m2005(self, imf, metallicity):
        """
        Query the database for a Maraston 2005 SSP corresponding to the given
        initial mass function and metallicity.

        Parametres
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
            return SspM2005(result.imf, result.metallicity, result.age_grid,
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

    def get_filter_list(self):
        """Get the list of the filters in the database.

        Returns
        -------
        names, lamdba_eff : array, dictionnary
            names is the list of the filter names and lambda_eff is a
            dictionnary associating the effective wavelength (in nm) to the
            filter name
        """
        result = self.session.query(_Filter.name,
                                    _Filter.effective_wavelength).all()
        result = dict(result)
        return result.keys(), result

    def parse_filters(self):
        """Generator to parse the filter database."""
        for filt in self.session.query(Filter):
            yield Filter(filt.name, filt.description, filt.trans_type,
                         filt.trans_table, filt.effective_wavelength)

    def parse_ssp_m2005(self):
        """Generator to parse the Maraston 2005 SSP database."""
        for ssp in self.session.query(SspM2005):
            yield SspM2005(ssp.imf, ssp.metallicity, ssp.age_grid,
                           ssp.wavelength_grid, ssp.mass_table,
                           ssp.spec_table)
