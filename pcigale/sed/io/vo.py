# -*- coding: utf-8 -*-
# Copyright (C) 2014 Yannick Roehlly
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

from astropy.io.votable.tree import VOTableFile, Resource, Table, Field, Info


def save_sed_to_vo(sed, filename, mass=1.):
    """
    Save a SED object to a VO-Table file

    Parameters
    ----------
    sed : a pcigale.sed.SED object
        The SED to save.
    filename : string
        Name of the file to save the SED to.
    mass : float
        Mass of the galaxy in solar mass.

    """
    votable = VOTableFile()

    spectra_resource = Resource(id="Spectra")
    votable.resources.append(spectra_resource)

    # Total F_nu
    fnu_table = Table(votable, name="Fnu", id="Fnu")
    spectra_resource.tables.append(fnu_table)
    fnu_table.fields.extend([
        Field(votable, name="wavelength", datatype="double", unit="nm",
              ucd="em.wl"),
        Field(votable, name="F_nu", datatype="double", unit="mJy",
              ucd="phot.flux")
    ])
    fnu_table.create_arrays(len(sed.wavelength_grid))
    fnu_table.array["wavelength"] = sed.wavelength_grid
    fnu_table.array["F_nu"] = mass * sed.fnu

    # F_lambda contributions and total
    flambda_table = Table(votable, name="Flambda", id="Flambda")
    spectra_resource.tables.append(flambda_table)
    flambda_fields = [
        Field(votable, name="wavelength", datatype="double", unit="nm",
              ucd="em.wl"),
        Field(votable, name="F_lambda_total", datatype="double", unit="W/nm",
              ucd="phot.flux")]
    for name in sed.contribution_names:
        flambda_fields.append(Field(votable, name=name, datatype="double",
                                    unit="W/nm", ucd="phot.flux"))
    flambda_table.fields.extend(flambda_fields)
    flambda_table.create_arrays(len(sed.wavelength_grid))
    flambda_table.array["wavelength"] = sed.wavelength_grid
    flambda_table.array["F_lambda_total"] = mass * sed.luminosity
    for name in sed.contribution_names:
        flambda_table.array[name] = mass * sed.get_lumin_contribution(name)

    # SFH
    sfh_resource = Resource(id="Star_Formation_History")
    votable.resources.append(sfh_resource)
    sfh_table = Table(votable, name="SFH", id="SFH")
    sfh_resource.tables.append(sfh_table)
    sfh_table.fields.extend([
        Field(votable, name="time", datatype="double", unit="Myr",
              ucd="time.age"),
        Field(votable, name="SFR", datatype="double", unit="Msun/yr",
              ucd="phys.SFR")
    ])
    sfh_table.create_arrays(len(sed.sfh[0]))
    sfh_table.array["time"] = sed.sfh[0]
    sfh_table.array["SFR"] = mass * sed.sfh[1]

    # SED information to keywords
    votable.infos.append(Info(name="Galaxy mass in Msun", value=mass))
    votable.infos.append(Info(name="Redshift", value=sed.redshift))
    for name, value in sed.info.items():
        if name in sed.mass_proportional_info:
            votable.infos.append(Info(name=name, value=mass * value))
        else:
            votable.infos.append(Info(name=name, value=value))

    votable.set_all_tables_format('binary')
    votable.to_xml(filename)
