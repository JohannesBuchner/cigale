# -*- coding: utf-8 -*-
# Copyright (C) 2016, Universidad de Antofagasta
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Médéric Boquien

from collections import OrderedDict

from astropy.table import Table, Column

def save_sed_to_fits(sed, prefix, norm=1.):
    """
    Save a SED object to fits files

    Parameters
    ----------
    sed: a pcigale.sed.SED object
        The SED to save
    prefix: string
        Prefix of the fits file containing the path and the id of the model
    norm: float
        Normalisation factor of the SED

    """
    info = OrderedDict()
    for name in sed.mass_proportional_info:
        info[name] = str(norm * sed.info[name])
    else:
        info[name] = str(sed.info[name])

    table = Table(meta=info)
    table['wavelength'] = Column(sed.wavelength_grid, unit="nm")
    table['Fnu'] = Column(norm * sed.fnu, unit="mJy")
    table['L_lambda_total'] = Column(norm * sed.luminosity, unit="W/nm")
    for name in sed.contribution_names:
        table[name] = Column(norm * sed.get_lumin_contribution(name),
        unit="W/nm")
    table.write("{}_best_model.fits".format(prefix))

    table = Table(meta=info)
    table["time"] = Column(sed.sfh[0], unit="Myr")
    table["SFR"] = Column(norm * sed.sfh[1], unit="Msun/yr")
    table.write("{}_SFH.fits".format(prefix))
