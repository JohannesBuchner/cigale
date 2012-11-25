# -*- coding: utf-8 -*-
"""
Copyright (C) 2012 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

This script is used the build pcigale internal database containing:
- The various filter transmission tables;
- The Maraston 2005 single stellar population (SSP) data;
- The Dale and Helou 2002 infra-red templates.

"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../'))
import glob
import numpy as np
from scipy import interpolate
from pcigale.data import Database, Filter, SspM2005


filters_dir = os.path.join(os.path.dirname(__file__), 'filters/')
m2005_dir = os.path.join(os.path.dirname(__file__), 'maraston2005/')
dh2002_dir = os.path.join(os.path.dirname(__file__), 'dh2002/')


def build_base():
    base = Database(writable=True)
    base.upgrade_base()

    print('#' * 78)
    ########################################################################
    # Filter transmission tables insertion                                 #
    ########################################################################
    print("1- Importing filters...\n")
    for filter_file in glob.glob(filters_dir + '*.dat'):
        with open(filter_file, 'r') as filter_file_read:
            filter_name = filter_file_read.readline().strip('# \n\t')
            filter_type = filter_file_read.readline().strip('# \n\t')
            filter_description = filter_file_read.readline().strip('# \n\t')
        filter_table = np.genfromtxt(filter_file)
        # The table is transposed to have table[0] containing the wavelength
        # and table[1] containing the transmission.
        filter_table = filter_table.transpose()
        # We convert the wavelength from Å to nm.
        filter_table[0] *= 0.1

        print("Importing %s... (%s points)" % (filter_name,
                                               filter_table.shape[1]))

        new_filter = Filter(filter_name, filter_description,
                            filter_type, filter_table)

        # We normalise the filter and compute the effective wavelength.
        new_filter.normalise()

        base.add_filter(new_filter)
    print("\nDONE\n")
    print('#' * 78)

    ########################################################################
    # Maraston 2005 SSP insertion                                          #
    ########################################################################
    print("2- Importing Maraston 2005 SSP\n")

    # Age grid (1My to 13.7Gy with 1My step)
    age_grid = np.arange(1e-3, 13.701, 1e-3)

    # Transpose the table to have access to each value vector on the first
    # axis
    kroupa_mass = np.genfromtxt(m2005_dir + 'stellarmass.kroupa').transpose()
    salpeter_mass = \
        np.genfromtxt(m2005_dir + '/stellarmass.salpeter').transpose()

    for spec_file in glob.glob(m2005_dir + '*.rhb'):

        print("Importing %s..." % spec_file)

        spec_table = np.genfromtxt(spec_file).transpose()
        metallicity = spec_table[1, 0]

        if 'krz' in spec_file:
            imf = 'kr'
            mass_table = np.copy(kroupa_mass)
        elif 'ssz' in spec_file:
            imf = 'ss'
            mass_table = np.copy(salpeter_mass)
        else:
            raise ValueError('Unknown IMF!!!')

        # Keep only the actual metallicity values in the mass table
        # we don't take the first column which contains metallicity
        mass_table = mass_table[1:, mass_table[0] == metallicity]

        # Interpolate the mass table over the new age grid
        mass_table = interpolate.interp1d(mass_table[0], mass_table)(age_grid)

        # Remove the age column from the mass table
        mass_table = np.delete(mass_table, 0, 0)

        # Remove the metallicity column from the spec table
        spec_table = np.delete(spec_table, 1, 0)

        # Convert the wavelength from Å to nm
        spec_table[1] = spec_table[1] * 0.1

        # For all ages, the lambda grid is the same.
        lambda_grid = np.unique(spec_table[1])

        # Creation of the age vs lambda flux table
        tmpList = []
        for wavelength in lambda_grid:
            [age_grid_orig, lambda_grid_orig, flux_orig] = \
                spec_table[:, spec_table[1, :] == wavelength]
            flux_orig = flux_orig * 10 * 1.e-7  # From erg/s^-1/Å to W/nm
            flux_regrid = interpolate.interp1d(age_grid_orig,
                                               flux_orig)(age_grid)

            tmpList.append(flux_regrid)
        flux_age = np.array(tmpList)

        base.add_ssp_m2005(SspM2005(imf, metallicity, age_grid,
                                    lambda_grid, mass_table, flux_age))

    print("\nDONE\n")
    print('#' * 78)

    ########################################################################
    # Dale and Helou 2002 templates insertion                              #
    ########################################################################
    print("3- Importing Dale and Helou 2002 templates\n")

    # Getting the alpha grid for the templates
    dhcal = np.genfromtxt(dh2002_dir + 'dhcal.dat')
    alpha_grid = dhcal[:, 1]

    # Getting the lambda grid for the templates (we checked that all share the
    # same grid).
    first_template = np.genfromtxt(dh2002_dir + 'irdh01.spec', skip_header=1)
    lambda_grid = first_template[:, 0] * 0.1  # Convert Å to nm

    templates = []

    for i in range(len(alpha_grid)):
        filename = dh2002_dir + 'irdh' + ("%02d" % (i + 1)) + '.spec'
        print("Importing %s..." % filename)
        table = np.genfromtxt(filename, skip_header=1)[:, 1]  # Luminosity
                                                              # column
        # The table give the luminosity density in Lsun/Å normalised to 1 Lsun
        # over the full spectrum. As we converted the wavelengths to nm, we
        # must multiply the density per 10 to keep the normalisation.
        table = table * 10
        templates.append(table)

    templates = np.array(templates)

    data = (alpha_grid, lambda_grid, templates)

    base.add_dh2002_infrared_templates(data)

    print("\nDONE\n")
    print('#' * 78)

    base.session.close_all()

if __name__ == '__main__':
    build_base()
