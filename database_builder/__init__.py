# -*- coding: utf-8 -*-
# Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""
This script is used the build pcigale internal database containing:
- The various filter transmission tables;
- The Maraston 2005 single stellar population (SSP) data;
- The Dale and Helou 2002 infra-red templates.

"""
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../'))
import glob
import io
import itertools
import numpy as np
from scipy import interpolate
import scipy.constants as cst
from pcigale.data import (Database, Filter, SspM2005, SspBC03, AgnFritz2006,
                          DALE2014, DL2007)


def read_bc03_ssp(filename):
    """Read a Bruzual and Charlot 2003 ASCII SSP file

    The ASCII SSP files of Bruzual and Charlot 2003 have se special structure.
    A vector is stored with the number of values followed by the values
    separated by a space (or a carriage return). There are the time vector, 5
    (for Chabrier IMF) or 6 lines (for Salpeter IMF) that we don't care of,
    then the wavelength vector, then the luminosity vectors, each followed by
    a 52 value table, then a bunch of other table of information that are also
    in the *colors files.

    Parameters
    ----------
    filename : string

    Returns
    -------
    time_grid: numpy 1D array of floats
              Vector of the time grid of the SSP in Myr.
    wavelength: numpy 1D array of floats
                Vector of the wavelength grid of the SSP in nm.
    spectra: numpy 2D array of floats
             Array containing the SSP spectra, first axis is the wavelength,
             second one is the time.

    """

    def file_structure_generator():
        """Generator used to identify table lines in the SSP file

        In the SSP file, the vectors are store one next to the other, but
        there are 5 informational lines after the time vector. We use this
        generator to the if we are on lines to read or not.
        """
        if "chab" in filename:
            bad_line_number = 5
        else:
            bad_line_number = 6
        yield("data")
        for i in range(bad_line_number):
            yield("bad")
        while True:
            yield("data")

    file_structure = file_structure_generator()
    # Are we in a data line or a bad one.
    what_line = file_structure.next()
    # Variable conting, in reverse order, the number of value still to
    # read for the read vector.
    counter = 0

    time_grid = []
    full_table = []
    tmp_table = []

    with open(filename) as file_:
        # We read the file line by line.
        for line in file_:
            if what_line == "data":
                # If we are in a "data" line, we analyse each number.
                for item in line.split():
                    if counter == 0:
                        # If counter is 0, then we are not reading a vector
                        # and the first number is the length of the next
                        # vector.
                        counter = int(item)
                    else:
                        # If counter > 0, we are currently reading a vector.
                        tmp_table.append(float(item))
                        counter -= 1
                        if counter == 0:
                            # We reached the end of the vector. If we have not
                            # yet store the time grid (the first table) we are
                            # currently reading it.
                            if time_grid == []:
                                time_grid = tmp_table[:]
                            # Else, we store the vector in the full table,
                            # only if its length is superior to 250 to get rid
                            # of the 52 item unknown vector and the 221 (time
                            # grid length) item vectors at the end of the
                            # file.
                            elif len(tmp_table) > 250:
                                full_table.append(tmp_table[:])

                            tmp_table = []

            # If at the end of a line, we have finished reading a vector, it's
            # time to change to the next structure context.
            if counter == 0:
                what_line = file_structure.next()

    # The time grid is in year, we want Myr.
    time_grid = np.array(time_grid, dtype=float)
    time_grid = time_grid * 1.e-6

    # The first "long" vector encountered is the wavelength grid. The value
    # are in Ångström, we convert it to nano-meter.
    wavelength = np.array(full_table.pop(0), dtype=float)
    wavelength = wavelength * 0.1

    # The luminosities are in Solar luminosity (3.826.10^33 ergs.s-1) per
    # Ångström, we convert it to W/nm.
    luminosity = np.array(full_table, dtype=float)
    luminosity = luminosity * 3.826e27
    # Transposition to have the time in the second axis.
    luminosity = luminosity.transpose()

    # In the SSP, the time grid begins at 0, but not in the *colors file, so
    # we remove t=0 from the SSP.
    return time_grid[1:], wavelength, luminosity[:, 1:]


def build_filters(base):
    filters_dir = os.path.join(os.path.dirname(__file__), 'filters/')
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


def build_m2005(base):
    m2005_dir = os.path.join(os.path.dirname(__file__), 'maraston2005/')

    # Age grid (1 Myr to 13.7 Gyr with 1 Myr step)
    age_grid = np.arange(1, 13701)

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
            imf = 'krou'
            mass_table = np.copy(kroupa_mass)
        elif 'ssz' in spec_file:
            imf = 'salp'
            mass_table = np.copy(salpeter_mass)
        else:
            raise ValueError('Unknown IMF!!!')

        # Keep only the actual metallicity values in the mass table
        # we don't take the first column which contains metallicity
        mass_table = mass_table[1:, mass_table[0] == metallicity]

        # Interpolate the mass table over the new age grid. We multiply per
        # 1000 because the time in Maraston files is given in Gyr.
        mass_table = interpolate.interp1d(mass_table[0] * 1000,
                                          mass_table)(age_grid)

        # Remove the age column from the mass table
        mass_table = np.delete(mass_table, 0, 0)

        # Remove the metallicity column from the spec table
        spec_table = np.delete(spec_table, 1, 0)

        # Convert the wavelength from Å to nm
        spec_table[1] = spec_table[1] * 0.1

        # For all ages, the lambda grid is the same.
        lambda_grid = np.unique(spec_table[1])

        # Creation of the age vs lambda flux table
        tmp_list = []
        for wavelength in lambda_grid:
            [age_grid_orig, lambda_grid_orig, flux_orig] = \
                spec_table[:, spec_table[1, :] == wavelength]
            flux_orig = flux_orig * 10 * 1.e-7  # From erg/s^-1/Å to W/nm
            age_grid_orig = age_grid_orig * 1000  # Gyr to Myr
            flux_regrid = interpolate.interp1d(age_grid_orig,
                                               flux_orig)(age_grid)

            tmp_list.append(flux_regrid)
        flux_age = np.array(tmp_list)

        # Use Z value for metallicity, not log([Z/H])
        metallicity = {-1.35: 0.001,
                       -0.33: 0.01,
                       0.0: 0.02,
                       0.35: 0.04}[metallicity]

        base.add_ssp_m2005(SspM2005(imf, metallicity, age_grid,
                                    lambda_grid, mass_table, flux_age))


def build_bc2003(base):
    bc03_dir = os.path.join(os.path.dirname(__file__), 'bc03//')

    # Time grid (1 Myr to 20 Gyr with 1 Myr step)
    time_grid = np.arange(1, 20000)

    # Metallicities associated to each key
    metallicity = {
        "m22": 0.0001,
        "m32": 0.0004,
        "m42": 0.004,
        "m52": 0.008,
        "m62": 0.02,
        "m72": 0.05
    }

    for key, imf in itertools.product(metallicity, ["salp", "chab"]):
        base_filename = bc03_dir + "bc2003_lr_" + key + "_" + imf + "_ssp"
        ssp_filename = base_filename + ".ised_ASCII"
        color3_filename = base_filename + ".3color"
        color4_filename = base_filename + ".4color"

        print("Importing %s..." % base_filename)

        # Read the desired information from the color files
        color_table = []
        color3_table = np.genfromtxt(color3_filename).transpose()
        color4_table = np.genfromtxt(color4_filename).transpose()
        color_table.append(color4_table[6])        # Mstar
        color_table.append(color4_table[7])        # Mgas
        color_table.append(10 ** color3_table[5])  # NLy
        color_table.append(color3_table[1])        # B4000
        color_table.append(color3_table[2])        # B4_VN
        color_table.append(color3_table[3])        # B4_SDSS
        color_table.append(color3_table[4])        # B(912)

        color_table = np.array(color_table)

        ssp_time, ssp_wave, ssp_lumin = read_bc03_ssp(ssp_filename)

        # Regrid the SSP data to the evenly spaced time grid.
        color_table = interpolate.interp1d(ssp_time, color_table)(time_grid)
        ssp_lumin = interpolate.interp1d(ssp_time,
                                         ssp_lumin)(time_grid)

        base.add_ssp_bc03(SspBC03(
            imf,
            metallicity[key],
            time_grid,
            ssp_wave,
            color_table,
            ssp_lumin
        ))


def build_dh2002(base):
    dh2002_dir = os.path.join(os.path.dirname(__file__), 'dh2002/')

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
        # The table give the luminosity density in Lsun/Å normalised to 1 Lsun
        # over the full spectrum. As we converted the wavelengths to nm, we
        # must multiply the density per 10 to keep the normalisation.
        table = table * 10
        templates.append(table)

    templates = np.array(templates)

    data = (alpha_grid, lambda_grid, templates)

    base.add_dh2002_infrared_templates(data)

def build_dale2014(base):

    dh2002_dir = os.path.join(os.path.dirname(__file__), 'dh2002/')
    dale2014_dir = os.path.join(os.path.dirname(__file__), 'dale2014/')

    # Getting the alpha grid for the templates
    d14cal = np.genfromtxt(dh2002_dir + 'dhcal.dat')
    alpha_grid = d14cal[:, 1]

    # Getting the lambda grid for the templates and convert from microns to nm.
    first_template = np.genfromtxt(dale2014_dir + 'spectra.0.00AGN.dat')
    wave = first_template[:, 0] * 1E3

    # Getting the stellar emission and interpolate it at the same wavelength grid
    stell_emission_file = np.genfromtxt(dale2014_dir + 'stellar_SED_age13Gyr_tau10Gyr.spec')
    # A -> to nm
    wave_stell = stell_emission_file[:,0] * 0.1
    # W/A -> W/nm
    stell_emission = stell_emission_file[:,1] * 10
    stell_emission_interp = np.interp(wave,wave_stell,stell_emission)

    # The models are in nuFnu and contain stellar emission.
    # We convert this to W/nm and remove the stellar emission.

    # Emission from dust heated by SB
    fraction = 0.0
    filename = dale2014_dir + "spectra.0.00AGN.dat"
    print("Importing {}...".format(filename))
    datafile = open(filename)
    data = "".join(datafile.readlines())
    datafile.close()

    for al in range(1,len(alpha_grid),1):
        lumin_with_stell = np.genfromtxt(io.BytesIO(data.encode()), usecols=(al))
        lumin_with_stell = pow(10,lumin_with_stell) / wave
        constant = lumin_with_stell[7] / stell_emission_interp[7]
        lumin = lumin_with_stell - stell_emission_interp * constant
        lumin[lumin<0] = 0
        lumin[wave<2E3] = 0
        norm = np.trapz(lumin, x = wave)
        lumin = lumin/norm

        base.add_dale2014(DALE2014(fraction, alpha_grid[al-1], wave, lumin))

    # Emission from dust heated by AGN - Quasar template
    fraction = 1.0
    filename = dale2014_dir + "spectra.1.00AGN.dat"
    print("Importing {}...".format(filename))
    datafile = open(filename)
    data = "".join(datafile.readlines())
    datafile.close()

    for al in range(1,len(alpha_grid),1):
        lumin_quasar = np.genfromtxt(io.BytesIO(data.encode()), usecols=(al))
        lumin_quasar = pow(10,lumin_quasar) / wave
        lumin_quasar[lumin_quasar<0] = 0
        norm = np.trapz(lumin_quasar, x = wave)
        lumin_quasar = lumin_quasar/norm

        base.add_dale2014(DALE2014(fraction, alpha_grid[al-1], wave, lumin_quasar))


def build_dl2007(base):
    dl2007_dir = os.path.join(os.path.dirname(__file__), 'dl2007/')

    qpah = {
        "00": 0.47,
        "10": 1.12,
        "20": 1.77,
        "30": 2.50,
        "40": 3.19,
        "50": 3.90,
        "60": 4.58
    }

    umaximum = ["1e3", "1e4", "1e5", "1e6"]
    uminimum = ["0.10", "0.15", "0.20", "0.30", "0.40", "0.50", "0.70",
                "0.80", "1.00", "1.20", "1.50", "2.00", "2.50", "3.00",
                "4.00", "5.00", "7.00", "8.00", "10.0", "12.0", "15.0",
                "20.0", "25.0"]

    # Here we obtain the wavelength beforehand to avoid reading it each time.
    datafile = open(dl2007_dir + "U{}/U{}_{}_MW3.1_{}.txt".format(umaximum[0],
                                                                  umaximum[0],
                                                                  umaximum[0],
                                                                  "00"))
    data = "".join(datafile.readlines()[-1001:])
    datafile.close()

    wave = np.genfromtxt(io.BytesIO(data.encode()), usecols=(0))
    # For some reason wavelengths are decreasing in the model files
    wave = wave[::-1]
    # We convert wavelengths from μm to nm
    wave *= 1000.

    # The models are in Jy cm² sr¯¹ H¯¹. We convert this to W/nm.
    conv = 4. * np.pi * 1e-30 / cst.m_p * cst.c / (wave * wave) * 1e9

    for model in sorted(qpah.keys()):
        for umin in uminimum:
            filename = dl2007_dir + "U{}/U{}_{}_MW3.1_{}.txt".format(umin,
                                                                     umin,
                                                                     umin,
                                                                     model)
            print("Importing {}...".format(filename))
            datafile = open(filename)
            data = "".join(datafile.readlines()[-1001:])
            datafile.close()
            lumin = np.genfromtxt(io.BytesIO(data.encode()), usecols=(2))
            # For some reason fluxes are decreasing in the model files
            lumin = lumin[::-1]
            # Conversion from Jy cm² sr¯¹ H¯¹ to W/nm
            lumin *= conv

            base.add_dl2007(DL2007(qpah[model], umin, umin, wave, lumin))
            for umax in umaximum:
                filename = dl2007_dir + "U{}/U{}_{}_MW3.1_{}.txt".format(umin,
                                                                         umin,
                                                                         umax,
                                                                         model)
                print("Importing {}...".format(filename))
                datafile = open(filename)
                data = "".join(datafile.readlines()[-1001:])
                datafile.close()
                lumin = np.genfromtxt(io.BytesIO(data.encode()), usecols=(2))
                # For some reason fluxes are decreasing in the model files
                lumin = lumin[::-1]

                # Conversion from Jy cm² sr¯¹ H¯¹ to W/nm
                lumin *= conv

                base.add_dl2007(DL2007(qpah[model], umin, umax, wave, lumin))


def build_fritz2006(base):
    fritz2006_dir = os.path.join(os.path.dirname(__file__), 'fritz2006/')

    model_list = np.genfromtxt(fritz2006_dir + "fritz.dat")

    for model_line in model_list:

        (model_nb, agn_type, r_ratio, tau,
         beta, gamma, theta, psy) = model_line

        # Convert some floats to int
        model_nb = int(model_nb)
        agn_type = int(agn_type)

        wave, lumin = np.genfromtxt("{}AGN_fritz{}.spec".format(fritz2006_dir,
                                                                model_nb),
                                    skip_header=1).transpose()

        # Convert the wavelength from Å to nm
        wave = wave * 0.1

        # Convert the luminosity from erg/s^-1/Å to W/nm
        lumin = lumin * 10 * 1.e-7

        agn = AgnFritz2006(model_nb, agn_type, r_ratio, tau, beta, gamma,
                           theta, psy, wave, lumin)

        base.add_fritz2006_agn(agn)


def build_base():
    base = Database(writable=True)
    base.upgrade_base()

    print('#' * 78)
    print("1- Importing filters...\n")
    build_filters(base)
    print("\nDONE\n")
    print('#' * 78)

    print("2- Importing Maraston 2005 SSP\n")
    build_m2005(base)
    print("\nDONE\n")
    print('#' * 78)

    print("3- Importing Bruzual and Charlot 2003 SSP\n")
    build_bc2003(base)
    print("\nDONE\n")
    print('#' * 78)

    print("4- Importing Dale and Helou (2002) templates\n")
    build_dh2002(base)
    print("\nDONE\n")
    print('#' * 78)

    print("5- Importing Draine and Li (2007) templates\n")
    build_dl2007(base)
    print("\nDONE\n")
    print('#' * 78)

    print("6- Importing Fritz et al. (2006) models\n")
    build_fritz2006(base)
    print("\nDONE\n")
    print('#' * 78)

    print("7- Importing Dale et al (2014) templates\n")
    build_dale2014(base)
    print("\nDONE\n")
    print('#' * 78)

    base.session.close_all()


if __name__ == '__main__':
    build_base()
