# -*- coding: utf-8 -*-
# Copyright (C) 2015 Institute of Astronomy
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Médéric Boquien

import argparse
from astropy.table import Table, Column
import astropy.units as u
import matplotlib.pyplot as plt
import multiprocessing as mp
import numpy as np
from pcigale.data import Database, Filter
import sys

__version__ = "0.1-alpha"


def list_filters():
    """Print the list of filters in the pcigale database.
    """
    with Database() as base:
        filters = {name: base.get_filter(name) for name in
                   base.get_filter_names()}

    name = Column(data=[filters[f].name for f in filters], name='Name')
    description = Column(data=[filters[f].description for f in filters],
                         name='Description')
    wl = Column(data=[filters[f].effective_wavelength for f in filters],
                name='Effective Wavelength', unit=u.nm, format='%d')
    filter_type = Column(data=[filters[f].trans_type for f in filters],
                         name='Type')
    samples = Column(data=[filters[f].trans_table[0].size for f in filters],
                     name="Points")

    t = Table()
    t.add_columns([name, description, wl, filter_type, samples])
    t.sort(['Effective Wavelength'])
    t.pprint(max_lines=-1, max_width=-1)


def add_filters(fnames):
    """Add filters to the pcigale database.
    """
    with Database(writable=True) as base:
        for fname in fnames:
            with open(fname, 'r') as f_fname:
                filter_name = f_fname.readline().strip('# \n\t')
                filter_type = f_fname.readline().strip('# \n\t')
                filter_description = f_fname.readline().strip('# \n\t')
            filter_table = np.genfromtxt(fname)
            # The table is transposed to have table[0] containing the
            # wavelength and table[1] containing the transmission.
            filter_table = filter_table.transpose()
            # We convert the wavelength from Å to nm.
            filter_table[0] *= 0.1

            print("Importing {}... ({} points)".format(filter_name,
                                                       filter_table.shape[1]))

            new_filter = Filter(filter_name, filter_description, filter_type,
                                filter_table)

            # We normalise the filter and compute the effective wavelength.
            # If the filter is a pseudo-filter used to compute line fluxes, it
            # should not be normalised.
            if not filter_name.startswith('PSEUDO'):
                new_filter.normalise()
            else:
                new_filter.effective_wavelength = np.mean(
                    filter_table[0][filter_table[1] > 0]
                )

            base.add_filter(new_filter)


def del_filters(fnames):
    """Delete filters from the pcigale database
    """
    with Database(writable=True) as base:
        names = base.get_filter_names()
        for fname in fnames:
            if fname in names:
                base.del_filter(fname)
                print("Removing filter {}".format(fname))
            else:
                print("Filter {} not in the database".format(fname))


def worker_plot(fname):
    """Worker to plot filter transmission curves in parallel

    Parameters
    ----------
    fname: string
        Name of the filter to be plotted
    """
    with Database() as base:
        _filter = base.get_filter(fname)
    plt.clf()
    plt.plot(_filter.trans_table[0], _filter.trans_table[1], color='k')
    plt.xlim(_filter.trans_table[0][0], _filter.trans_table[0][-1])
    plt.minorticks_on()
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('Relative transmission')
    plt.title("{} filter".format(fname))
    plt.tight_layout()
    plt.savefig("{}.pdf".format(fname))


def plot_filters(fnames):
    """Plot the filters provided as parameters. If not filter is given, then
    plot all the filters.
    """
    if len(fnames) == 0:
        with Database() as base:
            fnames = base.get_filter_names()
    with mp.Pool(processes=mp.cpu_count()) as pool:
        pool.map(worker_plot, fnames)


def main():

    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(help="List of commands")

    list_parser = subparsers.add_parser('list', help=list_filters.__doc__)
    list_parser.set_defaults(parser='list')

    add_parser = subparsers.add_parser('add', help=add_filters.__doc__)
    add_parser.add_argument('names', nargs='+', help="List of file names")
    add_parser.set_defaults(parser='add')

    del_parser = subparsers.add_parser('del', help=del_filters.__doc__)
    del_parser.add_argument('names', nargs='+', help="List of filter names")
    del_parser.set_defaults(parser='del')

    plot_parser = subparsers.add_parser('plot', help=plot_filters.__doc__)
    plot_parser.add_argument('names', nargs='*', help="List of filter names")
    plot_parser.set_defaults(parser='plot')

    if len(sys.argv) == 1:
        parser.print_usage()
    else:
        args = parser.parse_args()
        if args.parser == 'list':
            list_filters()
        elif args.parser == 'add':
            add_filters(args.names)
        elif args.parser == 'del':
            del_filters(args.names)
        elif args.parser == 'plot':
            plot_filters(args.names)
