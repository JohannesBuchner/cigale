# -*- coding: utf-8 -*-
# Copyright (C) 2012, 2015 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

from distutils.command.build import build

from setuptools import find_packages, setup


class custom_build(build):
    def run(self):
        # Build the database.
        import database_builder
        database_builder.build_base()

        # Proceed with the build
        build.run(self)

entry_points = {
    'console_scripts': ['pcigale = pcigale:main',
                        'pcigale-plots = pcigale_plots:main',
                        'pcigale-filters = pcigale_filters:main']
}

setup(
    name="pcigale",
    version="0.8.0",
    packages=find_packages(exclude=["database_builder"]),

    install_requires=['numpy', 'scipy', 'sqlalchemy', 'matplotlib',
                      'configobj', 'astropy'],

    entry_points=entry_points,

    cmdclass={"build": custom_build},
    package_data={'pcigale': ['data/data.db'],
                  'pcigale_plots': ['data/CIGALE.png']},

    author="The CIGALE team",
    author_email="cigale@lam.fr",
    description="Python Code Investigating Galaxy Emission",
    license="CeCILL-V2",
    keywords="astrophysics, galaxy, SED fitting"
)
