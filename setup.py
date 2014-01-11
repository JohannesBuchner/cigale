# -*- coding: utf-8 -*-
# Copyright (C) 2012, 2013 Centre de données Astrophysiques de Marseille
# Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt
# Author: Yannick Roehlly

from setuptools import setup, find_packages
from distutils.command.build import build


class custom_build(build):
    def run(self):
        # Build the database.
        import database_builder
        database_builder.build_base()

        # Proceed with the build
        build.run(self)

entry_points = {
    'console_scripts': ['pcigale = pcigale:main']
}

setup(
    name="pcigale",
    version="0.1a",
    packages=find_packages(exclude=["database_builder"]),

    install_requires=['numpy', 'scipy', 'sqlalchemy', 'atpy', 'matplotlib',
                      'configobj', 'progressbar', 'pyfits', 'astropy'],

    entry_points=entry_points,

    include_package_data=True,
    cmdclass={"build": custom_build},
    package_data={'': ['*.db']},

    author="Yannick Roehlly",
    author_email="yannick.roehlly@oamp.fr",
    description="Python Code Investigating Galaxy Emission",
    license="CeCILL-V2",
    keywords="astrophysics, galaxy, SED fitting"

)
