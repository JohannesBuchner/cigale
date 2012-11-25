# -*- coding: utf-8 -*-
"""
Copyright (C) 2012 Centre de données Astrophysiques de Marseille
Licensed under the CeCILL-v2 licence - see Licence_CeCILL_V2-en.txt

@author: Yannick Roehlly <yannick.roehlly@oamp.fr>

"""

from setuptools import setup, find_packages
from distutils.command.build import build



class custom_build(build):
    def run(self):
        # Build the database.
        import database_builder
        database_builder.build_base()

        # Proceed with the build
        build.run(self)

setup(
    name = "pcigale",
    version = "0.1a",
    packages = find_packages(exclude=["database_builder"]),

    install_requires = ['numpy', 'scipy', 'sqlalchemy'],
    requires = ['numpy','scipy', 'sqlalchemy'],

    include_package_data=True,
    cmdclass={"build": custom_build},
    package_data = {'': ['*.db']},

    author = "Yannick Roehlly",
    author_email = "yannick.roehlly@oamp.fr",
    description = "Python Code Investigating Galaxy Emission",
    license = "CeCILL-V2",
    keywords = "astrophysics, galaxy, SED fitting"
)
