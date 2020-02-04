#!/usr/bin/env python
# -*- coding: utf-8 -*-

import setuptools

setuptools.setup(
    name="ktkit",
    version="0.0.0",
    url="https://github.com/samstudio8/ktkit",

    description="A command line tool for playing with kraken2 output files",
    long_description="",

    author="Sam Nicholls",
    author_email="sam@samnicholls.net",

    maintainer="Sam Nicholls",
    maintainer_email="sam@samnicholls.net",

    packages=setuptools.find_packages(),
    include_package_data=True,

    install_requires=[
    ],

    entry_points = {
        "console_scripts": [
            "ktkit=ktkit:cli",
        ]
    },

    classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
    ],

)
