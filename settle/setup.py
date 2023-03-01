#! /usr/bin/env python
"""
Setup for settle, the beans edition
"""
import os
from setuptools import setup, Extension

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...


def read(fname):
    """Read a file"""
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def get_version():
    """Get the version number of settle"""
    from . import __version__
    return __version_


# What is the reasonable version of numpy to start with?
reqs = ['numpy>=1.16']

extensions = [Extension('libsettle',
                       ["src/library.c"],
                       depends=["src/library.h"],
                       include_dirs=["src"],
              ),
]

setup(
    name="settle",
    #    packages=find_packages(exclude='test'),

    description="Computes ignition conditions for Type I X-ray bursts using a multi-zone model of the Neutron star accreting layer",
    long_description=read('README.rst'),
    licence='MIT',
    author="Andrew Cumming",
    author_email='Andrew.Cumming@gmail.com',
    # include_package_data=True,
    keywords='settle',

    python_requires='!=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*',
    classifiers=[
        # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable"
        # as the current state of your package
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
)
