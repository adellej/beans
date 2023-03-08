#! /usr/bin/env python
"""
Setup for settle, the beans edition
"""
import os
import subprocess
from setuptools import setup
# , find_packages, Extension

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...


def read(fname):
    """Read a file"""
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def get_version():
    """Get the version number of pySettle"""
    build_install_libsettle()
    import pySettle
    return pySettle.__version__


# The C/C++ coded shared library had it's own makefike
# to build. Hence we call that rather than duplicating that
# within setuptools


def build_install_libsettle():
    """compile and install C/C++ library libsettle"""
    result = subprocess.run("make -C libsettle install",
                            shell=True,
                            check=True,
                            capture_output=True)
    print(result.stdout, result.stderr)
    if result.returncode == 0:
        print("libsettle built OK.")
    else:
        print("libsettle bulid FAILED!")

    return result.returncode


# What is the reasonable minimal version of numpy to start with?
reqs = ['numpy>=1.16']

setup(
    name="pySettle",
    packages=['pySettle'],
    # packages=find_packages(exclude='test'),
    version=get_version(),

    description="Computes ignition conditions for Type I X-ray bursts\
                 using a multi-zone model of the Neutron star accreting layer",
    long_description=read('README.rst'),
    license='MIT',
    # the author of the original settle packege
    # https://github.com/andrewcumming/settle
    author="Andrew Cumming",
    author_email='andrew.cumming@mcgill.ca',
    # include_package_data=True,
    keywords='settle',

    install_requires=['numpy>=1.16'],

    python_requires='>=3.6',
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
