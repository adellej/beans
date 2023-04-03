#! /usr/bin/env python
"""
Setup for settle, the beans edition
"""
import os
import subprocess
import sysconfig
from setuptools import setup, Extension
# , glob
# , find_packages

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...


def read(fname):
    """Read a file"""
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def get_version():
    """Get the version number of pySettle"""
    # Not used now - building the settle dynamic library using setuptools
    # build_install_libsettle()
    print_platform()
    import pySettle
    return pySettle.__version__


# The C/C++ coded shared library had it's own makefike
# to build. This function is left here, for the case
# there ever is a need to use this more transparent build method
# rather than setuptools internal approach


def build_install_libsettle():
    """compile and install C/C++ library libsettle
       using a custom Makefile
    """
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


def print_platform():
    system_platform = sysconfig.get_platform()
    print("========================")
    print("Platform is " + system_platform)
    print("========================")


# What is the reasonable minimal version of numpy to start with?
reqs = ['numpy>=1.16']

setup(
    name="pySettle",
    # prefer to set the package name explicitly, lets KISS!
    packages=['pySettle'],
    # the following does not work - will not exclude test procedures
    # used MANIFEST.in instead.
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

    install_requires=reqs,

    python_requires='>=3.6',

    classifiers=[
        # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable"
        # as the current state of your package
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],

    ext_modules=[
        Extension(
            "settle",
            # this compiles all *.c and *.cc source files
            # and links them into dybamic library
            #   glob.glob("libsettle/*.c*")
            # let's be explicit to have full control
            sources=["libsettle/gasdev.c", "libsettle/nrutil.c",
                     "libsettle/root.c", "libsettle/eos.cc",
                     "libsettle/odeint.cc", "libsettle/settle.cc",
                     "libsettle/spline.cc", "libsettle/useful.cc"],
            extra_compile_args=["-Ofast",
                                "-Wno-unused-but-set-variable",
                                "-Wno-unused-parameter",
                                "-Wno-unused-variable",
                                "-Wno-unused-result"]
        )],
)
