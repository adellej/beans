#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages
import re


def get_property(prop, project):
    result = re.search(r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop),
                       open(project + '/__init__.py').read())
    return result.group(1)


def get_version():
    """Get the version number of BEANSp"""
    # ## the original inspired by Paul's Aegean package
    # ## does not work with siplified imports trick in __init__.py
    import beansp
    return beansp.__version__


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

# add all libraries
requirements = ["numpy>=1.16",
                "matplotlib",
                "emcee>=3.0",
                "corner",
                "astropy",
                "scipy",
                "tables",
                "chainconsumer",
                "h5py>=2.10.0",
                "pySettle"]

package_name = 'beansp'

setup(
    author=get_property('__author__', package_name),
    author_email=get_property('__email__', package_name),
    name=package_name,
    python_requires='>=3.6.0',
    version=get_property('__version__', package_name),
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    description="Bayesian Estimation of Accreting Neutron Stars parameters",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='beans',
    packages=find_packages(include=[package_name]),
    package_data={package_name: ['data/*']},
    test_suite='tests',
    url='https://github.com/adellej/beans',
    zip_safe=False,
)
