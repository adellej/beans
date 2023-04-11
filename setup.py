#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages
# , find_namespace_packages
# , Extension
# import glob


def get_version():
    """Get the version number of BEANSp"""
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

setup(
    author="Adelle Goodwin",
    author_email='adelle.goodwin@monash.edu',
    python_requires='>=3.6.0',
    version=get_version(),
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
    description="Bayesian parameter Estimation of Accreting Neutron Stars",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='beans',
    name='beansp',
    packages=find_packages(include=['beansp']),
    package_data={'beansp': ['data/*']},
    test_suite='tests',
    url='https://github.com/adellej/beans',
    zip_safe=False,
)
