#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages, Extension
import glob


def get_version():
    """Get the version number of pySettle"""
    import beansp
    return beansp.__version__


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()


# add all libraries 
requirements = ["numpy", "matplotlib"]


setup(
    author="Adelle Goodwin",
    author_email='adelle.goodwin@monash.edu',
    python_requires='>=3.5.0',
    version=get_version(),
    version='0.8.1',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="Bayesian parameter Estimation of Accreting Neutron Stars",
    install_requires=requirements, #just list the requirements here
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='beans',
    name='beans',
    packages=find_packages(include=['beans', 'beans.*']),
    test_suite='tests',
    url='https://github.com/adellej/beans',
    zip_safe=False,
)
