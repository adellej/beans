=====
beans
=====


.. .. image:: https://img.shields.io/pypi/v/beans.svg
..         :target: https://pypi.python.org/pypi/beans

.. .. image:: https://img.shields.io/travis/adellej/beans.svg
..         :target: https://travis-ci.org/adellej/beans

.. .. image:: https://readthedocs.org/projects/beans/badge/?version=latest
..         :target: https://beans.readthedocs.io/en/latest/?badge=latest
..         :alt: Documentation Status




Bayesian parameter Estimation of Accreting Neutron Stars
--------------------------------------------------------

* Free software: MIT license
* Documentation: https://beans-7.readthedocs.io/en/latest/


Features
--------

This software uses a Markov Chain Monte Carlo approach to match observations of an accreting neutron star in outburst with a simple ignition model to predict unobservable parameters such as neutron star mass, radius, surface gravity, distance and inclination of the source, and accreted fuel composition. The code is all written in Python 3, except for settle which is a c++ code with a python wrapper. It makes use of Dan Foreman-Mackey's python implementation of MCMC, emcee, available here - https://github.com/dfm/emcee.

Credits
-------

Software written by Adelle Goodwin. See Goodwin et al. (2019) - https://arxiv.org/pdf/1907.00996.

This code makes use of settle (written by Andrew Cumming), was based on code written by Duncan Galloway, and uses Dan Foreman-Mackey's python implementation of emcee, MCMC. 

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
