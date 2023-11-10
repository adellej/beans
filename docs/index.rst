Welcome to beansp's documentation!
======================================

This software uses a Markov Chain Monte Carlo approach to match observations of an accreting neutron star in outburst with a simple ignition model to constrain parameters including the neutron star mass, radius, surface gravity, distance and system inclination, and accreted fuel composition. 

The code is written in Python 3, except for settle which is a C++ code
with a python wrapper (now implemented as
pySettle,
available via pypi or at
https://github.com/adellej/pysettle.

beansp makes use of Dan Foreman-Mackey's python implementation of MCMC, emcee, available at https://github.com/dfm/emcee.

*"The p is silent"*

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   readme
   installation
   build
   usage
   modules
   contributing

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
