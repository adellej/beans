Welcome to beansp's documentation!
======================================

This software uses a Markov Chain Monte Carlo approach to match observations of an accreting neutron star in outburst with a simple ignition model to constrain parameters including the neutron star mass, radius, surface gravity, distance and system inclination, and accreted fuel composition. 

The code is written in Python 3, except for settle which is a C++ code with a python wrapper. It makes use of Dan Foreman-Mackey's python implementation of MCMC, emcee, available at https://github.com/dfm/emcee.

Credits
=======

Software written by Adelle Goodwin; for a full description see Goodwin et al. (2019, https://doi.org/10.1093/mnras/stz2638 or preprint at https://arxiv.org/pdf/1907.00996).

The algorithm is based on code written by Duncan Galloway, and depends on pySettle (https://github.com/adellej/pysettle), which was forked from the original settle written by Andrew Cumming.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   readme
   installation
   usage
   modules
   contributing
   history

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
