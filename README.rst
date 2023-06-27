======
BEANSp
======


.. .. image:: https://img.shields.io/pypi/v/beans.svg
..         :target: https://pypi.python.org/pypi/beans

.. .. image:: https://img.shields.io/travis/adellej/beans.svg
..         :target: https://travis-ci.org/adellej/beans

.. .. image:: https://readthedocs.org/projects/beans/badge/?version=latest
..         :target: https://beans.readthedocs.io/en/latest/?badge=latest
..         :alt: Documentation Status



Bayesian Estimation of Accreting Neutron Star parameters
-----------------------------------------------------------------

* Free software: MIT license
* Documentation: https://beans-7.readthedocs.io/en/latest/
* Repo: https://github.com/adellej/beans


Features
--------

This software uses a Markov Chain Monte Carlo approach to match observations of an accreting neutron star in outburst with a simple ignition model to constrain parameters including the neutron star mass, radius, surface gravity, distance and system inclination, and accreted fuel composition. 

The code is written in Python 3, except for settle which is a C++ code with a python wrapper. It makes use of Dan Foreman-Mackey's python implementation of MCMC, emcee, available at https://github.com/dfm/emcee.

Credits
-------

Software written by Adelle Goodwin; for a full description see Goodwin et al. (2019, https://doi.org/10.1093/mnras/stz2638 or preprint at https://arxiv.org/pdf/1907.00996).

The algorithm is based on code written by Duncan Galloway, and depends on pySettle (https://github.com/adellej/pysettle), which was forked from the original settle written by Andrew Cumming.

Package installation and usage
------------------------------
BEANSp is on pyPI (https://pypi.org/project/beansp/) so installation is easy - either straight or in virtual environment:

   .. code-block::
   
      pip install beansp
  
   .. ::
   
   .. code-block::
   
      from beansp.beans import Beans 

(Please refer to `this simple test script <https://github.com/adellej/beans/blob/master/tests/test_sft_beans.py>`_ as an example.)

Build and installation from this github repository
--------------------------------------------------

Please refer to `build instructions <https://github.com/adellej/beans/blob/master/BUILD.rst>`_.


