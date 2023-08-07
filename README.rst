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

Software written by Adelle Goodwin and Duncan Galloway; for a full description see Goodwin et al. (2019, https://doi.org/10.1093/mnras/stz2638 or preprint at https://arxiv.org/pdf/1907.00996).

The original algorithm is described in Galloway & Cumming (2006, 
https://iopscience.iop.org/article/10.1086/507598).

pySettle was forked from the original settle written by Andrew Cumming and
available at https://github.com/andrewcumming/settle

Package installation and usage
------------------------------
BEANSp is on pyPI (https://pypi.org/project/beansp/) so installation is
easy - either system-wide, or in virtual environment:

.. code-block:: console
   
    pip install beansp
  
You can then import the main ``Beans`` module as follows:
   
.. code-block:: python
   
    from beansp import Beans 

(Please refer to `this simple test script <https://github.com/adellej/beans/blob/master/tests/test_sft_beans.py>`_ as an example.)

Build and installation from this github repository
--------------------------------------------------

Please refer to `build instructions <https://github.com/adellej/beans/blob/master/BUILD.rst>`_.


