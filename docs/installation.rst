.. highlight:: shell

============
Installation
============

From source
------------

The source for beans can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/adellej/beans

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/adellej/beans/tarball/master

.. Once you have a copy of the source, you can install it with:

.. .. code-block:: console

..     $ python setup.py install


.. _Github repo: https://github.com/adellej/beans
.. _tarball: https://github.com/adellej/beans/tarball/master

The first thing to do is check that you have all of the required python dependencies installed. These are listed below (to install either use conda or pip):

- emcee v3.0 or above
- matplotlib
- numpy 
- corner
- random
- math
- astropy
- scipy
- tables
- chainconsumer
- multiprocessing
- os
- time
- h5py


Once you have downloaded the source code, and have all of the dependencies, navigate to beans/settle and you will need to compile settle. Settle can be compiled by typing the following (provided you have gcc):

.. code-block:: console

    $ c++ -Ofast -o settle *.c *.cc -lm

This will create a file called libsettle.so which is a precompiled binary that is used by beans to run settle. 

Now that you have compiled settle I recommend you run the test suite to check you have all the required dependencies and the code is operating as expected. To do this navigate to beans/ and type:

.. code-block:: console

    $ pytest

If the tests all pass then you are good to go!

