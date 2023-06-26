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
- h5py v2.10.0 or above. To install h5py v2.10.0 on Ubuntu you will need to use pip wheels, e.g. $ pip install h5py==2.10.0
- pytest 
- pickle
-pathlib

Once you have downloaded the source code, and have all of the dependencies, navigate to beans/settle and you will need to compile settle. Compiling settle requires a different command depending if you are using Mac or Linux. For Mac type:

.. code-block:: console

    $ make mac

For Linux type:

.. code-block:: console

    $ make linux


This will create a file called libsettle.so which is a precompiled binary that is used by beans to run settle. The makefile may need to be edited if either of these methods do not compile the code.

I also recommend that you add BEANS to your python path so that you can access it from any location, and so that the test suite will work. You can do this by opening your bashrc or bash_profile file and adding the following line:

.. code-block:: console

    export PYTHONPATH=[PATH TO BEANS]:${PYTHONPATH}

where you replace [PATH TO BEANS] with the path to where you have placed the BEANS directory on your computer. For example, on my computer this becomes:

.. code-block:: console

    export PYTHONPATH=${HOME}/BEANS/beans:${PYTHONPATH}


Now that you have compiled settle I recommend you run the test suite to check you have all the required dependencies and the code is operating as expected. To do this navigate to beans/ and type:

.. code-block:: console

    $ pytest

If the tests all pass then you are good to go!

