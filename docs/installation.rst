.. highlight:: shell

============
Installation
============

BEANSp is on `pyPI`_ so installation is easy - either system-wide or in virtual environment:

.. code-block:: console
    pip install beansp

And then just import

.. code-block:: console

    from beansp import Beans

(Please refer to `this simple test script`_ as an example.)

.. _pypi: https://pypi.org/project/beansp
.. _this simple test script: https://github.com/adellej/beans/blob/master/tests/test_sft_beans.py

From source
------------

The source for beans can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/adellej/beans

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/adellej/beans/tarball/master

.. _Github repo: https://github.com/adellej/beans
.. _tarball: https://github.com/adellej/beans/tarball/master

Once you have a copy of the source, you can install it following the same
`build instructions`_ as for pySettle.

.. _build instructions: https://github.com/adellej/pysettle/blob/master/BUILD.rst

Testing
-------

Once you have compiled settle we recommend you run the test suite to check you have all the required dependencies and the code is operating as expected. To do this navigate to the top-level directory and type:

.. code-block:: console

    $ pytest

If the tests all pass then you are good to go!

