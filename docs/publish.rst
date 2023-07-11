Publish package on PyPI
----------------------------------------

.. code-block::

   python3 -m pip install twine

.. ::

**Test PyPI** : for testing that all works, but not yet really publishing to a place where all the world is searching for python packages.

.. code-block::

   python3 -m twine upload --repository testpypi dist/*

.. ::

**The Real PyPI**

.. code-block::

   python3 -m twine upload dist/*

.. ::


`BEANSp on PyPI:  https://pypi.org/project/beansp <https://pypi.org/project/beansp>`_

