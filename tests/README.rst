===========
beans tests
===========

Bayesian Estimation of Accreting Neutron Star parameters (BEANSp)

Test suite
----------

Short functional test (SFT) for settle
======================================

*Note: this is only end-to-end positive test which compares one specific settle solver run with expected results and prints Passed/Failed.*

Here is how to run short functional settle test on a linux box. (tested on Ubuntu 20.04LTS)
  
1. Create/activate conda environment for beans:

.. sourcecode::
   
   # create conda env if does not exist
   conda create --name beans python==3.8
   conda activate beans
   pip install -r requirements.txt # assuming we are in the beans git repo checkout folder
   
   \... or only the "activate" line if such an envirinment already does exist.

2. compile & install settle lib (goes to ``/usr/local/lib``, requires sudo pernissions) and run the SFT in one command

.. code::

   cd settle
   make test

   This prints a bunch of lines with numbers that are the settle solver results plus binary result of one settle run (comparing with expected values) - either "PASSED" or "FAILED".


Performance test (mecnhmark)
============================


