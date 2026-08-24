===============
Building BEANSp
===============

Build and installation from this github repository
--------------------------------------------------

#. Clone the beans repository

   .. code-block:: console
    
      git clone https://github.com/adellej/beans
      cd beans
   

#. [optional] Create and activate a clean conda environment

   The example below will create and activate an environment named `beans` with all requirements. Tested with Python 3.13, but beans should work with Python 3.9 onwards.

   The first command below will remove an existing environment if needed, e.g. to start from scratch

   .. code-block:: console
    
      conda remove -n beans --all
      conda env create -f environment.yml
      conda activate beans

#. Install/upgrade pip, build and local install

   .. code-block:: console
  
      python3 -m pip install --upgrade pip
      python3 -m pip install --upgrade build

      python3 -m build
      python3 -m pip install .

   *Note: when working on the code, in case of doubts that recent changes got propagated, uninstall & purge the installed module _before_* ``pip install`` *to ensure the installed version has all the recent modifications.*

   .. code-block:: console
     
      python3 -m pip -v uninstall beansp
      python3 -m pip -v cache purge

After completing these steps, once the enviroment is activated, beansp
should just work from every directory.  Import as follows:

.. code-block:: python
   
      from beansp import Beans,beans

..
    removed from the code block above:
      # test build & local install
      # The "-e" install does not seem to be reliable for re-install on Linux
      #       - keeps pulling some old build from somewhere middlewhere.
      #         python -m pip install -e .*
      # This is more reliable:


Testing
-------

Once you have compiled settle we recommend you run the test suite to check you have all the required dependencies and the code is operating as expected. To do this navigate to the top-level directory and type:

.. code-block:: console

    pytest

Run short functional test (SFT) manually
----------------------------------------

.. code-block:: console

   cd tests
   python ./test_sft_beans.py
 


If the tests all pass then you are good to go!
