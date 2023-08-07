===============
Building BEANSp
===============

Build and installation from this github repository
--------------------------------------------------

#. Clone the beans repository

   .. code-block:: console
    
      git clone https://github.com/adellej/beans
      cd beans
   

#. Create and activate a clean conda environment

   The example is for python 3.8, but should work for any version 3.6 to 3.11 as well.

   .. code-block:: console
    
      # optional: remove existing environment if needed - to start from scratch
      conda remove -n beans-3.8 --all
      # create blank conda environment
      conda create --name beans-3.8 python==3.8.*
      conda activate beans-3.8

      
#. Install/upgarde pip, build and local install

   .. code-block:: console
  
      python3 -m pip install --upgrade pip
      python3 -m pip install --upgrade build

   .. code-block:: console
  
      # test build & local install
      # The "-e" install does not seem to be reliable for re-install on Linux
      #       - keeps pulling some old build from somewhere middlewhere.
      #         python -m pip install -e .*
      # This is more reliable:
      python3 -m build
      python3 -m pip install .

   .. ::
   
   *Note: when workinng on the code, in case of doubts that recent changes got propagated, uninstall & purge the installed module _before_* ``pip install`` *to ensure the installed version has all the recent modifications.*

   .. code-block:: console
     
      python3 -m pip -v uninstall beansp
      python3 -m pip -v cache purge

   After this, in that enviroment, beansp just works from every directorty, providing the conda environment is activated.
   Imports like:

   .. code-block:: python
   
      from beansp.beans import Beans 

   (See `test_sft_beans.py <tests/test_sft_beans.py>`_.)


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
