==========
beans test
==========

Bayesian Estimation of Accreting Neutron Star parameters (BEANSp)
-----------------------------------------------------------------

Test suite
==========

Short functional test (SFT) for settle
--------------------------------------

*Note: this is only end-to-end positive test which compares one specific settle solver run with expected results and prints Passed/Failed.*

Here is how to run short functional settle test on a linux box. (tested on Ubuntu 20.04LTS)

1. compile and install settle lib (goes to /usr/local/lib).

::
   cd beans/settle
   make install
   
Create/activate conda environment for beans:

::
   # create conda env if does not exist
   conda create --name beans python==3.8
   conda activate beans
   pip install -r requirements.txt
   
... or only the "activate" line if already does exist


Performance test (mecnhmark)
----------------------------


_Cookiecutter: https://github.com/audreyr/cookiecutter
_`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
