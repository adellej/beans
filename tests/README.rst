===========
beans tests
===========

Bayesian Estimation of Accreting Neutron Star parameters (BEANSp)

Performance test (benchmark) for settle
=======================================

Jupyter notebook test runs settle on 60 lines of input data from a file. Repeats the run a number of times in a loop. It prints duration of the solver run for each step, overall and average. The "settle sum time" means accumulated time spent calling settle solver, stripped of the main loop python code overhead, but including the overhead of the python wrapper around the settle C/C++ library interface function. Reference times are TBD, nothing to compare with at the moment.

*Note: execution times depend on overall system load.*

#. Compile & install settle lib (goes to ``/usr/local/lib``, requires sudo pernissions)

   .. code-block::

      cd settle
      make install
   
#. Use the same conda environment as for beans, just add jupyter

   .. code-block::

      conda activate beans
      pip install jupyter

#. Open and run the jupyter motebook

   .. code-block::

      cd ../tests/benchmark
      jupyter notebook beans_benchmark.ipynb


