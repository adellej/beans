=====
Usage
=====

To use beans, there are a few steps you need to follow ::

1. Collect all of the required data in the correct format and put it in beans/data folder.
2. Setup all of the initial conditions in initialise.py
3. Run the code!

============
Observations
============

This code works best with an uninterrupted train of X-ray bursts from a single outburst of a source. The observations need to be in ascii format and put in the beans/data/ folder. There are 3 types of observations, with only 2 required. These are persistent, burst, and satellite gtis.

Persistent observations:

Burst observations:

Satellite gtis:

Once you have collected the required data in the correct format and placed it in the beans/data/ folder, you can move on to initialisation.

==============
Initialisation
==============

All initialisation is done in initialise.py. This is the only code that you need to edit (unless you want to change/edit things). The paramters you need to enter are listed below.

- **ndim, nwalkers**
ndim is the dimension of your parameter space (will be 10 unless you add extra parameters to theta). nwalkers is the number of walkers you want the MCMC algorithm to use. Something around 200 should be fine. If you are having convergence issues try doubling the number of walkers - check out the emcee documentation for more information.

- **nsteps**
This is the maximum number of steps the MCMC algorithm will take. Every 100 steps the code checks the autocorrelation time for convergence and will terminate the run if things are converged. So you can set nsteps to something quite large (maybe 10000), but if things are not converging the code will take a very long time to run.

- **theta**
This sets the initial location of your walkers in parameter space. Theta is each of the parameters we care about:

.. code-block:: console

    theta = X, Z, Q_b, f_a, f_E, r1, r2, r3

So an example set of starting conditions would be:
.. code-block:: 

    theta = 0.5, 0.015, 0.2, 2.1, 3.5, 0.108, 0.90, 0.5, 1.4, 11.2

See parameters for a description of each of the parameters.

- **run_id**
This is a unique identifier for each code run you do. Should be a string, and it will set the location that the chains and analysis are saved. E.g. if i were modelling SAX J1808.4--3658 I would choose something like run_id = "1808/test1". I recommend using a folder for each different source you are modelling. 

- **threads**
This is required because emcee runs in parallel, so needs to know how many threads (or how many cores your computer has) that it can run on. This is usually 4 for a standard computer.

- **numburstssim**
This needs to be an integer value of half the number of bursts you want to simulate. I.e. simulate this many from the reference burst in either direction. Don't forget to account for missed bursts!

- **numburstsobs**
Number of observed bursts in your dataset

- **ref_ind**
Index of the reference burst (should be middle of predicted burst train - don't forgot python indexing starts at 0). This burst will not be simulated but will be used as a reference to predict the other bursts.

- **gti_checking**
This is an option to turn on gti checking. 1 for on, 0 for off. If this is on, the code will check that each modelled burst train predicts bursts that were not observed ones to fall in satellite observing gaps. 

- **obsname**
Path to observation data file. Should be a string, e.g. '/Users/adelle/Documents/beans/data/1808_obs.txt'. 

- **burstname**
Path to burst data file. Should be a string, e.g. '/Users/adelle/Documents/beans/data/1808_bursts.txt'

- **gtiname**
Path to gti data file. Should be a string, e.g. '/Users/adelle/Documents/beans/data/1808_gti.txt'

- **bc**
Bolometric correction to apply to the persistent flux measurements. If they are already bolometric fluxes just set this to 1.0.

- **restart**
If your run crashes and you would like to restart from the save file of a previous run with the run_id set above, set this to True. Can also be used if your max step number was not high enough and the chains did not converge before the run finished if you want to start where it finished last time. If this is a new run, set this to False.

================
Running the Code
================

Please note that the code can take a long time (~week) to run, depending on the number of bursts in the burst train, and the number of steps you choose to use. So I recommend running it on a desktop you know is not going switch off and using terminal software such as tmux or similar. 

Once you have filled out the required parameters in initialise.py and put all of the required data files in beans/data/, you are ready to go. To run the code just type:

.. code-block:: console

    $ python beans.py

This will print some text to the terminal and if all is well you will see a progress bar appear which will give you an idea of how long the run is going to take. When you see "Complete! Chains are converged" this means the run finished, and the chains were converged. When you see "Complete! WARNING max number of steps reached but chains are not converged." This means the run finished but reached the maximum number of steps (nsteps) without converging. 

=====================
Analysing the Results
=====================

The output of the MCMC algorithm is saved in hdf5 format, and will be located in whichever folder you chose when you set **run_id**. For initial analysis of the chains you can run:

.. code-block:: console
    $ python analyse.py

And it will create a plot of the posterior distributions of your parameters. 
