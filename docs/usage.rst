==========
Parameters
==========

There are 10 important parameters that the code will predict, and 3 important parameters that are observed (and you give the code). These are:

Observed:

- **time**
The observed times of each burst in your burst train.

- **fluence**
The observed fluence of each burst in your burst train. Fluence is the integrated burst flux.

- **alpha**
The observed alpha of each burst in your burst train. Alpha is the ratio between fluence and the persistent flux (i.e. nuclear to gravitational energy).

Predicted:

These are the *theta* parameters. 

- **X** (or x_0)
Hydrogen mass fraction of the accreted fuel.

- **Z**
CNO metallicity of the accreted fuel.

- **Q_b**
Base heating (MeV/nucleon)

- **f_a** 
Nuisance parameter in likelihood.

- **f_E**
Nuisance parameter in likelihood.

- **r1** 
Scaling factor.

- **r2** 
Scaling factor between observed and predicted alphas.

- **r3**
Scaling factor between observed and predicted fluences. 

For equations describing the scaling factors see Goodwin et al. (2019) Equations 12, 13, and 14 - https://arxiv.org/pdf/1907.00996 


=====
Usage
=====

To use beans, there are a few steps you need to follow ::

1. Collect all of the required data in the correct format and put it in beans/data folder.
2. Choose all of the initial conditions and initialise the class object Beans. 
3. Run the code!


Observations
------------

This code works best with an uninterrupted train of X-ray bursts from a single outburst of a source. The observations need to be in ascii format and put in the beans/data/ folder. There are 3 types of observations, with only 2 required. These are persistent, burst, and satellite gtis.

**Persistent observations:**
Ascii format, with columns in the following order:
start time (MJD), stop time (MJD) persistent flux measurements (in 3-25keV 1e-9 erg/cm^2/s), pflux error

**Burst observations:**
Ascii format, with columns in the following order:
time (MJD), fluence (in 1e-9 erg/cm^2/s) fluence error, alpha, alpha error

**Satellite gtis:**
These are the satellite telescope "good time intervals" (gtis). This is information on when the telescope was actually observing the source, and for how long, and when it was looking at other things, or was at a point in its orbit where it could not observe. The gtis should be available from the raw telescope data. The file format should be a tab-separated file with 2 columns: start time of obs, stop time of obs (both in MJD).

Once you have collected the required data in the correct format and placed it in the beans/data/ folder, you can move on to initialisation.


Initialisation
--------------

All initialisation is done by calling the class object Beans. The parameters you need to enter are listed below. Example iniitlisation would be something like:

.. code-block:: console

    from beans import Beans
    
    B = Beans(ndim=10, nwalkers=200, nsteps=100, run_id="1808/test1", obsname='1808_obs.txt', burstname='1808_bursts.txt', gtiname='1808_gti.txt', theta= (0.5, 0.015, 0.2, 2.1, 3.5, 0.108, 0.90, 0.5, 1.4, 11.2), numburstssim=3, numburstsobs=4, bc=2.21, ref_ind=1, gti_checking=0, threads = 4, restart=False)

This should print some information to the terminal that will tell you if reading in the observation data and testing the model was successful. If there are no errors here, move on to running the code. 


- **ndim, nwalkers**
ndim is the dimension of your parameter space (will be 10 unless you add extra parameters to theta). nwalkers is the number of walkers you want the MCMC algorithm to use. Something around 200 should be fine. If you are having convergence issues try doubling the number of walkers - check out the emcee documentation for more information.

- **nsteps**
This is the maximum number of steps the MCMC algorithm will take. Every 100 steps the code checks the autocorrelation time for convergence and will terminate the run if things are converged. So you can set nsteps to something quite large (maybe 10000), but if things are not converging the code will take a very long time to run.

- **theta**
This sets the initial location of your walkers in parameter space. Theta is each of the parameters we care about:

.. code-block:: console

    theta = X, Z, Q_b, f_a, f_E, r1, r2, r3, M, R

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


Running the Code
----------------

Please note that the code can take a long time (~week) to run, depending on the number of bursts in the burst train, and the number of steps you choose to use. So I recommend running it on a desktop you know is not going to switch off and using terminal software such as tmux or similar. 

Once you have filled out the required parameters in initialise.py and put all of the required data files in beans/data/, you are ready to go. Running the code is done with the following command:

.. code-block:: console

    from beans import Beans

    B = Beans(ndim=10, nwalkers=200, nsteps=100, run_id="1808/test1", obsname='1808_obs.txt', burstname='1808_bursts.txt', gtiname='1808_gti.txt', theta= (0.5, 0.015, 0.2, 2.1, 3.5, 0.108, 0.90, 0.5, 1.4, 11.2), numburstssim=3, numburstsobs=4, bc=2.21, ref_ind=1, gti_checking=0, threads = 4, restart=False)

    B.do_run()
    

This will print some text to the terminal and if all is well you will see a progress bar appear which will give you an idea of how long the run is going to take. When you see "Complete! Chains are converged" this means the run finished, and the chains were converged. When you see "Complete! WARNING max number of steps reached but chains are not converged." This means the run finished but reached the maximum number of steps (nsteps) without converging. 


Analysing the Results
---------------------

The output of the MCMC algorithm is saved in hdf5 format, and will be located in whichever folder you chose when you set **run_id**. For initial analysis of the chains you can run:

.. code-block:: console

    from beans import Beans

    B = Beans(ndim=10, nwalkers=200, nsteps=100, run_id="1808/test1", obsname='1808_obs.txt', burstname='1808_bursts.txt', gtiname='1808_gti.txt', theta= (0.5, 0.015, 0.2, 2.1, 3.5, 0.108, 0.90, 0.5, 1.4, 11.2), numburstssim=3, numburstsobs=4, bc=2.21, ref_ind=1, gti_checking=0, threads = 4, restart=False)

    B.do_analysis()

And it will create a plot of the posterior distributions of your parameters. Further analysis is up to you. 

The interesting model information is saved in the "blobs" part of the sampler. This is where the parameters for each model run that was executed by emcee are saved (the output of the generate_burst_train routine). Unfortunately to save in HDF5 format this dictionary had to be converted to a string, so it needs to be turned back into a dictionary when you read in the save file. Have a look in the function do_analysis for an example on how to do this. 

**Checking Chain Convergence**

There are 2 main methods of checking the convergence and behaviour of your MCMC chains. One is the autocorrelation time, which emcee conveniently calculates for you, and the other is the acceptance fraction. Goodman and Weare (2010) provide a good discussion on what these are and why they are important. Running analyse.py will print these to the terminal for you to check. 

**Obtaining Parameter Constraints**

The posterior distributions are the true constraints on your parameters that MCMC gives you. However, you may wish to obtain numbers with uncertainties to report for the parameters. There are a few ways this can be done, you could choose to take the maximum likelihood value, or you could take the middle value of the distributions. The analysis code in do_analysis does this one way, but you should always check multiple methods and see if the results are significantly different. The "predicted" parameters are Xpred, Zpred, basepred, dpred, cosipred, xippred, xibpred, masspred, radiuspred, gravitypred, redshiftpred, and the central values of these and 1 sigma uncertainties are saved in the text file (run_id)_parameterconstraints_pred.txt. The "observed" parameters are time, fluence, and alpha. These are arrays that contain an entry for each of the predicted bursts. These will be as long as the numburstssim you chose in the initialisation. The time array has 1 extra element than the ebs and alphas because ebs and alphas do not include predictions for the reference burst (with index tref). 