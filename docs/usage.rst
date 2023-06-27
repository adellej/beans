=====
Usage
=====

To use beans, there are a few steps you need to follow ::

1. Collect all of the required data in the correct format and put it in beans/data folder.
2. Choose all of the initial conditions and initialise the class object Beans. 
3. Run the code!


Observations
------------

The code has two modes of operation: "burst train" mode, and "ensemble"
mode.

"Burst train" mode takes an uninterrupted (as much as possible) sequence of X-ray bursts, e.g. from a single outburst of a source. The measured properties of the bursts and the persistent flux covering the outburst are supplied in ascii format files, usually in the beans/data/ folder. There are 3 types of input data files, with only 2 required: persistent flux, burst parameters, and satellite gtis.

"Ensemble" mode takes multiple sets of regular burst measurements, for example as provided in `Galloway et al. (2017)`_.  The burst measurements, along with the persistent flux, are all supplied in the burst parameter file, and the other two files are not required

.. _Galloway et al. (2017): https://ui.adsabs.harvard.edu/abs/2017PASA...34...19G

**Persistent flux history:**
Required for "burst train" mode, set via the ``obsname`` parameter; set to
``None`` for "ensemble" mode runs. 
Ascii format, with columns in the following order:
start time (MJD), stop time (MJD) persistent flux measurements (in 1e-9 erg/cm^2/s), pflux error

**Burst observations:**
Required for both modes, and set via the ``burstname`` parameter; for "burst train" mode, ascii format, with columns in the following order:
time (MJD), fluence (in 1e-6 erg/cm^2/s) fluence error, alpha, alpha
error. For "ensemble" mode, the columns are 
time (MJD), fluence (in 1e-6 erg/cm^2/s) fluence error, alpha, alpha
error, (bolometric) persistent flux (1e-9 erg/cm^2/s), persistent flux error, recurrence time (hr), recurrence time error



**Satellite gtis:**
These are the satellite telescope "good time intervals" (GTIs), specifying
intervals when the telescope was actually observing the source, and so
(for example) can be used to rule out the presence of predicted bursts
within those intervals in which bursts were not observed. The GTIs are
only used for "burst train" mode, and the file should be specified via the
``gtiname`` parameter. The GTIs should be available from the raw telescope data. The file format should be a tab-separated file with 2 columns: start time of obs, stop time of obs (both in MJD).

Once you have collected the required data in the correct format and placed it in the beans/data/ folder, you can move on to initialisation.


Initialisation
--------------

Initialisation is done by instantiating a Beans object. The parameters you
need to specify are listed below. Example iniitlisation would be something like:

.. code-block:: console

    from beansp import Beans
    
    B = Beans(nwalkers=200, nsteps=100, run_id="1808/test1", obsname='1808_obs.txt', burstname='1808_bursts.txt', gtiname='1808_gti.txt', theta= (0.5, 0.015, 0.2, 2.1, 3.5, 0.108, 0.90, 0.5, 1.4, 11.2), numburstssim=3, bc=2.21, ref_ind=1, threads = 4, restart=False)

The code should display some information to the terminal that will tell you if reading in the observation data and testing the model was successful. If there are no errors here, move on to running the code. 


- **nwalkers**
``nwalkers`` is the number of walkers you want the MCMC algorithm to use. Something around 200 should be fine. If you are having convergence issues try doubling the number of walkers - check out the emcee documentation for more information.

- **nsteps**
The desired number of steps the MCMC algorithm will take. Every 100 steps the code checks the autocorrelation time for convergence and will terminate the run if things are converged. So you can set nsteps to something quite large (maybe 10000), but if things are not converging the code will take a very long time to run.

- **theta**
This sets the initial location of your walkers in parameter space.  ``theta`` includes each of the input parameters to the model:

.. code-block:: console

    theta = X, Z, Q_b, f_a, f_E, r1, r2, r3, M, R

So an example set of starting conditions would be:

.. code-block:: 

    theta = 0.5, 0.015, 0.2, 2.1, 3.5, 0.108, 0.90, 0.5, 1.4, 11.2

See parameters for a description of each of the parameters.

- **run_id**
A string identifier to label each code run you do. 
It can include the location that the chains and analysis are saved. E.g.
if I were modelling SAX J1808.4--3658 I would choose something like
``run_id = "1808/test1"``. If the package is installed as recommended, you
can run the code from within the directory in which you wish to store the
output

- **threads**
This is required because emcee runs in parallel, so needs to know how many threads (or how many cores your computer has) that it can run on. This is usually 4 for a standard computer.

- **ref_ind**
Index of the adopted reference burst, for "burst train" mode. In this mode the code simulates the burst train both forward and backward in time, so the reference burst should be in the middle of predicted burst train; don't forgot Python indexing starts at 0. This burst will not be simulated but will be used as a reference to predict the times of other bursts.

- **numburstssim**
In "burst train" mode, this is the number of bursts to simulate *in each
direction*. I.e. set to roughly half the number of bursts you want to
simulate, to cover your entire observed train. Don't forget to account for missed bursts!

- **obsname**
Path to observation data file. Should be a string, e.g.  '/Users/adelle/Documents/beans/data/1808_obs.txt'. Set to ``None`` to trigger an "ensemble" run

- **burstname**
Path to burst data file. Should be a string, e.g. '/Users/adelle/Documents/beans/data/1808_bursts.txt'

- **gtiname**
Path to GTI data file. Should be a string, e.g.
'/Users/adelle/Documents/beans/data/1808_gti.txt'. Set to ``None`` to skip
GTI checking

- **bc**
Bolometric correction to apply to the persistent flux measurements, in "burst train" mode. If they are already bolometric estimates just set this to 1.0.

- **restart**
If your run is interrrupted and you would like to restart from the save file of a previous run with the ``run_id`` set above, set this to True.  Can also be used if your max step number was not high enough and the chains did not converge before the run finished if you want to start where it finished last time. If this is a new run, set this to ``False``.


Running the Code
----------------

Once you have initialised the ``Beans`` object and ensured all the data is
available, you are ready to go. Running the code is done with the following command:

.. code-block:: console

    B.do_run()
    

If all is well you will see a progress bar appear which will give you an idea of how long the run is going to take. 

When you see ``Complete! Chains are converged`` this means the run finished, and the chains were converged. 

When you see ``Complete! WARNING max number of steps reached but chains
are not converged.`` This means the run finished but reached the maximum
number of steps ``nsteps`` without converging. 


Analysing the Results
---------------------

The output of the MCMC algorithm is saved in HDF5 format, and will be
located in whichever folder you chose when you set ``run_id``. For initial analysis of the chains you can run:

.. code-block:: console

    B.do_analysis()

And it will create a plot showing the estimated autocorrelation times
throughout the run, as well as the posterior distributions of your
parameters. 

Typically you will omit the initial "burn-in" phase and only use the
walker positions in the later part of the run; you can specify how many
steps to skip with the ``burnin`` parameter.

The model predictions at each step are saved in the "blobs" part of the sampler, which are used together with the parameter values to display the various plots below. For compatibility with the HDF5 format the model prediction dictionary must be converted to a string, and so it needs to be turned back into a dictionary item-by-item (e.g. with ``eval``) when you read in the save file. 

Several other options are possible for built-in analysis, and can be
specified via the ``options`` keyword to ``do_analysis``, which accepts a
list of strings, specifying one or more of:

``autocor``
  plot estimates of the autocorrelation times for each parameter, as a function of timestep

``chain``
  plot the first 300 iterations of the chains

``posteriors``
  show a "corner" plot giving the distirbution of the raw posteriors of the model parameters

``mrcorner``
  show a "corner" plot with just the neutron star parameters, *M*, *R*, *g* and *1+z*

``fig6``
  replicate Figure 6 from `Goodwin et al. (2019) <https://doi.org/10.1093/mnras/stz2638>`, a "corner" plot with *xi_b*, *xi_p*, *d*, *Q_b*, *Z*

``fig8``
  replicate Figure 8 from `Goodwin et al. (2019) <https://doi.org/10.1093/mnras/stz2638>`, plotting *xi_b* vs. *xi_p* and models (where available, via the `concord <https://github.com/outs1der/concord>` repository) for comparison',

``comparison``
  plot the observed and predicted burst times and fluences

You can choose to display the figures for each analysis, or save to a PDF
by specifying ``savefig=True`` in the call to ``do_analysis``.

**Checking Chain Convergence**

There are two main methods of checking the convergence and behaviour of your MCMC chains. One is the autocorrelation time, which ``emcee`` conveniently calculates for you, and the other is the acceptance fraction. Goodman and Weare (2010) provide a good discussion on what these are and why they are important. Running ``analyse.py`` will print these to the terminal for you to check. 

**Obtaining Parameter Constraints**

The model parameter posterior distributions are the most detailed
constraints on your parameters provided by the  MCMC algorithm. However,
you may wish to summarise by giving central values with uncertainties to
report for the parameters. There are a few ways this can be done; e.g. 
take the maximum likelihood value and the upper and lower limits
encompassing the desired confidence fraction, or you could take the 50th
percentile value of
the distributions. The analysis code in ``do_analysis`` does this one way,
but you should always check multiple methods and see if the results are
significantly different.

The walker position are converted to give two additional parameters,
distance *d*, burst emission anisotropy *xi_b*, and persistent emission
anisotropy *xi_p*.
The central values of these and 1 sigma
uncertainties are saved in the text file
``(run_id)_parameterconstraints_pred.txt``.

The  model predictions include the burst time, fluence, and alpha, which are stored as arrays containing an entry for each of the predicted bursts. These arrays will include as many elements as are chosen via the ``numburstssim`` parameter on initialisation.  The time array has 1 extra element than the fluence and alpha arrays, because the latter parameters do not include predictions for the reference burst (with index ``ref_ind``). 
