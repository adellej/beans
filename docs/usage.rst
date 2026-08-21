===========
Basic usage
===========

To use beans, there are a few steps you need to follow:

1. Collect all of the required observational data in the correct format
and put it in ``beans/data`` folder.
2. Choose all of the initial conditions and initialise the class object :class:`beansp.Beans`.
3. Run the code!
4. Analyse the results


1. Observational data
---------------------

The code has two principal modes of operation: "burst train" mode, and
"ensemble" mode.

"Burst train" mode takes an uninterrupted (as much as possible) sequence
of X-ray bursts, e.g. from a single outburst of a source. The measured
properties of the bursts and the persistent flux covering the outburst are
supplied in ASCII format files, usually in the ``beans/data`` folder. There
are 3 types of input data files:
burst properties, persistent flux history, and (optional) satellite gtis.

"Ensemble" mode takes multiple sets of regular burst measurements, for
example as provided in `Galloway et al. (2017)`_.  The burst measurements,
along with the persistent flux, are all supplied in the burst properties
file, and no other files are required.

.. _Galloway et al. (2017): https://ui.adsabs.harvard.edu/abs/2017PASA...34...19G

**Burst observations:**
Required for both modes, and set via the ``burstname`` parameter. For
"burst train" mode, ASCII format, tab-separated columns in the following order:
time (MJD), bolometric burst fluence & error (in units of
:math:`10^{-6}\ \rm{erg\,cm^{-2}}`),
alpha, alpha error. An example distributed with the code is
``data/1808_bursts.txt``.

.. note:: The burst observations file for "train" mode can optionally include (before the alpha columns) the peak
    burst flux & error (in 
    units of :math:`10^{-9}\ \rm{erg\,cm^{-2}\,s^{-1}}`), and a flag
    indicating the presence of photospheric radius-expansion (1=non-PRE,
    2=PRE); these columns
    are for use with the ``pflux=True`` option. 

For "ensemble" mode, the columns are
time (MJD, epoch/arbitrary), fluence & error (in :math:`10^{-6}\ \rm{erg\,cm^{-2}}`),
alpha & error, (bolometric) persistent flux & error
(:math:`10^{-9}\ \rm{erg\,cm^{-2}\,s^{-1}}`),
recurrence time & error (hr).
An example distributed with the code is
``data/1826_bursts.txt``.

**Persistent flux history:**
Required for "burst train" mode, set via the ``obsname`` parameter; set to
``None`` for "ensemble" mode runs.
ASCII format, with tab-separated columns in the following order:
start time (MJD), stop time (MJD) persistent flux & error (in 
:math:`10^{-9}\ \rm{erg\,cm^{-2}\,s^{-1}}`).
An example distributed with the code is
``data/1808_obs.txt``.

.. note:: It is important that the persistent flux measurements fully cover the burst train.


Once you have collected the required data in the correct format
you can move on to initialisation.


2. Initialisation
-----------------

Initialisation is done by instantiating a :class:`beansp.Beans` object (a "bean", why
not). The parameters you might normally
need to specify are listed below.

Example initialisation would be something like:

.. code-block:: python

    from beansp import Beans

    B = Beans(nwalkers=200, nsteps=100, run_id="test1", 
        obsname='1808_obs.txt', burstname='1808_bursts.txt', 
        theta= (0.58, 0.013, 0.4, 3.5, 1.0, 1.0, 1.5, 11.8, 1.0), 
        numburstssim=3, bc=2.21, ref_ind=1, threads = 4, restart=False)

The above example will create an object labeled ``test1`` for a "train"
mode run (based on the presence of both ``obsname`` and ``burstname``
files) and read in
observations and bursts from the nominated files, in the current
directory.
The code should display some information to the terminal that will tell you if reading in the observation data and testing the model was successful. 

.. note:: At the commencement of a run the code will save the settings in a file
    ``<run_id>.ini`` which serves as a record, along with the ``.h5`` file
    that stores the chains and model evaluations. If you're restarting or
    re-running the same experiment you can use the ``config_file`` option to
    replicate:

.. code-block:: python

    B = Beans(config_file="test1.ini", prior=beans.prior_1808,corr=beans.corr_goodwin19)
    B.restart = True
    B.nsteps = 1000
 
In the example given above we're reading in all the parameters from the
previous ``test1`` run, but updating the number of steps, to go longer this
time (for example).

Note also the specification of the ``prior`` and ``corr`` functions; these
settings replicate the baseline run for SAX J1808.4-3658.
Each of the initialisation parameters are described in more detail in
:ref:`parameter-description`.

A good way to get started with the code is to use one of the supplied
``.ini`` files in the ``data`` subdirectory to initialise a run and work
through a simulation and analysis cycle.

The MCMC walkers will initially be distributed around the starting point given
by the ``theta`` parameter, which includes each of the input parameters to
the model:

.. code-block:: python

    theta = X, Z, Q_b, d, xi_b, xi_p [, M, R [, f_t] ]

(the square brackets indicate that ``M, R`` and ``f_t`` are optional):

| ``X`` - Hydrogen mass fraction
| ``Z`` - CNO metallicity
| ``Q_b`` - "base" flux (MeV/nucleon)
| ``d`` - source distance (kpc)
| ``xi_b`` - anisotropy term for burst emission
| ``xi_p`` - anisotropy term for persistent emission
| ``M`` - neutron star mass (:math:`M_\odot`)
| ``R`` - neutron star radius (km)
| ``f_t`` - systematic error on burst times

So an example set of starting conditions would be:

.. code-block:: python

    theta = 0.36, 0.016, 0.8, 2.4, 1.0, 1.4, 2.1, 12.2

The ``f_t`` (systematic error on the times) in this case is not included.

..
     see parameters for a description of each of the parameters.

You can visualise the predictions with the starting conditions using 
the :meth:`beansp.Beans.plot` method.

If there are no errors or other issues here, move on to running the code.

3. Running the Code
-------------------

Once you have initialised the :class:`beansp.Beans` object and ensured all the data is
available, you are ready to go. Running the code is done with the following command:

.. code-block:: python

    B.do_run()


If all is well you will see a progress bar appear which will give you an idea of how long the run is going to take.
The output of the MCMC algorithm is saved in HDF5 format, and will be
located in whichever folder you chose when you set ``run_id``.

When you see ``Complete! Chains are converged`` this means the run finished, and the chains were converged.

When you see ``Complete! WARNING max number of steps reached but chains
are not converged.`` This means the run finished but reached the maximum
number of steps ``nsteps`` without converging.

.. _analysis:

4. Analysing the Results
------------------------

Once the run is complete and/or you've reloaded the Beans object from an
``.ini`` file, you can examine the progress as follows:

.. code-block:: python

    B.do_analysis('chain') # or
    B.chain()

which will create a "trace" plot showing the evolution of the walkers
throughout the run.

Typically you will omit the initial "burn-in" phase and only use the
walker positions in the later part of the run; you can specify how many
steps to skip with the ``burnin`` parameter. Commonly we include (say) the
last 1000 steps, indicated by a negative value, i.e. ``burnin=-1000``.

You can visualise the posterior distributions with the ``posteriors``
option:

.. code-block:: python

    B.do_analysis('posteriors') # or
    B.posteriors()

For each of the analysis options you can add ``savefig=True`` or
``savefig=<file name>`` to save the plot to a file. For the posteriors
option, if you know the input parameters, you can indicate them via the
``truths`` parameter.

The model predictions at each step are saved in the "blobs" part of the sampler, which are used together with the parameter values to display the various plots below. 

..
    For compatibility with the HDF5 format the model prediction dictionary must be converted to a string, and so it needs to be turned back into a dictionary item-by-item (e.g. with ``eval``) when you read in the save file.

Several other options are possible for built-in analysis, and can be
specified via the ``options`` keyword to :meth:`beansp.Beans.do_analysis`,
which accepts a list of strings, specifying one or more of:

``autocor``
  plot estimates of the autocorrelation times for each parameter, as a function of timestep

``chain``
  plot the walker positions

``posteriors``
  show a "corner" plot giving the distirbution of the raw posteriors of the model parameters

``mrcorner``
  show a "corner" plot with just the neutron star parameters, *M*, *R*, *g* and *1+z*

``comparison``
  plot the observed and predicted burst times and fluences

``fig6``
  replicate Figure 6 from `Goodwin et al. (2019) <https://doi.org/10.1093/mnras/stz2638>`_, a "corner" plot with *xi_b*, *xi_p*, *d*, *Q_b*, *Z*

``fig8``
  replicate Figure 8 from `Goodwin et al. (2019) <https://doi.org/10.1093/mnras/stz2638>`_, plotting *xi_b* vs. *xi_p* and models (where available, via the `concord <https://github.com/outs1der/concord>`_ repository) for comparison',

You can choose to display the figures for each analysis, or skip with
``show=False`` e.g. if you're running in a window which doesn't allow
graphical display

**Obtaining Parameter Constraints**

The model parameter posterior distributions are the most detailed
constraints on your parameters provided by the  MCMC algorithm. However,
you may wish to summarise by giving percentile values with uncertainties to
report for the parameters. There are a few ways this can be done; e.g.
take the maximum likelihood value and the upper and lower limits
encompassing the desired confidence fraction, or you could take the 50th
percentile value of
the distributions. The analysis code in :meth:`beansp.Beans.do_analysis`
does this one way, but you should always check multiple methods and see if
the results are significantly different.

Call :meth:`beansp.Beans.write_param_uncert` to save the 
central values of these and 1 sigma
uncertainties in the text file
``(run_id)_parameterconstraints_pred.txt``.
This method will also generate a LaTeX table skeleton that can be pasted
into a paper (for example).

The  model predictions include the burst time, fluence, and alpha, which are stored as arrays containing an entry for each of the predicted bursts. These arrays will include as many elements as are chosen via the ``numburstssim`` parameter on initialisation.  The time array has 1 extra element than the fluence and alpha arrays, because the latter parameters do not include predictions for the reference burst (with index ``ref_ind``).
