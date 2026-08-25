=======
History
=======

3.0.0 (2026-08-25)
------------------

This release includes a host of small improvements and bug fixes, as described below, and is intended to accompany the paper "Constraining neutron-star properties with ensembles of thermonuclear bursts: application to SRGA J144459.2−604207" by D. K. Galloway et al. (submitted to PASA; available at http://arxiv.org/abs/2608.21196).

* fixed a bug with a spurious bolometric correction term in flux_to_mdot, which affected all runs from 2023 September through to 2026 June (versions 2.2.0-2.68.0) for which the bolometric correction was >1. This bug unfortunately affected the SAX J1808.4-3658 runs in the 2024 paper, but have now been revised in the 2026 paper.
* fixed a bug in sim_data that would mix up the model-generated alpha values; also now allows generation of simulated ensemble mode data (including from a burst train)
* added a new burst train-generating routine, punkt_train, which allows gaps and "independent" burst trains (i.e. where the recurrence time for each burst depends on the previous observed burst (not the previous simulated burst)
* added the routine grid_interp and accompanying prior, prior_grid, for use with the multiepoch_mcmc model grids of Johnston et al. 2020 (see https://doi.org/10.1093/mnras/staa1054 ). Model is now stored with the beans object and will be saved in the .ini file for checking on restart
* added preliminary bilby implementation with run_bilby.py
* added compare method to plot multiple posteriors from different runs simultaneously
* replaced the systematic errors on fluence and alpha, with an optional parameter on time (f_t)
* modified the burst_table method to also generate more accurate "observed" alpha values incorporating the persistent flux model
* added a more extensive set of data files for the various sources on which the code has been applied, including the new estimates of the alpha values
* added rst (reStructuredText) option for write_param_uncert method, to allow parameter tables to be included in documentation
* added archive method to save (pickle) the Beans objects, and a restore method to recover them; permits deleting the .h5 files to save space
* extended the prune method which now permits adding additional parameters (e.g. M, R if the existing run was performed with those parameters fixed) and walkers (if you want to expand the run, e.g. for a "production" run for publication). Resulting arrays or pickle files can be used to continue/restart runs, via the pos option to do_run
* updated plot method options with "page" and "column" layout, intended for inclusion into papers with suitable font sizes etc. Can now also specify an alternate save filename by passing to savefig
* added the option to also include a likelihood component based on the peak burst fluxes; can list the peak burst fluxes in the burst data file, along with the PRE status; then likelihood is based on the comparison of the peak fluxes with the model-generated (isotropic) Eddington flux
* lots of small changes and updates to ensure correct build and running with Python 3.13; added environment.yml
* upgraded to the new version of ChainConsumer and updated all the function calls etc.

2.62.0 (2025-01-11)
-------------------

This version has been updated to work with Python 3.13, including using
the newest version of ChainConsumer for the plots. 

2.25.0 (2024-03-16)
-------------------

This release includes a host of small improvements, as described below, and is intended to accompany the paper "Inferring system parameters from the bursts of the accretion-powered pulsar IGR J17498–2921" by D. K. Galloway et al. (in prep).

* the model parameter ratios have been replaced with physical parameters including distance, and burst and persistent flux anisotropy, explicitly
* can now specify spline interpolation for the fluxes, with an optional smoothing parameter
* depending upon the number of parameters provided in the initialisation stage, the comparisons can be performed with "canonical" values for neutron star mass and radius, and no systematic errors
* the alphas and fluences can also be excluded from the comparison
* the burst model matching has been improved to work more reliably, and the analysis scripts are more robust in cases of inhomogeneous model predictions
* can now choose explicitly the starting points for the walkers, and a prune function is provided to fine-tune the walker ensemble.
* prior to running the chains, the code checks if the provided positions correspond to valid models
* the burst_table method will generate a table including all burst parameters, including alphas inferred from the model predictions

1.0.0 (2023-05-18)
------------------

This release presents a substantial update of the code (now known as beansp) with the help of software engineers at ADACS (https://adacs.org.au)

The settle routine has been significantly updated, tested to build and run on both linux and mac platforms, and packaged separately, now available on pypi.

The beans package has been updated extensively, with bug fixes, more reliable operations, and additional comments. New installation and testing procedures have been established and the package is also available on pypi

0.1.0 (2019-09-19)
------------------

First release on PyPI.
