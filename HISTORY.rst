=======
History
=======

2.62.0 (2025-01-11)
------------------

This version has been updated to work with Python 3.13, including using
the newest version of ChainConsumer for the plots. 

2.25.0 (2024-03-16)
------------------

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
