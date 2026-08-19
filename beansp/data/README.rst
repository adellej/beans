=================
BEANSp input data
=================

This subdirectory contains data and input files used to generate results
with BEANSp

SAX J1808.4-3658
----------------

The original application of the distant ancestor of the present code,
based on 4 bursts and persistent flux measurements during the 2002
transient outburst of this 401-Hz accretion-powered millisecond pulsar.

The latest set of results were reported in `Galloway et al. (2026)`_, as part
of an attempt to better establish the use of the alpha-values in the
likelihood. Below are the relevant files:

.. _Galloway et al. (2026): http://arxiv.org/abs/26xx.yyyyy

* 1808_base.ini - configuration file for the "base" run
* 1808_bursts_newalpha.txt - updated input data file for bursts, including improved alpha measurements and peak burst fluxes
* 1808_obs.txt - input data file for persistent flux measurements

Previous analysis reported by 
`Galloway et al. (2024) <https://doi.org/10.1093/mnras/stae2422>`_ used
the original set of files and were affected by a bug in the translation of
persistent flux to accretion rate. This bug mainly affected the distance. 
The burst data file includes earlier alpha measurements, but these were
not used in the likelihoods for those runs. Below are the relevant files

* 1808_bursts.txt - input data file for bursts
* 1808_gti.txt - input data file giving GTI intervals

IGR J17498-2921
---------------

These files were used for the application to the 8-burst sample
accumulated from the 2011 outburst of this (also, coincidentally) 401-Hz
accretion-powered millisecond pulsar.

The latest set of results were also reported in `Galloway et al. (2026)`_.

* 17498_base.ini - configuration file for the new "base" run, now including the alpha values
* 17498_bursts_newalpha.txt - updated input data file for bursts, with improved alpha estimates
* 17498_obs_falanga12.txt - input data file for (bolometric) persistent flux measurements, adapted from `Falanga et al. (2012) <https://doi.org/10.1051/0004-6361/201219582>`_

Earlier, results were reported in 
`Galloway et al. (2024) <https://doi.org/10.1093/mnras/stae2422>`_ on runs
made excluding the alpha-values from the likelihood.
Since the input flux values were already bolometric, these results were
not affected by the persistent-flux to mdot bug.

* 17498_base22p.ini - configuration file for the "base" run
* 17498_bursts_8.txt - input data file for bursts

GS 1826-24
----------

Data file for three epochs of burst observations from the "Clocked
Burster", adapted from the sample provided by `Galloway et al. (2017)
<https://doi.org/10.1017/pasa.2017.12>`_.

This file is intended for use with BEANSp in "ensemble" mode

* 1826_base.ini - configuration file for the ensemble run, using the grid_mr model which in turn relies on the interpolation scheme of `Johnston et al. (2020) <http://doi.org/10.1093/mnras/staa1054>`_
* 1826_bursts.txt - example "ensemble" mode data

SRGA J144459.2–604207
---------------------

Data files for "ensemble" mode application to the 448 Hz millisecond
pulsar, with averages of up to 14 epochs of daily data from the 2024
outburst, as reported in `Galloway et al. (2026)`_.

* 144459_base.ini - configuration file for the "base" run on the 5-epoch selected data
* 144459_ensbl_day.txt - full set of daily ensemble data, with 14 epochs
* 144459_ensbl_selected.txt - reduced set of 5 epochs intended to improve agreement with the later epochs (with longer recurrence times)

As part of the most recent study, we also performed analysis on a set of
simulated data, to verify operation of the code. 

* 144459_sim_ini - configuration file for the simulation run
* 144459_ensbl_sim1.txt - simulated data with 5 epochs matching roughly the properties of SRGA J144459.2–604207


Miscellaneous files
-------------------

* mr_prior_fit.txt - used by mrprior.py to generate the mass, radius prior used in the prior_1808 function
* example_1808_pfluxfromminbar.png - screenshot demonstrating how to extract flux measurements from the `MINBAR web interface <http://burst.sci.monash.edu>`_
