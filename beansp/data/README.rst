======
BEANSp input data
======

This subdirectory contains data and input files used to generate results
with BEANSp

SAX J1808.4-3658
----------------

These files have been used for the original application, with 4 bursts
observed during the 2002 transient outburst of this 401-Hz
accretion-powered millisecond pulsar.

The latest set of results were reported in `Galloway et al. (2024) <https://doi.org/10.1093/mnras/stae2422>`_

* 1808_base.ini - configuration file for the "base" run
* 1808_bursts.txt - input data file for bursts
* 1808_gti.txt - input data file giving GTI intervals
* 1808_obs.txt - input data file for persistent flux measurements

IGR J17498-2921
---------------

These files were used for the application to the 8-burst sample
accumulated from the 2011 outburst of this (also, coincidentally) 401-Hz
accretion-powered millisecond pulsar.

The latest set of results were reported in `Galloway et al. (2024) <https://doi.org/10.1093/mnras/stae2422>`_

* 17498_base22p.ini - configuration file for the "base" run
* 17498_bursts_8.txt - input data file for bursts
* 17498_obs_falanga12.txt - input data file for (bolometric) persistent flux measurements, adapted from `Falanga et al. (2012) <https://doi.org/10.1051/0004-6361/201219582>`_

GS 1826-24
----------

Data file for three epochs of burst observations from the "Clocked
Burster", adapted from the sample provided by `Galloway et al. (2017)
<https://doi.org/10.1017/pasa.2017.12>`_.

This file is intended for use with BEANSp in "ensemble" mode

1826_bursts.txt

Miscellaneous files
-------------------

* mr_prior_fit.txt - used by mrprior.py to generate the mass, radius prior used in the prior_1808 function
* example_1808_pfluxfromminbar.png - screenshot demonstrating how to extract flux measurements from the `MINBAR web interface <http://burst.sci.monash.edu>`_
