======
settle
======

Settling solver - the beans edition
-----------------------------------------------------------------

* Forked from settle project by Andrew Cumming
  https://github.com/andrewcumming/settle
* Repo: https://github.com/adellej/beans/settle


Features
--------

This code computes ignition conditions for Type I X-ray bursts using a multi-zone model of the accreting layer (including hot CNO hydrogen burning, but not helium burning), but a one-zone ignition criterion. For more details, see Cumming & Bildsten (2000).

Credits
-------

Rotational Evolution During Type I X-Ray Bursts, Andrew Cumming, Lars Bildsten (2000) - https://arxiv.org/abs/astro-ph/0004347

Package installation and usage
------------------------------

activate conda environment based on python 3.8 and requirements.txt in the root of beans repo, then do

    cd settle

    # upgrade pip and build
    python3 -m pip install --upgrade pip
    python3 -m pip install --upgrade build

    # test build & local install
    python -m pip --verbose install -e .

After this, in that enviroment, pySettle just works from everywhere.
imports like

    from pySettle import settler as se

(See beans/settle/test/try.py)

