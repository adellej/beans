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

Activate conda environment based on python 3.8 (or later) and requirements.txt in the root of beans repo

    | # remove existing environment if needed - to start from scratch
    | conda remove -n settle-3.8 --all
    | # create blank conda environment (has numpy in it by defeult now with python 3.8)
    | conda create --name settle-3.8 python==3.8
    | conda activate settle-3.8

then do

    | cd settle
    | # upgrade pip and build
    | python3 -m pip install --upgrade pip
    | python3 -m pip install --upgrade build


    | # test build & local install
    | # The "-e" install does not seem to be reliable for re-install - keeps pulling some old build from some where.
    | # Do not use:        python -m pip --verbose install -e .
    | # This is more ereliable:
    | python3 -m build --sdist
    | python3 -m pip install .

After this, in that enviroment, pySettle just works from every directorty, providing the conda environment is activated.
Imports like:

    | from pySettle import settler as se

(See beans/settle/test/test_settle_sft.py)

Run short functional test (SFT) manually
----------------------------------------

    | cd test
    | python ./test_settle_sft.py
    

