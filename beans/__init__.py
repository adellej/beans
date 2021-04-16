# -*- coding: utf-8 -*-

"""Top-level package for beans."""

__author__ = """Adelle Goodwin"""
__email__ = 'adelle.goodwin@monash.edu'
__version__ = '0.1.0'


from beans import Beans
train = 1

from beans import Beans

B = Beans(ndim=10, nwalkers=200, nsteps=100, run_id="1808/test1", obsname='1808_obs.txt', burstname='1808_bursts.txt', gtiname='1808_gti.txt', theta= (0.5, 0.015, 0.2, 2.1, 3.5, 0.108, 0.90, 0.5, 1.4, 11.2), numburstssim=3, numburstsobs=4, bc=2.21, ref_ind=1, gti_checking=0, threads = 4, restart=False)

B.do_run()
