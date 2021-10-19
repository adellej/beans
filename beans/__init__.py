# -*- coding: utf-8 -*-

"""Top-level package for beans."""

__author__ = """Adelle Goodwin"""
__email__ = 'adelle.goodwin@monash.edu'
__version__ = '0.1.0'


from beans import Beans

B = Beans(ndim=10, nwalkers=450, nsteps=2000, run_id="",obsname='none', burstname='../data/1826_bursts.txt', gtiname='../data/1808_gti.txt', theta= (0.52, 0.025, 3.0, 1.5, 0.5, 0.42, 0.4, 0.4, 1.8, 6.0), numburstssim=3, numburstsobs=3, bc=1.5, ref_ind=1, gti_checking=0, restart=False,train=0)

B.do_run()
