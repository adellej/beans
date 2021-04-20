# -*- coding: utf-8 -*-

"""Top-level package for beans."""

__author__ = """Adelle Goodwin"""
__email__ = 'adelle.goodwin@monash.edu'
__version__ = '0.1.0'


from beans import Beans

B = Beans(ndim=10, nwalkers=450, nsteps=2500, run_id="1808/test1",obsname='none', burstname='../data/1826_bursts.txt', gtiname='../data/1808_gti.txt', theta= (0.5, 0.02, 1.0, 2.1, 3.5, 0.5, 0.50, 0.5, 1.7, 11.2), numburstssim=1, numburstsobs=3, bc=2.21, ref_ind=1, gti_checking=0, restart=False,train=0)


B.do_run()
