""" tests likelihood and prior functions for the emcee initialisation """

from beansp.mrprior import mr_prior
import numpy as np
from beansp.get_data import get_obs
import pathlib

# def test_likelihood():

#     return

# def test_lnZprior():

#     return


def test_fit_file():
    
    M = 1.4
    R = 11.2

    return mr_prior(M, R)

# def test_lnprior(theta):

#     return


# def test_lnprob(theta, x, y, yerr):

#     return
