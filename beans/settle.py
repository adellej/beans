""" Runs settle based on the input parameters and extracts the recurrence
time, fluence and alphas """

from pySettle import settler as se
import numpy as np


def settle(base, z, x_0, mdot, cfac, mass, radius):

    # initialize settle interface
    settl = se.Settle()

    # run settle:
    res = settl.full(F=base, M=mdot, X=x_0, Z=z, C=0, R=radius, Ma=mass)

    # extract results for comparison with obs
    result = np.recarray(
        (1,), dtype=[("tdel", np.float64), ("E_b", np.float64), ("alpha", np.float64)]
    )
    # assign elements
    result.tdel = res[1] * cfac
    result.E_b = res[2] * cfac
    result.alpha = res[0]
    result.mdot = mdot

    return result
