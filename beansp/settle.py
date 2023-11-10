from pySettle import settler as se
import numpy as np


def settle(base, z, x_0, mdot, mass, radius, corr=None):
    """
    Runs settle based on the input parameters and extracts the recurrence
    time, fluence and alphas
    Can also do various rescaling here, e.g. via the cfac parameter

    :param base: base flux [MeV/nucleon]
    :param z: accreted CNO metallicity
    :param x_0: accreted H-fraction
    :param mdot: accretion rate as a fraction of Eddington (as defined by
      settle)
    :param mass: neutron star mass [M_sun]
    :param radius: neutron star radius [km]
    :param corr: correction function to modify the settle output (replaces the old cfac option)

    :returns: 3-element array with recurrence time [hr], burst energy
      [1e39 erg], and alpha, all in the observer frame
    """

    # initialize settle interface
    settl = se.Settle()

    # run settle, returning a 3-element tuple with alpha, trec [hr], and
    # burst energy [1e39 erg], all values in the observer frame

    res = settl.full(F=base, M=mdot, X=x_0, Z=z, C=0, R=radius, Ma=mass)

    if corr is not None:
        res = corr(res, F=base, M=mdot, X=x_0, Z=z, R=radius, Ma=mass)

    # extract results for comparison with obs
    result = np.recarray(
        (1,), dtype=[("tdel", np.float64), ("E_b", np.float64), ("alpha", np.float64)]
    )
    # assign elements
    result.tdel = res[1]
    result.E_b = res[2]
    result.alpha = res[0]
    result.mdot = mdot

    return result
