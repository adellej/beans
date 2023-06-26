""" Runs settle based on the input parameters and extracts the recurrence
time, fluence and alphas """

from pySettle import settler as se
import numpy as np


def settle(base, z, x_0, mdot, cfac, mass, radius, **kwargs):
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
    # Scaling factor for Burst Fluence
    Fscalef = 0.03940475 * result.E_b ** 2 - 0.04471508 * result.E_b - 0.22456746

    # Function used to create fitting polyomial for Recurrence scaling
    def f(coff, a, b, c):
        # 3rd degree
        # form = [1, a, b, c, a**2, a*b, a*c, b**2, b*c, c**2, a**3, a**2*b, a**2*c, a*b**2, a*b*c, a*c**2, b**3, b**2*c, b*c**2, c**3]
        # 2nd degree
        form = [1, a, b, c, a ** 2, a * b, a * c, b ** 2, b * c, c ** 2]
        # 1st degree
        # form = [1,a,b,c]
        result = np.dot(coff, form)
        return result

    # 3rd order polynomial model
    coff = [0.00000000e+00,  5.21062320e+00,  2.00153869e-01,  2.41105856e+01,
     -7.20744408e+00,  1.36160251e+00,  2.11588927e+00, -1.99451593e-02,
     -5.79474586e-01, -4.19491308e+01]


    intercept = -4.836294980805745
    Rscalef = (f(coff, x_0, result.tdel, mdot) + intercept)
    if kwargs['scaling']:
        E_b = result.E_b - Fscalef
        if result.mdot > 0.095:
            tdel = result.tdel - Rscalef
        else:
            if x_0 < 0.55:
                tdel = result.tdel
            else:
                tdel = result.tdel * 0.8 * (0.08 / mdot)
        alpha = result.alpha * (tdel / result.tdel) / (E_b / result.E_b)
        result.E_b = E_b
        result.tdel = tdel
        result.alpha = alpha
    return result
