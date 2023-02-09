""" Runs settle based on the input parameters and extracts the recurrence
time, fluence and alphas """

import settler as se
import numpy as np


def settle(base, z, x_0, mdot, cfac, mass, radius,**kwargs):
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
        #Scaling factor for Burst Fluence
        Fscalef = 0.03940475*result.E_b**2 -0.04471508*result.E_b - 0.22456746
        #Function used to create fiting polyomial for Recurrence scaling
        def f(coff, a, b, c):
            form = [1, a, b, c, a ** 2, a * b, a * c, b ** 2, b * c, c ** 2, a ** 3, a ** 2 * b, a ** 2 * c, a * b ** 2,
                    a * b * c, a * c ** 2, b ** 3, b ** 2 * c, b * c ** 2, c ** 3]
            result = np.dot(coff, form)
            return result
        #3rd order polynomial model
        coff = [0.00000000e+00,  2.96702341e+01, - 4.55179814e+01,  4.61933743e+01,
                - 6.08323104e+00, - 2.52189701e+02, - 1.55928328e+02,  7.09004523e+02,
                1.18553796e+02, - 7.26483562e+01,  4.60036188e+00, - 1.54955426e+02,
                2.79606110e+01,  3.48199814e+03, - 8.00558342e+01,  2.29525806e+02,
                8.87439121e+01, - 6.17308299e+03,  1.26544256e+03, - 4.86702849e+01]
        intercept = -7.1002920818991395
        Rscalef = (f(coff,x_0,z,mdot)+intercept)

        if kwargs['scaling']:
            E_b = result.E_b - Fscalef
            tdel = result.tdel - Rscalef
            result.E_b = E_b
            result.tdel = tdel
            alpha = result.alpha*(tdel/result.tdel)/(E_b/result.E_b)
            result.E_b = E_b
            result.alpha = alpha
        return result
