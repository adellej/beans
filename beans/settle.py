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
        #Function used to create fitting polyomial for Recurrence scaling
        def f(coff, a, b, c):
            form = [1, a, b, c, a ** 2, a * b, a * c, b ** 2, b * c, c ** 2, a ** 3, a ** 2 * b, a ** 2 * c, a * b ** 2,
                    a * b * c, a * c ** 2, b ** 3, b ** 2 * c, b * c ** 2, c ** 3]
            result = np.dot(coff, form)
            return result
        #3rd order polynomial model
        coff = [0.00000000e+00,  6.47786754e+00, - 6.62057703e-01, - 4.94660078e-01,
                - 1.34437019e+01,  2.88410912e+00,  1.28543957e+01, - 1.88454844e-02,
                - 7.08280358e-01,  3.13151175e+00,  4.91757270e+00, - 1.24429141e+00,
                2.68123025e+01,  2.41494411e-02, - 7.97613423e+00, - 6.93040409e+01,
                - 1.87299631e-03,  5.86797810e-01,  1.08885313e+01,  1.93331966e+00]
        intercept = -1.3301585711542072
        Rscalef = (f(coff,x_0,result.tdel,mdot)+intercept)
        if kwargs['scaling']:
            E_b = result.E_b - Fscalef
            if result.mdot >0.095:
                tdel = result.tdel - Rscalef
            else:
                if x_0 < 0.55:
                    tdel = result.tdel
                else:
                    tdel = result.tdel*0.8*(0.08/mdot)
            alpha = result.alpha*(tdel/result.tdel)/(E_b/result.E_b)
            result.E_b = E_b
            result.tdel = tdel
            result.alpha = alpha
        return result
