""" test to see if runmodel functions as expected """
import sys, os
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '../beans')
from beansp.burstrain import *
from beansp.run_model import runmodel
import numpy as np
from beansp import Beans
from beansp.get_data import get_obs
import pathlib


def test_run_model():

    numburstssim = 3 # this needs to be an integer value of half the number of bursts you want to simulate. I.e. simulate this many from the reference burst in either direction. Don't forget to account for missed bursts!
    numburstsobs = 4  # number of observed bursts in your dataset

    # Ordering of parameters in tuple is (from v2.0.0 onwards):
    # X, Z, Q_b, dist, xi_b, xi_p, mass, radius, f_a, f_E
    # theta = 0.5, 0.015, 0.2, 2.1, 3.5, 0.108, 0.90, 0.5, 1.4, 11.2
    theta = 0.5, 0.02, 0.36, 4.0, 0.56, 1.1, 1.4, 11.8, 1.0, 1.0

    path_to_data_file_obs = (
            pathlib.Path(__file__).resolve().parent.parent / "beansp/data" / "1808_obs.txt"
        )

    path_to_data_file_bursts = (
            pathlib.Path(__file__).resolve().parent.parent / "beansp/data" / "1808_bursts.txt"
        )

    path_to_data_file_gti = (
            pathlib.Path(__file__).resolve().parent.parent / "beansp/data" / "1808_gti.txt"
        )


    # bstart0, bstart, fluen, fluene, obs, obs_err, pflux, pfluxe, tobs, st, et = get_obs(ref_ind=1, bc=2.21, obsname=path_to_data_file_obs, burstname=path_to_data_file_bursts, gtiname=path_to_data_file_gti)
    B = Beans(ref_ind=1, bc=2.21, obsname=path_to_data_file_obs, burstname=path_to_data_file_bursts)#, gtiname=path_to_data_file_gti)

    ref_ind = 1
    # tref = bstart[ref_ind]
    tref = B.tref
    gti_checking =0
    # this is for emcee:
    y = B.y
    yerr = B.yerr
    # y = obs
    # yerr = obs_err
    x = 0 # in our case, we do not require x (independent variables), however for input into MCMC we need to define a x

    # updated call to include the train, numburstsobs parameter, and add the
    # full model which is now also output
    # test, valid, full_model = runmodel(theta,y,tref,bstart,pflux, pfluxe, tobs,numburstssim,len(bstart),ref_ind,gti_checking, 1, st, et)
    test, valid, full_model = runmodel(theta,B)
    
    # updated here with the new ("physical" parameters), v2.0.0 and later
    exp = [7.68277452e-02, 2.71329835, 3.86210558, # burst times (d)
           5.67417971, 5.86224830, 6.32795291, 7.25595091, # fluences
           95.2263036, 101.222238, 108.805858] # alphas

    result = np.allclose(test, exp)#, rtol=1e-4)

    assert result

    return
