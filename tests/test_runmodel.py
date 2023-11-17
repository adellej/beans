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

    # Ordering of parameters in tuple is (from v2.0.0 onwards):
    # X, Z, Q_b, dist, xi_b, xi_p, mass, radius, f_a, f_E
    # theta = 0.5, 0.015, 0.2, 2.1, 3.5, 0.108, 0.90, 0.5, 1.4, 11.2
    # theta = 0.5, 0.02, 0.36, 4.0, 0.56, 1.1, 1.4, 11.8, 1.0, 1.0
    # Another mod to make more sensible outputs for v2.2.0 onwards,
    # adjusted the flux to mdot conversion
    theta = 0.5, 0.02, 0.36, 2.25, 0.56, 1.1, 1.4, 11.8, 1.0, 1.0

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

    # updated call to include the train, numburstsobs parameter, and add the
    # full model which is now also output
    # test, valid, full_model = runmodel(theta,y,tref,bstart,pflux, pfluxe, tobs,numburstssim,len(bstart),ref_ind,gti_checking, 1, st, et)
    test, valid, full_model = runmodel(theta, B)
    print ('test={}'.format(test))
    
    # updated here with the new ("physical" parameters), v2.0.0 and later
    # and new mdot-to-flux conversion, v2.2.0 and later
    # and updates to pySettle 1.3.0, removing the 0.65 factor

    exp = [ 0.14350014,  3.19096834,  5.22761251,    # burst times (d)
           16.74947202, 16.80537819, 18.85328316, 22.82975017,  # fluences
           58.41027422, 64.38414044, 70.46505538] # alphas

    result = np.allclose(test, exp)#, rtol=1e-4)

    assert result

    return
