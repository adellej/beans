""" test to see if runmodel functions as expected """
import sys, os
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '../beans')
from beansp.burstrain import *
from beansp.run_model import runmodel
import numpy as np
from beansp import Beans, beans
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
    B = Beans(corr=None, config_file=None, interp='linear', 
        fluen=True, alpha=True, numburstssim=3, ref_ind=1, bc=2.21,
        obsname=path_to_data_file_obs, burstname=path_to_data_file_bursts,# gtiname=path_to_data_file_gti)
        test_model=False)

    Bs = Beans(corr=None, config_file=None, interp='spline', 
        fluen=True, alpha=True, numburstssim=3, ref_ind=1, bc=2.21,
        obsname=path_to_data_file_obs, burstname=path_to_data_file_bursts,# gtiname=path_to_data_file_gti)
        test_model=False)

    Bc = Beans(corr=beans.corr_goodwin19, config_file=None, interp='spline', 
        fluen=True, alpha=True, numburstssim=3, ref_ind=1, bc=2.21,
        obsname=path_to_data_file_obs, burstname=path_to_data_file_bursts,# gtiname=path_to_data_file_gti)
        test_model=False)

    # updated call to include the train, numburstsobs parameter, and add the
    # full model which is now also output
    # test, valid, full_model = runmodel(theta,y,tref,bstart,pflux, pfluxe, tobs,numburstssim,len(bstart),ref_ind,gti_checking, 1, st, et)
    test, valid, full_model = runmodel(theta, B)
    print ('\nRunning with theta={}\n  result={}'.format(theta,test))
    
    # updated here with the new ("physical" parameters), v2.0.0 and later
    # and new mdot-to-flux conversion, v2.2.0 and later
    # and updates to pySettle 1.3.0, removing the 0.65 factor
    # and updates to the way the time offset is calculated, 2.11.x and later
    # and the fix to the mdot_Edd function, v2.21.1

    exp = [  0.3371498 ,  3.84590496,  6.76410934,    # burst times (d)
           17.94856258, 17.35869927, 20.12729269, 27.51898359,  # fluences
           60.40744506, 66.85451224, 74.21961878] # alphas

    result = np.allclose(test, exp, rtol=1e-3)
    assert result

    # some additional tests with different options

    theta_nosym = theta[:8]
    test, valid, full_model = runmodel(theta_nosym, B)
    print ('\nRunning with theta={}\n  result={}'.format(theta_nosym,test))

    # without the systematic errors, should be no difference to the model

    result = np.allclose(test, exp, rtol=1e-3)
    assert result

    # now with canonical mass and radius

    theta_nomr = theta[:6]
    test, valid, full_model = runmodel(theta_nomr, B)
    print ('\nRunning with theta={}\n  result={}'.format(theta_nomr,test))

    # with the canonical mass and radius there are slight differences

    exp_nomr = [ 0.55146657,  3.64700296,  5.89032248,
                 14.92333594, 14.96612381, 17.06069978, 21.61310762,
                 62.70263335, 69.5733182 , 76.85360513 ]
    result = np.allclose(test, exp_nomr, rtol=1e-3)
    assert result

    # and try the spline interpolation

    test, valid, full_model = runmodel(theta, Bs)
    print ('\nRunning with spline interpolation, theta={}\n  result={}'.format(theta,test))

    # which gives slight differences again in the predictions

    exp_spline = [ 0.37045536,  3.85703407,  6.82600001, 
                   16.91872503, 17.25034895, 20.18742875, 27.78212457,
                   60.04565281, 66.9551983 , 74.37555756 ]
    result = np.allclose(test, exp_spline, rtol=1e-3)
    assert result

    # add the correction

    test, valid, full_model = runmodel(theta, Bc)
    print ('\nRunning with goodwin19 correction, theta={}\n  result={}'.format(theta,test))

    # which gives not-so-slight differences in the predictions

    exp_corr = [ 0.41762404,  3.19702031,  4.48089481,
                 10.99004728, 11.32612194, 12.51036817, 14.61326653,
                 60.62210726, 65.20925343, 70.08695522 ]
    result = np.allclose(test, exp_corr, rtol=1e-3)
    assert result

    return
