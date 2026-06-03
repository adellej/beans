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

    # Ordering of parameters in tuple is (from v2.53.0dev onwards):
    # X, Z, Q_b, dist, xi_b, xi_p, mass, radius, f_t
    # time systematic error term f_t replaced the (deprecated) systematic
    # errors on the fluences and alphas
    # theta = 0.5, 0.015, 0.2, 2.1, 3.5, 0.108, 0.90, 0.5, 1.4, 11.2
    # theta = 0.5, 0.02, 0.36, 4.0, 0.56, 1.1, 1.4, 11.8, 1.0, 1.0

    # Another mod to make more sensible outputs for v2.2.0 onwards,
    # adjusted the flux to mdot conversion
    # theta = 0.5, 0.02, 0.36, 2.25, 0.56, 1.1, 1.4, 11.8, 1.0

    # And another, after fixing a bug with the flux_to_mdot in 2.68.1
    theta = 0.58, 0.013, 0.4, 3.5, 1.0, 1.0, 1.5, 11.8, 1.0

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

    # linear interpolation, no correction

    B = Beans(corr=None, config_file=None, interp='linear', 
        fluen=True, alpha=True, numburstssim=3, ref_ind=1, bc=2.21,
        obsname=path_to_data_file_obs, burstname=path_to_data_file_bursts,# gtiname=path_to_data_file_gti)
        test_model=False)

    # spline interpolation, no correction

    Bs = Beans(corr=None, config_file=None, interp='spline', 
        fluen=True, alpha=True, numburstssim=3, ref_ind=1, bc=2.21,
        obsname=path_to_data_file_obs, burstname=path_to_data_file_bursts,# gtiname=path_to_data_file_gti)
        test_model=False)

    # spline interpolation, with corr_goodwin19 correction

    Bc = Beans(corr=beans.corr_goodwin19, config_file=None, interp='spline', 
        fluen=True, alpha=True, numburstssim=3, ref_ind=1, bc=2.21,
        obsname=path_to_data_file_obs, burstname=path_to_data_file_bursts,# gtiname=path_to_data_file_gti)
        test_model=False)

    # punkt_train, spline interpolation, with corr_goodwin19 correction
    # cannot currently run this without error, as the first (backwards)
    # simulation step cannot be performed

    Bp = Beans(corr=beans.corr_goodwin19, config_file=None, interp='spline', 
        fluen=True, alpha=True, bc=2.21, continuous=False, maxgap=2, 
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
    # added the reference burst time for v.2.32.0

    exp = [0.0837645978, 2.30514000, 3.34174103, 4.80458512,   # burst times (d)
           4.13526712, 4.14180596, 4.23639823, 4.51396060,     # fluences
           87.4369606, 99.9951515, 115.241809]                 # alphas

    result = np.allclose(test, exp, rtol=1e-3)
    assert result

    # some additional tests with different options
    # excluding (time) systematic error:

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

    exp_nomr = [ 0.46491736,  2.30514   ,  3.0527464 ,  4.04051396, 
                 3.89453988, 3.80815407,  3.76096016,  3.83814064,
                 76.18839969, 83.26951394, 95.6513451 ]
    result = np.allclose(test, exp_nomr, rtol=1e-3)
    assert result

    # and try the spline interpolation

    test, valid, full_model = runmodel(theta, Bs)
    print ('\nRunning with spline interpolation, theta={}\n  result={}'.format(theta,test))

    # which gives slight differences again in the predictions

    exp_spline = [ 0.11768532,   2.30514   ,   3.34364842,   4.82542518,
                   4.13300978,   4.13817863,   4.23738724,   4.52798589,
                   86.26987599, 100.07861092, 115.79040439]
    result = np.allclose(test, exp_spline, rtol=1e-3)
    assert result

    # add the correction

    test, valid, full_model = runmodel(theta, Bc)
    print ('\nRunning with goodwin19 correction, theta={}\n  result={}'.format(theta,test))

    # which gives not-so-slight differences in the predictions

    exp_corr = [ 0.86491876,   2.30514   ,   2.92768344,   4.73361625,
                 2.68662256,   2.69279975,   2.72743122,   2.96238376,
                 87.70798574,  96.33965845, 116.90341011]
    result = np.allclose(test, exp_corr, rtol=1e-3)
    assert result

    # additional tests of data for IGR 17498-2921

    path_to_data_file_obs = (
            pathlib.Path(__file__).resolve().parent.parent / "beansp/data" / "17498_obs_falanga12.txt")

    path_to_data_file_bursts = (
            pathlib.Path(__file__).resolve().parent.parent / "beansp/data" / "17498_bursts_8.txt")

    B17498 = Beans(config_file=None, corr=beans.corr_goodwin19, 
        interp='spline', smooth=0.1, 
        fluen=True, alpha=False, numburstssim=11, ref_ind=4, bc=1.0,
        obsname=path_to_data_file_obs, burstname=path_to_data_file_bursts,# gtiname=path_to_data_file_gti)
        test_model=False)

    theta_17498_ex = (1.77403857e-01, 1.67999877e-03, 2.03328677,
       6.38098467, 1.18448559, 1.21979245, 2.38079878, 1.19305420e+01)
    test, valid, full_model = runmodel(theta_17498_ex, B17498)
    print ('\nRunning for IGR J17498-2921 with theta={}\n  result={}'.format(
        theta_17498_ex,test))

    # These values copied from the last sample of the base22p run

    exp_times = [1.3539, 2.1199, 2.7784, 3.4019, 4.0217, 4.6417, 5.2664, 5.8978, 6.5352, 7.1762, 7.8321, 8.5168, 9.2304, 9.9863, 10.8276, 11.7921, 12.8447, 13.9178, 15.0281, 16.4121, 18.1007, 20.3002]
    exp_e_b = [2.8379, 2.7342, 2.6919, 2.6874, 2.6876, 2.6931, 2.7013, 2.7086, 2.7129, 2.7312, 2.7605, 2.7925, 2.8314, 2.9044, 2.9946, 3.0455, 3.0562, 3.0746, 3.2326, 3.4798, 3.8977]
    fpred = np.array(exp_e_b)*B17498.fluen_fac/theta_17498_ex[4]/theta_17498_ex[3]**2
    exp_alpha = [249.107, 243.861, 242.141, 241.961, 241.969, 242.192, 242.519, 242.806, 242.977, 243.732, 245.103, 246.568, 248.645, 253.028, 259.692, 264.639, 265.82, 267.986, 284.209, 298.512, 315.251]

    result1 = np.allclose(exp_times, full_model['time'],rtol=1e-4)
    result2 = np.allclose(fpred, full_model['fluen'],rtol=1e-4)
    result3 = np.allclose(exp_alpha, full_model['alpha'],rtol=1e-4)

    assert result1 & result2 & result3

    return
