""" Tests the data reading in routine """
import numpy as np
from beansp import Beans
from beansp.get_data import get_obs
import pathlib

def test_get_obs():

    # A couple of tests here to make sure the get_obs function works

    data_path = pathlib.Path(__file__).resolve().parent.parent / "beansp/data"
    print (data_path)

    # basic run with 1808 data, no peak fluxes

    path_to_data_file_obs = ( data_path / "1808_obs.txt")
    path_to_data_file_bursts = ( data_path / "1808_bursts.txt")
    path_to_data_file_gti = ( data_path / "1808_gti.txt")

    # bstart0, bstart, fluen, fluene, obs, obs_err, pflux, pfluxe, tobs, st, et = get_obs(ref_ind=1, bc=2.21, obsname=path_to_data_file_obs, burstname=path_to_data_file_bursts, gtiname=path_to_data_file_gti)
    B = Beans(obsname=path_to_data_file_obs,
        burstname=path_to_data_file_bursts, gtiname=path_to_data_file_gti,
        alpha=True, test_model=False)
    # this next step is redundant since the __init__ method is just
    # calling the same routine
    # get_obs(B)

    bstart_err = B.bstart_err

    # modified these arrays from v2.32.0 onwards, now includes the
    # reference burst time
    # obsexp = [0.41362, 3.18427, 4.42677, 2.6196, 2.6486, 2.9904, 3.4599, 107.0, 105.7, 121.1] 
    ## obserrexp = [0.0069, 0.0069, 0.0069, 0.0207, 0.0178, 0.0173, 0.0221, 2.0, 1.9, 2.3]
    # obserrexp = [bstart_err, bstart_err, bstart_err, 0.0207, 0.0178, 0.0173, 0.0221, 2.0, 1.9, 2.3]
    obsexp = [0.41362, 2.30514, 3.18427, 4.42677, 2.6196, 2.6486, 2.9904, 3.4599, 107.0, 105.7, 121.1] 
    obserrexp = [bstart_err, bstart_err, bstart_err, bstart_err, 0.0207, 0.0178, 0.0173, 0.0221, 2.0, 1.9, 2.3]

    result = np.allclose(obsexp, B.y)
    assert result
    result = np.allclose(obserrexp, B.yerr)
    assert result

    path_to_data_file_bursts = ( data_path / "1808_bursts_pflux.txt")

    # bursts with peak fluxes

    B2 = Beans(obsname=path_to_data_file_obs,
         burstname=path_to_data_file_bursts, gtiname=path_to_data_file_gti,
         alpha=False, test_model=False)

    result = np.allclose(obsexp[:8], B2.y)
    assert result

    bpflux = [215.11, 229.38, 232.42, 231.96]
    result = np.allclose(bpflux, B2.bpflux)

    # ensemble mode

    path_to_ensemble_file = ( data_path / "1826_bursts.txt")

    B3 = Beans(obsname=None, burstname=path_to_ensemble_file, test_model=False)

    obsexp_e = [5.14 , 4.177, 3.53 , 1.102, 1.126, 1.18 ]
    obserrexp_e = [0.07 , 0.01 , 0.004, 0.011, 0.016, 0.04 ]

    result = np.allclose(obsexp_e, B3.y)
    assert result

    result = np.allclose(obserrexp_e, B3.yerr)
    assert result

    return


if __name__ == "__main__":
    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
