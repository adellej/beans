""" Tests the data reading in routine """
import numpy as np
from beansp import Beans
from beansp.get_data import get_obs
import pathlib

def test_get_obs():

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
    B = Beans(obsname=path_to_data_file_obs, burstname=path_to_data_file_bursts, gtiname=path_to_data_file_gti)
    get_obs(B)

    bstart_err = B.bstart_err

    obsexp = [0.41362, 3.18427, 4.42677, 2.6196, 2.6486, 2.9904, 3.4599, 107.0, 105.7, 121.1] 
    # obserrexp = [0.0069, 0.0069, 0.0069, 0.0207, 0.0178, 0.0173, 0.0221, 2.0, 1.9, 2.3]
    obserrexp = [bstart_err, bstart_err, bstart_err, 0.0207, 0.0178, 0.0173, 0.0221, 2.0, 1.9, 2.3]

    result = np.allclose(obsexp, B.y)
    assert result
    result = np.allclose(obserrexp, B.yerr)
    assert result

    return


if __name__ == "__main__":
    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
