""" Tests the data reading in routine """
import numpy as np
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


    bstart0, bstart, fluen, obs, obs_err, pflux, pfluxe, tobs, st, et = get_obs(ref_ind=1, bc=2.21, gti_checking=1, obsname=path_to_data_file_obs, burstname=path_to_data_file_bursts, gtiname=path_to_data_file_gti)

    obsexp = [0.0, 2.770649999998568, 4.013149999998859, 2.6196, 2.6486, 2.9904, 3.4599, 107.0, 105.7, 121.1] 
    obserrexp = [0.0069, 0.0069, 0.0069, 0.0207, 0.0178, 0.0173, 0.0221, 2.0, 1.9, 2.3]

    result = np.allclose(obsexp, obs)
    assert result
    result = np.allclose(obserrexp, obs_err)
    assert result

    return


if __name__ == "__main__":
    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
