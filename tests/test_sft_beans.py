#!/usr/bin/env python3

# Short Functional Test (SFT) of BEANSp
#
# Beans.do_run() does not return anything,
# hence we are testing for random run-time errors

import os

from beansp.beans import Beans


def test_SFT():
    print("BEANSp Short Functional Test (SFT)...")
    local_data_dir = 'data'
    data_dir = os.path.join(os.path.dirname(
        os.path.abspath(__file__)), local_data_dir)
    if not os.path.isdir(data_dir):
        os.makedirs(data_dir)

    # directory/filename_base
    # output data files written to this directory, must exist,
    # otherwise the test fails!
    written_data_filename_base = os.path.join(data_dir, 'dummy')

    # MCU note: 1) test fails if there is no data and restart is set to True
    #           2) test fails if there is data and restart is set to False
    # For automated github action test, 1) is the case, all si good.
    # It's up to the tester to manage other cases when this test is
    # run manually.
    B = Beans(nwalkers=20,
              nsteps=10,
              run_id=written_data_filename_base,
              gti_checking=0,
              restart=False)
    B.do_run(plot=False, analyse=False)


if __name__ == "__main__":
    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
