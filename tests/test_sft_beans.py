#!/usr/bin/env python3

# Short Functional Test (SFT) of BEANSp
#
# Beans.do_run() does not return anything,
# hence we are testing for random run-time errors

import os

from beansp.beans import Beans

local_data_dir = 'data'

data_dir = os.path.join(os.path.dirname(
    os.path.abspath(__file__)), local_data_dir)

if not os.path.isdir(data_dir):
    os.makedirs(data_dir)

# directory/filename_base
# no file is actually created, but the directory must exist,
# otherwise the test fails!
written_data_filename_base = os.path.join(data_dir, 'dummy')

B = Beans(nwalkers=20,
          nsteps=10,
          run_id=written_data_filename_base,
          gti_checking=0,
          restart=True)
B.do_run(analyse=False)
