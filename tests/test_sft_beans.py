#!/usr/bin/env python3

# Short Functional Test (SFT) of BEANSp
#
# Beans.do_run() does not return anything,
# hence we are testing for random run-time errors

from beansp.beans import Beans

B = Beans(nwalkers=20,
          nsteps=10,
          run_id="../beansp/data/1808/dummy",
          gti_checking=0)
B.do_run(analyse=False)
