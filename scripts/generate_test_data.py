# This script will produce sample outputs for both the "ensemble" and burst
# train mode, for testing
#
# Run this script with
# python generate_test_data.py

from beansp import Beans
from beansp import beans

print('Generating test data with beansp v{}...'.format(beans.__version__))

# ensemble mode run (on GS 1826-24 data)

Be = Beans(obsname=None, run_id='ensemble', nwalkers=100, nsteps=1000, burstname='beansp/data/1826_bursts.txt', theta= (0.5, 0.02, 1.0, 6, 1.0, 1.0, 1.7, 11.2), alpha=False, threads=4, test_model=True, restart=False)

Be.do_run() 

# burst train mode run (on SAX J1808.4-3658 data)

Bt = Beans(run_id='train', nwalkers=100, nsteps=1000, burstname='beansp/data/1808_bursts.txt', obsname='beansp/data/1808_obs.txt', prior=beans.prior_1808, alpha=False, threads=4, test_model=True, restart=False)

Bt.do_run()

print('... done')
