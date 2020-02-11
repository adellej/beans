

from beans import Beans

B = Beans(run_id="17511/test3", theta= (0.64, 0.02, 0.15, 2.1, 3.5, 0.8, 0.9, 0.1, 1.4, 11.2),  numburstssim=4, numburstsobs=7, bc=1.0, ref_ind=4, nsteps=1000, obsname='/Users/adelle/Documents/MCMC_burstcode_beans/beans/data/17498_obs.txt', burstname='/Users/adelle/Documents/MCMC_burstcode_beans/beans/data/17498_bursts.txt', gtiname=None)

B.do_run()

#B.do_analysis()
