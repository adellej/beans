

from beans import Beans

B = Beans(run_id="17511/test1", theta= (0.7, 0.02, 0.2, 2.1, 3.5, 0.108, 0.90, 0.5, 1.4, 11.2),  numburstssim=15, numburstsobs=17, bc=1.0, ref_ind=9, nsteps=1000, obsname='/Users/adelle/Documents/MCMC_burstcode_beans/beans/data/17511_obs.txt', burstname='/Users/adelle/Documents/MCMC_burstcode_beans/beans/data/17511_bursts.txt', gtiname=None)

B.do_run()
