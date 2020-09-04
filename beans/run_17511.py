

from beans import Beans

B = Beans(run_id="17511/test3", theta= (0.6, 0.02, 0.2, 2.1, 3.5, 0.67, 0.5, 0.04, 1.4, 11.2),  numburstssim=12, numburstsobs=17, bc=1.0, ref_ind=8, nsteps=1000, obsname='/Users/adelle/Documents/MCMC_burstcode_beans/beans/data/17511_obs.txt', burstname='/Users/adelle/Documents/MCMC_burstcode_beans/beans/data/17511_bursts.txt', gtiname=None)

#B.do_run()

B.do_analysis()
