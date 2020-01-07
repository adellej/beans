


from beans import Beans


B = Beans(run_id='1808/test5', nsteps=1000, obsname='/Users/adelle/Documents/MCMC_burstcode_beans/beans/data/1808_obs.txt', burstname='/Users/adelle/Documents/MCMC_burstcode_beans/beans/data/1808_bursts.txt', gtiname='/Users/adelle/Documents/MCMC_burstcode_beans/beans/data/1808_gti.txt')

#B.do_run() 

B.do_analysis()

# import emcee
# import numpy as np 
# import h5py

# reader = emcee.backends.HDFBackend(filename='1808/test2')

# print(reader.get_chain())
