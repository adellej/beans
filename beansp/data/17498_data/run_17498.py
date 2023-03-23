

from beans import Beans

B = Beans(ndim=10, nwalkers=200, nsteps=200, run_id="17498_4bursts", 
obsname='/home/tom//Documents/Thermonuclear bursts project/beans install/beans/data/17498_obs_new.txt', 
burstname='/home/tom//Documents/Thermonuclear bursts project/beans install/beans/data/17498_bursts_new.txt', 
gtiname='/home/tom//Documents/Thermonuclear bursts project/beans install/beans/data/17498_gti.txt', 
theta=(0.64, 0.01, 0.15, 2.1, 3.5, 0.60, 0.90, 0.0714, 1.4, 11.2), numburstssim=3, numburstsobs=4, 
bc=1, ref_ind=1, gti_checking=0, threads=8, restart=False)

B.do_run()
#B.do_analysis()


#theta=(0.64, 0.01, 0.15, 2.1, 3.5, 0.78, 0.90, 0.083, 1.4, 11.2), numburstssim=7, numburstsobs=7
#/home/tom/Documents/Thermonuclear bursts project/beans install/beans/data/