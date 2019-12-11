

from beans import Beans

B = Beans(ndim=10, nwalkers=200, nsteps=1000, run_id="1808/tomtest1", 
obsname='/home/tom//Documents/Thermonuclear bursts project/beans install/beans/data/1808_obs.txt', 
burstname='/home/tom//Documents/Thermonuclear bursts project/beans install/beans/data/1808_bursts.txt', 
gtiname='/home/tom//Documents/Thermonuclear bursts project/beans install/beans/data/1808_gti.txt', 
theta= (0.5, 0.015, 0.2, 2.1, 3.5, 0.108, 0.90, 0.5, 1.4, 11.2), numburstssim=3, numburstsobs=4, 
bc=2.21, ref_ind=1, gti_checking=0, threads = 4, restart=False)

B.do_run()

#/home/tom/Documents/Thermonuclear bursts project/beans install/beans/data/