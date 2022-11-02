

from beans import Beans

B = Beans(run_id="17498/test3", theta= (0.64, 0.02, 0.15, 2.1, 3.5, 0.7, 0.9, 0.1, 1.4, 11.2), \
 numburstssim=4, numburstsobs=7, bc=1.0, ref_ind=5, nsteps=10, \
 obsname='/home/adelleg/BEANS_new/beans/data/17498_obs.txt',\
  burstname='/home/adelleg/BEANS_new/beans/data/17498_bursts.txt', gtiname=None)

B.do_run()

#B.do_analysis()
