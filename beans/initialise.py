
"""Initialisation. This is the only place where you need to enter parameters."""

## Python packages required:

# -------------------------------------------------------------------------#
## load local  modules
from .burstrain import *
from .run_model import runmodel
from .get_data import get_obs
from .mrprior import mr_prior
from .get_data import *
from .run_emcee import runemcee

# -------------------------------------------------------------------------#
# Begin the emcee initialisation:
# -------------------------------------------------------------------------#

def init(ndim, nwalkers, theta, run_id, threads, numburstssim, numburstsobs, ref_ind, gti_checking, obsname, burstname, gtiname,bc,restart):

    # -------------------------------------------------------------------------#
    # END OF PARAMETERS YOU NEED TO SET
    # -------------------------------------------------------------------------#

    # -------------------------------------------------------------------------#
    # Now let the code do it's thing
    # -------------------------------------------------------------------------#

    # get the data:
    if gti_checking == 1:
        tref, bstart, fluen, obs, obs_err, pflux, pfluxe, tobs, st, et = get_obs(ref_ind,bc,gti_checking,obsname, burstname, gtiname)
    else:
        tref, bstart, fluen, obs, obs_err, pflux, pfluxe, tobs = get_obs(ref_ind,bc,gti_checking,obsname, burstname, gtiname)
        st, et = None, None

    if burstname is not None:
        # tref = bstart[ref_ind]

        # this is for emcee:
        y = obs
        yerr = obs_err
        x = 0 # in our case, we do not require x (independent variables), however for input into MCMC we need to define a x
    else:
        # x, y, yerr, tref = None, None, None, 0.0
        x, y, yerr = None, None, None

    return x, y, yerr, tref, bstart, pflux, pfluxe, tobs, fluen, st, et
