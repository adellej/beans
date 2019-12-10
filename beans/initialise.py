
"""Initialisation. This is the only place where you need to enter parameters."""

## Python packages required:
import matplotlib.pyplot as plt
import numpy as np
import emcee
import corner
import random
import math
import subprocess
from astropy.io import ascii
import pickle
from matplotlib.ticker import MaxNLocator
import sys
import idlsave
from scipy.stats.kde import gaussian_kde
import scipy.stats as stats
import matplotlib.mlab as mlab
import tables
from scipy.interpolate import interp1d
from chainconsumer import ChainConsumer
from multiprocessing import Pool
import os
import time

# -------------------------------------------------------------------------#
## load local  modules
from settle import settle
from burstrain import *
from run_model import runmodel
from get_data import get_obs
from mrprior import mr_prior
from get_data import *
from run_emcee import runemcee

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
    if gti_checking ==1:
        bstart, fluen, obs, obs_err, pflux, pfluxe, tobs, st, et = get_obs(ref_ind,bc,gti_checking,obsname, burstname, gtiname)
    else:
        bstart, fluen, obs, obs_err, pflux, pfluxe, tobs = get_obs(ref_ind,bc,gti_checking,obsname, burstname, gtiname)

    tref = bstart[ref_ind]

    # this is for emcee:
    y = obs
    yerr = obs_err
    x = 0 # in our case, we do not require x (independent variables), however for input into MCMC we need to define a x

    print("So the data that emcee will use is obs, obs_err:")
    print(y, yerr)

    return x, y, yerr, tref, bstart, pflux, pfluxe, tobs, fluen