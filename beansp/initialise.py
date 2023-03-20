
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
# load local modules
from .settle import settle
from .burstrain import *
from .run_model import runmodel
from .get_data import get_obs
from .mrprior import mr_prior
from .get_data import *
from .run_emcee import runemcee

# -------------------------------------------------------------------------#
# Begin the emcee initialisation:
# -------------------------------------------------------------------------#

def init(ref_ind, gti_checking, obsname, burstname, gtiname, bc):
    """
    Function to read in and initialise the input parameters. Previously 
    we also passed ndim, nwalkers, theta, run_id, threads, numburstssim,
'   numburstsobs & restart, but these are not used

    :param ref_ind:
    :param gti_checking:
    :param obsname:
    :param burstname:
    :param gtiname:
    :param bc:
    :return: 
    """

    # -------------------------------------------------------------------------#
    # END OF PARAMETERS YOU NEED TO SET
    # -------------------------------------------------------------------------#

    # -------------------------------------------------------------------------#
    # Now let the code do it's thing
    # -------------------------------------------------------------------------#

    # get the data:
    if gti_checking:
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
