
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

def init():
    # -------------------------------------------------------------------------#
    # THESE ARE THE ONLY PARAMETERS YOU NEED TO ENTER
    # -------------------------------------------------------------------------#
    # nwalkers and nsteps are the number of walkers and number of steps for emcee to do:
    ndim, nwalkers = 10, 100
    nsteps = 100
    # Set starting value for each theta parameter:
    # Recall odering, theta: X, Z, Q_b, f_a, f_E, r1, r2, r3
    theta = 0.5, 0.015, 0.2, 2.1, 3.5, 0.108, 0.90, 0.5, 1.4, 11.2
    # Set run id, should update text file to describe what this run is doing.
    run_id = "1808/test1"
    # Set number of threads for emcee to use (e.g. number of cores your computer has)
    threads = 4
    os.environ["OMP_NUM_THREADS"] = f"{threads}"

    numburstssim = 3 # this needs to be an integer value of half the number of bursts you want to simulate. I.e. simulate this many from the reference burst in either direction. Don't forget to account for missed bursts!
    numburstsobs = 4  # number of observed bursts in your dataset

    # Define index of reference burst (should be middle of predicted burst train). This burst will not be simulated but will be used as a reference to predict the other bursts.
    ref_ind = 1
    # Option to turn on gti time checking (1 for on, 0 for off):
    gti_checking = 0

    # set location of your observation data files:
    obsname='/Users/adelle/Documents/MCMC_burstcode_beans/beans/data/1808_obs.txt'
    burstname='/Users/adelle/Documents/MCMC_burstcode_beans/beans/data/1808_bursts.txt'
    gtiname='/Users/adelle/Documents/MCMC_burstcode_beans/beans/data/1808_gti.txt'

    # and finally the bolometric correction to apply to your persistent flux (1.0 if they are already bolometric fluxes):
    bc = 2.21

    # if your run crashed and you would like to restart from a previous run, with run_id above, set this to True:
    restart = False

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

    return ndim, nwalkers, nsteps, run_id, theta, x, y, yerr, tref, bstart, pflux, pfluxe, tobs, numburstssim, numburstsobs, bc, ref_ind, gti_checking, fluen, restart