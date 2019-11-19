"""Main module. This will do the sampling, save the chains, and analyse the results."""

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
from likelihood import lnlike, lnprior, lnZprior, lnprob
from initialise import init

ndim, nwalkers, nsteps, run_id, theta, x, y, yerr, tref, bstart, pflux, pfluxe, tobs, numburstssim, numburstsobs, bc, ref_ind, gti_checking, fluen, restart = init()

# test model works:
test, valid = runmodel(theta, y, tref, bstart, pflux, pfluxe, tobs, numburstssim, ref_ind, gti_checking)
print("# -------------------------------------------------------------------------#")
print("Testing the model works..")
print("result: ", test, valid)
# Ideally want to format below so the identity of which value goes with which burst is clear

print("# -------------------------------------------------------------------------#")
# Testing the various functions. Each of these will display the likelihood value, followed by the model-results "blob"
print("Testing the prior and likelihood functions..")
print("lnprior:", lnprior(theta))
print("lnlike:", lnlike(theta, x, y, yerr))
print("lnprob:", lnprob(theta, x, y, yerr))
print("# -------------------------------------------------------------------------#")
print(f"The theta parameters will begin at: {theta}")
print("# -------------------------------------------------------------------------#")
## Running the chain
# we use multiprocessing to speed things up. Emcee parameters are defined in runemcee module. 

sampler = runemcee(nwalkers, nsteps, ndim, theta, lnprob, x, y, yerr, run_id, restart) # this will run the chains and save the output as a h5 file
print(f"Complete!")

# -------------------------------------------------------------------------#
# Finally analyse and display the results:

# # Displaying results
# import numpy as np
# from chainconsumer import ChainConsumer

# ndim = 10

# sampler1 = sampler.chain[:, :, 0:10]
# samples = sampler1[:, :, :].reshape(-1, ndim)

# c = ChainConsumer()

# cc = c.add_chain(
#     samples, parameters=["X", "Z", "Qb", "fa", "fE", "r1", "r2", "r3", "M", "R"]
# ).configure(smooth=False)


# fig = cc.plotter.plot()
# fig.savefig("plots/walks{}.png".format(run_id))

# -------------------------------------------------------------------------#

# finally reset emcee sampler for next runs.

sampler.reset()
