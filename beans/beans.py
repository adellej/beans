"""Main module."""

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

# -------------------------------------------------------------------------#
## load local  modules
from settle import settle
from burstrain import *
from runmodel import runmodel
from get_data import get_obs
from mrprior import mr_prior
from likelihood import *

# -------------------------------------------------------------------------#
# Begin the emcee initialisation:
# -------------------------------------------------------------------------#

# nwalkers and nsteps are the number of walkers and number of steps for emcee to do:
ndim, nwalkers = 10, 100
nsteps = 2
# Set starting value for each theta parameter:
# Recall odering, theta: X, Z, Q_b, f_a, f_E, r1, r2, r3
theta = 0.5, 0.015, 0.2, 2.1, 3.5, 0.108, 0.90, 0.5, 1.4, 11.2
# Set run id, should update text file to describe what this run is doing.
run_id = "test1.0"
# Set number of threads for emcee to use (e.g. number of cores your computer has)
threads = 8

numburstssim = 3 # this needs to be an integer value of half the number of bursts you want to simulate. I.e. simulate this many from the reference burst in either direction. Don't forget to account for missed bursts!
numburstsobs = 4  # number of observed bursts in your dataset

obs, obs_err, gti = get_data(
    obsname="obs.txt", burstname="bursts.txt", gtiname="gtis.dat"
)
bstart = obs[0:4]
# Define index of reference burst (should be middle of predicted burst train). This burst will not be simulated but will be used as a reference to predict the other bursts.
ref_ind = 1
tref = bstart[ref_ind]

# Option to turn on gti time checking (1 for on, 0 for off):
gti_checking = 0

# Get start and end gtis:
st = np.zeros(len(gti))
et = np.zeros(len(gti))
for i in range(0, len(gti)):
    st[i] = gti[i][0]
    et[i] = gti[i][1]

y = obs
yerr = obs_err
x = 0 # in our case, we do not require x (independent variables), however for input into MCMC we need to define a x

print(y, yerr)

# test model works:
test, valid = runmodel(
    theta, y, tref, bstart, pflux, pfluxe, tobs, numburstssim, ref_ind, gti_checking
)
print("result: ", test, valid)
# Ideally want to format below so the identity of which value goes with which burst is clear

# Testing the various functions. Each of these will display the likelihood value, followed by the model-results "blob"
print("lnprior:", lnprior(theta))
print("lnlike:", lnlike(theta, x, y, yerr))
print("lnprob:", lnprob(theta, x, y, yerr))
print(theta)

## Running the chain

# This section now defines the initial walker positions and next defines the chain and runs emcee.

pos = [theta + 1e-4 * np.random.randn(ndim) for i in range(nwalkers)]

print("Ready to run", run_id, "with", len(pos), "walkers")

# Define the sampler for emcee. Threads argument should be set to how many cores your computer has.
# sampler.reset()
sampler = emcee.EnsembleSampler(
    nwalkers, ndim, lnprob, args=(x, y, yerr), threads=threads
)

print("Running sampler...")
# Option to run without progress bar:
# run = sampler.run_mcmc(pos, 1000)

# Option to run with progress bar:
width = 60

for i, result in enumerate(sampler.sample(pos, iterations=nsteps)):
    n = int((width + 1) * float(i) / nsteps)
    sys.stdout.write("\r[{0}{1}]".format("#" * n, " " * (width - n)))
sys.stdout.write("\n")

print("...sampler run " + run_id + " complete")

# -------------------------------------------------------------------------#
# Save multi-threaded sampler:


def __getstate__(sampler):
    sampler_dict = sampler.__dict__.copy()
    del sampler_dict["pool"]
    return sampler_dict


sampler_nopool = __getstate__(sampler)
pickle.dump(sampler_nopool, open(run_id + "_chain.p", "wb"))


def __setstate__(self, state):
    self.__dict__.update(state)


# -------------------------------------------------------------------------#
# -------------------------------------------------------------------------#
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
