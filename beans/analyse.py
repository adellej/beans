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
from likelihood import lnlike, lnZprior, lnprob
from initialise import init


def analyse(name):
# run_id = "chains_1808/test1"
    burnin = 10
    thin = 5
    ndim, nwalkers, nsteps, run_id, theta, x, y, yerr, tref, bstart, pflux, pfluxe, tobs, numburstssim, numburstsobs, bc, ref_ind, gti_checking, fluen, restart = init()

    # load in sampler:
    reader = emcee.backends.HDFBackend(filename=name)
    print(reader.iteration)
    #sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr), backend=reader)

    tau = reader.get_autocorr_time()
    burnin = int(2 * np.max(tau))
    thin = int(0.5 * np.min(tau))
    samples = reader.get_chain(discard=burnin, flat=True, thin=thin)
    log_prob_samples = reader.get_log_prob(discard=burnin, flat=True, thin=thin)
    log_prior_samples = reader.get_blobs(discard=burnin, flat=True, thin=thin)

    print("burn-in: {0}".format(burnin))
    print("thin: {0}".format(thin))
    print("flat chain shape: {0}".format(samples.shape))
    print("flat log prob shape: {0}".format(log_prob_samples.shape))
    print("flat log prior shape: {0}".format(log_prior_samples.shape))

    all_samples = np.concatenate(
        (samples, log_prob_samples[:, None], log_prior_samples[:, None]), axis=1
    )

    labels = list(map(r"$\theta_{{{0}}}$".format, range(1, ndim + 1)))
    labels += ["log prob", "log prior"]

    figure = corner.corner(all_samples, labels=labels)
    figure.savefig('corner.pdf')


    # model = generate_burst_train( base,
    #     z,
    #     x_0,
    #     r1,
    #     r2,
    #     r3,
    #     mass,
    #     radius,
    #     bstart,
    #     pflux,
    #     pfluxe,
    #     tobs,
    #     numburstssim,
    #     run=run
    #     double=double,
    #     debug=debug,
    # )

analyse(name="chains_1808/test1.h5")