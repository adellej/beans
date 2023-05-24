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
from scipy.stats import gaussian_kde
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
import ast
import matplotlib.gridspec as gridspec

# -------------------------------------------------------------------------#
## load local  modules
from .settle import settle
from .burstrain import *
from .run_model import runmodel
from .get_data import get_obs
from .mrprior import mr_prior
from .get_data import *
from .run_emcee import runemcee


def get_param_uncert_obs(param_array, numburstssim, percentile=[16,50,84]):
    '''
    Get uncertainties on predicted parameters. This routine accepts a
    list-of-lists with (for example) the start times of bursts for each
    walker and each timestep:
      [ [ t1, t2, t3, t4... tn ], 
        [ t1, t2, t3, t4... tn ],
            .
            .
            .
    and calculates the 50th percentile value and +/- 1-sigma confidence
    limits for each of the t1, t2, t3, returning a list-of-tuples for each
    parameter, with the 50th percentile and the upper & lower uncertainties 
    for each parameter in each tuple.

    :param param_array: array listing the predicted values
    :param numburstssim: number of events per instance (2nd dimension of
      param_array)

    :return: list of tuples giving parameter centroid and limits
    '''

    plist = list()
    for i in range(0,numburstssim):
        plist = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
            zip(*np.percentile(param_array, percentile, axis=0)))
    plist2 = list()
    plist3 = list(plist)
    for i in range(0,numburstssim):
        plist2.append(plist3[i])

    return plist2


def get_param_uncert(param_array, percentile=[16,50,84]):
    '''
    Calculate uncertainties on individual (scalar) parameters.

    :param param_array: array of parameter values to calculate percentiles
    :param percentile: percentiles to calculate, normally 1-sigma
    :return: array of percentile values; by default, the 50th percentile,
      upper error and lower error
    '''

    v = np.percentile(param_array, percentile, axis=0)

    return np.array([v[1], v[2]-v[1], v[1]-v[0]])

