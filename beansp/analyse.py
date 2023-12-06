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


def get_param_uncert_obs(param_array, percentile=[16,50,84]):
    '''
    Get uncertainties on predicted parameters. This routine accepts a
    list-of-lists with (for example) the start times of bursts for each
    walker and each timestep:
      [ [ t1, t2, t3, t4... tn ], 
        [ t1, t2, t3, t4... tn ],
            .
            .
            .
    dimensions (nsteps*walkers, nbursts) and calculates the 50th
    percentile value and +/- 1-sigma confidence limits for each of the t1,
    t2, t3, returning a list-of-tuples for each parameter, with the 50th
    percentile and the upper & lower uncertainties for each parameter in
    each tuple.

    Now modified to handle inhomogeneous arrays, where the number of
    elements is not constant; this can arise where we have multiple
    regions of parameter space that are being explored simultaneously

    :param param_array: array of dimension (nsteps*walkers, nbursts) listing the predicted values

    :return: dict with one or more entries, each with key of the number of bursts and a list of tuples giving parameter centroid and limits
    '''

    try:
        numburstssim = np.shape(param_array)[1]
    except:
        # if the shape fails, then it indicates we have an inhomogeneous
        # set of model results

        numburstssim = [len(x) for x in param_array]
        result = dict()
        for _n in set(numburstssim):
            _param_array = [x for x in param_array if len(x) == _n]
            result[_n] = get_param_uncert_obs(_param_array, percentile)[_n]

        return result

    plist = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
        zip(*np.percentile(param_array, percentile, axis=0)))

    return {numburstssim: list(plist)}


def get_param_uncert_part(param_array, percentile=[16,50,84], partition=None):
    '''
    As for get_param_uncert_obs, but uses a "partition" array which
    divides the set of lists into discrete groups for analysis

    :param param_array: array of dimension (nsteps*walkers, nbursts) listing the predicted values
    :param partition: arbitrary array of same length as the first dimension of param_array

    :return: dict with one or more entries, each with key from the partition array and a list of tuples giving parameter centroid and limits
    '''

    if partition is None:
        return get_param_uncert_obs(param_array, percentile)

    if len(partition) != len(param_array):
        print ('** ERROR ** partition array is wrong shape ({} != {},?), ignoring'.format(len(partition),len(param_array)))
        partition = None

    # now loop over the discrete partition elements and generate the stats

    result = dict()
    for _n in set(partition):
        _param_array = [x for i, x in enumerate(param_array) if partition[i] == _n]
        # this is fucking stupid, would be unnecessary if you could index
        # a dict values
        result[_n] = tuple(get_param_uncert_obs(_param_array, percentile).values())[0]

    return result


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

