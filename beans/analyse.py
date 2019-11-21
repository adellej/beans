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
import ast
import matplotlib.gridspec as gridspec

# -------------------------------------------------------------------------#
## load local  modules
from settle import settle
from burstrain import *
from run_model import runmodel
from get_data import get_obs
from mrprior import mr_prior
from get_data import *
from run_emcee import runemcee


def get_param_uncert_obs1(param_array, numburstssim):
    # Get uncertainties on individual parameters:
    p1, p2, p3, p4 ,p5, p6, p7 = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(param_array, [16, 50, 84], axis=0)))
    # this will return 
    return p1, p2, p3, p4, p5, p6, p7
def get_param_uncert_obs(param_array, numburstssim):
    # Get uncertainties on individual parameters:
    p1, p2, p3, p4 ,p5, p6 = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(param_array, [16, 50, 84], axis=0)))
    # this will return 
    return p1, p2, p3, p4, p5, p6

def get_param_uncert(param_array):
    # Get uncertainties on individual parameters:
    p = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(param_array, [16, 50, 84], axis=0)))
    # this will return 
    return p
                                        