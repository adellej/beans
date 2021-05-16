import numpy as np

# -------------------------------------------------------------------------#
## load local  modules
from .burstrain import *
from .run_model import runmodel
from .get_data import get_obs
from .mrprior import mr_prior
from .get_data import *
from .run_emcee import runemcee


# def get_param_uncert_obs1(param_array, numburstssim):
#     # Get uncertainties on individual parameters:
#     p1, p2, p3, p4 ,p5, p6, p7 = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(param_array, [16, 50, 84], axis=0)))
#     # this will return
#     return p1, p2, p3, p4, p5, p6, p7

def get_param_uncert_obs(param_array, numburstssim):
    # Get uncertainties on individual parameters:
    plist = list()
    for i in range(0,numburstssim):
        plist = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(param_array, [16, 50, 84], axis=0)))
    plist2 = list()
    plist3 = list(plist)
    for i in range(0,numburstssim):
        plist2.append(plist3[i])
    return plist2

def get_param_uncert(param_array):
    # Get uncertainties on individual parameters:
    p = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(param_array, [16, 50, 84], axis=0)))
    return p
