import pathlib

import numpy as np
import tables
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import stats

# read in file
path_to_data_file = pathlib.Path(__file__).parent / "data" / "mr_prior_fit.txt"
data = np.loadtxt(path_to_data_file)
mu = data[:, 0]
sigma = data[:, 1]

# create grid points for mass
Marray = np.linspace(0.2, 3.0, 97)

# only go to 97 because last 3 mass grid points have very lowly populated posteriors for R and cannot get a good fit of mu and sigma.

# make prior function:
# required arrays: mu, sissgma and Marray
def mr_prior(M, R):
    # hard mass limits: M = 0.2-2.5, R = 9.5-16
    # exclude values outside of domain of interpolation:
    if M > 2.5 or M < 0.2:
        return -np.inf
    if R > 16 or R < 9.5:
        return -np.inf
    else:
        # extract closest mass grid coordinate:
        index = np.searchsorted(Marray, M)

        for i in [index]:
            y = stats.norm.pdf(R, mu[i], sigma[i])

            return y
