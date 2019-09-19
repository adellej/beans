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
# -------------------------------------------------------------------------#
## load local  modules
from settle import settle
from burstrain import *
from runmodel import runmodel

# -------------------------------------------------------------------------#
## Set location of settle:

#For Adelle
# adelle2:
#binpath = '/home/adelleg/Documents/project/Burst_code/settle2.0/settle'
#adelle1:
#binpath = '/home/adelle/burst_code/burst_code/settle2.0/settle'
# Mac:
#binpath = '/Users/adelle/Documents/burst_code/burst-code-new/settle2.0/settle'
# Oxygen:
#binpath = '/home/adelle/burst_code/settle2.0/settle'


# -------------------------------------------------------------------------#
# Now we can begin the emcee initialisation:
# -------------------------------------------------------------------------#
# Need tobs to do this:

#nwalkers and nsteps are the number of walkers and number of steps for emcee to do:
ndim, nwalkers = 10, 100
nsteps = 2
# Set starting value for each theta parameter:
# Recall odering, theta: X, Z, Q_b, f_a, f_E, r1, r2, r3
theta = 0.5, 0.015, 0.2, 2.1, 3.5, 0.108, 0.90, 0.5, 1.4, 11.2
# Set run id, should update text file to describe what this run is doing.
run_id='test1.0'
# Set number of threads for emcee to use (e.g. number of cores your computer has)
threads = 8

numburstssim = 3 #this needs to be an integer value of half the number of bursts you want to simulate. I.e. simulate this many from the reference burst in either direction

## To begin with the following data is required, in ascii format, with columns in the given order:
# burst observations (time (MJD), fluence, fluence error)
# observations (times (MJD), persistent flux measurements, pflux error)
# satellite gti data (start time, stop time (in MJD))
# predicted number of bursts, including missed bursts (this sets the number of bursts for the generate burst train routine to simulate)
# index of bursts that fall after missed bursts, and the predicted number of bursts missed

# We define y as the observational array and yerr the error of this array
s = idlsave.read('1808-match_split.dat')
obs_err = s.obs_err
obs = s.obs
tobs = s.tobs
pfluxe = s.pfluxe
tobs_err = s.tobs_err
gti = s.gti
bstart_err = s.bstart_err
pflux = s.pflux
bstart = s.bstart
fluen = s.fluen
alpha = s.alpha


# Define index of reference burst (should be middle of predicted burst train). This burst will not be simulated but will be used as a reference to predict the other bursts.
ref_ind = 1
tref = bstart[ref_ind]
n_burst = len(bstart)

# Option to turn on gti time checking (1 for on, 0 for off):
gti_checking = 0

tref = bstart[ref_ind]
n_burst = len(bstart)

# Get gtis:
st = np.zeros(len(gti))
et = np.zeros(len(gti))
for i in range (0, len(gti)):
    st[i] = gti[i][0]
    et[i] = gti[i][1]

y = obs
yerr = obs_err
x = 0 # in our case, we do not require x (independent variables), however for input into MCMC we need to define a x

print(y, yerr)

# test model works:
test, valid = runmodel(theta,y,tref,bstart,pflux, pfluxe, tobs,numburstssim,ref_ind,gti_checking)
print('result: ', test, valid)
# Ideally want to format below so the identity of which value goes with which burst is clear


# import data required for mass radius priors:

import numpy as np
import tables
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# read in file
file2 = tables.open_file('qlmxb_threep_base_all_out') # base model

table2 = file2.root.markov_chain0.data
markovobject2 = file2.get_node("/markov_chain0", "data")

# read data:
# create a dictionary so data format is useful
Data2 = {}
for name in markovobject2:
    Data2['{}'.format(name)] = name.read()
# get radius data
Radius2 = {}
for key in Data2:
    if "/markov_chain0/data/R_" in key:
        Radius2['{}'.format(key)] = (Data2['{}'.format(key)])

# now just get radii we want:
R_mcmc2 = []
for i in range(0,100):
    R_mcmc2.append(Radius2["/markov_chain0/data/R_{} (EArray(291827,)) ''".format(i)])
# exclude (mask) bad data where Rx has been set to zero as M_max is less than mass
for i in range(0,len(R_mcmc2)):
    R_mcmc2[i] = np.ma.array(R_mcmc2[i], mask=False)
    for j in range(0,len(R_mcmc2[0])):
        if R_mcmc2[i][j] == 0.0:
            R_mcmc2[i].mask[j] = True

# create grid points for mass
Marray = np.linspace(0.2,3.0,97) #only go to 97 because last 3 mass grid points have very lowly populated posteriors for R and cannot get a good fit of mu and sigma.

mu = []
sigma = []
for i in range(0,97):
    mu.append(np.mean(R_mcmc2[i]))
    sigma.append(np.std(R_mcmc2[i]))

# -------------------------------------------------------------------------#
## Now we define the functions that emcee requires
# define likelihood as a function of theta, x, y and yerr as this is what emcee expects as the inputs
from burstrain import *
def lnlike(theta, x, y ,yerr):

# define y = "data" parameters

    for x,i in zip([ x for x in range(0, len(bstart)-1) if x != ref_ind], [i for i in range(0, len(bstart)-1) if i != ref_ind]):
        globals()['t%s' % i] = y[x]
    for x,i in zip(range(len(bstart)-1, len(fluen)+len(bstart)-1),range(0,len(bstart))):
        globals()['Eb%s' % i] = y[x]
    for x,i in zip(range(len(fluen)+len(bstart)-1, len(y)),range(0, len(bstart-1))):
        globals()['a%s' % i] = y[x]

# define yerr as variance terms (errors) for our data parameters (listed in same order as for y)
# *note that we have to enter three time errors for the code to work however in reality the error should be the same for all of them (st0, st2 and st3 are really dummy parameters)

    for x,i in zip([ x for x in range(0, len(bstart)-1) if x != ref_ind], [i for i in range(0, len(bstart)-1) if i != ref_ind]):
        globals()['st%s' % i] = yerr[x]
    for x,i in zip(range(len(bstart)-1, len(fluen)+len(bstart)-1),range(0,len(bstart))):
        globals()['sEb%s' % i] = yerr[x]
    for x,i in zip(range(len(fluen)+len(bstart)-1, len(y)),range(0, len(bstart-1))):
        globals()['sa%s' % i] = yerr[x]

# define theta = model parameters, which we define priors for

    X, Z, Q_b, f_a, f_E, r1, r2, r3, mass, radius = theta

# Instead of treating s_t as a parameter, we just hardwire it here

    s_t = 10./1440.

# call model from IDL code defined as modeldata(base, z, x, r1, r2 ,r3)
    model, valid = runmodel(theta,y,tref,bstart,pflux, pfluxe, tobs,numburstssim,ref_ind,gti_checking)

    if not valid:
        return -np.inf, model

# multiplying by scaling factors to match with the data
    model[len(bstart)-1:len(fluen)+len(bstart)-1] *= r3
    model[len(fluen)+len(bstart)-1:len(y)] *= r2

# To simplify final likelihood expression we define inv_sigma2 for each data parameter that describe the error.
# The variance (eg sEb0) is underestimated by some fractional amount, f, for each set of parameters.

    sEb = yerr[len(bstart)-1:len(fluen)+len(bstart)-1]
    sa = yerr[len(fluen)+len(bstart)-1:len(yerr)]

    inv_sigma2 = []
    for i in range (0,len(bstart)-1):
        inv_sigma2.append(1.0/(s_t**2))
    for i in range(0,len(bstart)):
        inv_sigma2.append(1.0/((sEb[i]*f_E)**2))
    for i in range(0,len(bstart)-1):
        inv_sigma2.append(1.0/((sa[i]*f_a)**2))

# Final likelihood expression
    cpts = (y - (model))**2*inv_sigma2 - (np.log(inv_sigma2))

# Test if the result string is defined here. It is, so we return the selected elements of result
# instead of the downselection in model


    base=Q_b
    z=Z
    x=X
    r1=r1
    r2=r2
    r3=r3
    mass=mass
    radius=radius

    model2 = generate_burst_train(base,z,x,r1,r2,r3,mass, radius,bstart,pflux, pfluxe, tobs,numburstssim,run=run,double=double,debug=debug)

# Now also return the model
    return -0.5*np.sum( cpts ), model2

# -------------------------------------------------------------------------#
# This is the expression for the prior on the model parameters. These are pretty uninformative

# Define priors for theta

def lnZprior(z):
    # This beta function for the metallicity prior is from Andy Casey and is an approximation of the metallicity of a mock galaxy
    # at 2.5-4.5 kpc for the location of 1808. Assuming ZCNO = 0.01 is average value.
    from scipy import stats
    import numpy as np
    beta = stats.beta
    ZCNO = 0.01
    #if z == 0.:
    #    return -np.inf
    #else:
    return(np.log(beta(10,3).pdf((np.log10(z/ZCNO)+3)/3.75)/(3.75*np.log(10)*z)))

# make prior function:
# required arrays: mu, sigma and Marray
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

def lnprior(theta):
    import numpy as np
    X, Z, Q_b, f_a, f_E, r1, r2, r3, mass, radius = theta


    if 0.00001 < X < 0.76 and 0.00001 < Z < 0.056 and 0.000001 <= Q_b < 5.0 and 0 < f_a < 100 and 0 < f_E < 100 and   \
      0.005 < r1 < 1.0 and 0.005 < r2 < 3.0 and 0 < r3*1e3 < 1000 and 1.15 < mass < 2.5 \
      and 9 < radius < 17: #upper bound and lower bounds of each parameter defined here. Bounds were found by considering an estimated value for each parameter then giving reasonable limits.
        return 0.0 + lnZprior(Z) + mr_prior(mass, radius)
    return -np.inf

# -------------------------------------------------------------------------#
# Finally we combine the likelihood and prior into the overall lnprob function, called by emcee

# define lnprob, which is the full log-probability function
def lnprob(theta, x, y, yerr):
    import numpy as np
    lp = lnprior(theta)

# Now also returns the model, to accumulate along with the likelihoods

    like, model = lnlike(theta, x, y ,yerr)

    if (not np.isfinite(lp)) or (not np.isfinite(like)):
        return -np.inf, model

    return lp + like, model
# -------------------------------------------------------------- #
# Testing the various functions. Each of these will display the likelihood value, followed by the model-results "blob"
print("lnprior:",lnprior(theta))
print("lnlike:",lnlike(theta, x, y ,yerr))
print("lnprob:",lnprob(theta, x, y, yerr))
print(theta)

# -------------------------------------------------------------------------#
## Running the chain

# This section now defines the initial walker positions and next defines the chain and runs emcee.

pos = [theta + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

print('Ready to run',run_id,'with',len(pos),'walkers')

# Define the sampler for emcee. Threads argument should be set to how many cores your computer has.
#sampler.reset()
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr), threads=threads)

print('Running sampler...')
# Option to run without progress bar:
#run = sampler.run_mcmc(pos, 1000)

#Option to run with progress bar:
width = 60

for i, result in enumerate(sampler.sample(pos, iterations=nsteps)):
    n = int((width+1) * float(i) / nsteps)
    sys.stdout.write("\r[{0}{1}]".format('#' * n, ' ' * (width - n)))
sys.stdout.write("\n")

print('...sampler run '+run_id+' complete')

# -------------------------------------------------------------------------#
#Save multi-threaded sampler:

def __getstate__(sampler):
        sampler_dict = sampler.__dict__.copy()
        del sampler_dict['pool']
        return sampler_dict

sampler_nopool = __getstate__(sampler)
pickle.dump (sampler_nopool, open (run_id+"_chain.p", "wb"))


def __setstate__(self, state):
        self.__dict__.update(state)
# -------------------------------------------------------------------------#
# -------------------------------------------------------------------------#
# # Displaying results
import numpy as np
from chainconsumer import ChainConsumer
ndim = 10

sampler1 = sampler.chain[:,:,0:10]
samples = sampler1[:,:,:].reshape(-1,ndim)

c = ChainConsumer()

cc = c.add_chain(samples,parameters=["X", "Z", "Qb", "fa", "fE", "r1", "r2", "r3", "M", "R"]).configure(smooth=False)


fig = cc.plotter.plot()
fig.savefig('plots/walks{}.png'.format(run_id))

# -------------------------------------------------------------------------#

# finally reset emcee sampler for next runs.

sampler.reset()
