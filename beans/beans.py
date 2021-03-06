"""Main module. This has functions that do the sampling, save the chains, and analyse the results."""
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
from burstrain import generate_burst_train, next_burst, get_a_b, mean_flux
from run_model import runmodel
from get_data import get_obs
from mrprior import mr_prior
from get_data import get_obs
from run_emcee import runemcee
from analyse import get_param_uncert_obs, get_param_uncert
from initialise import init

class Beans:

    def __init__(self, ndim=10, nwalkers=200, nsteps=100, run_id="1808/test1", obsname='1808_obs.txt', burstname='1808_bursts.txt', gtiname='1808_gti.txt', theta= (0.44, 0.01, 0.18, 2.1, 3.5, 0.108, 0.90, 0.5, 1.4, 11.2), numburstssim=3, numburstsobs=4, bc=2.21, ref_ind=1, gti_checking=0, threads = 4, restart=False):

        from initialise import init
        from run_model import runmodel
        # Set up initial conditions:

        self.ndim = ndim
        self.nwalkers = nwalkers # nwalkers and nsteps are the number of walkers and number of steps for emcee to do
        self.nsteps = nsteps
        self.run_id = run_id # Where you want output to be saved and under what name
        self.theta = theta # Set starting value for each theta parameter, Recall odering, theta: X, Z, Q_b, f_a, f_E, r1, r2, r3
        self.threads = threads # Number of threads for emcee to use (e.g. number of cores your computer has)
        self.numburstssim = numburstssim # this needs to be an integer value of half the number of bursts you want to simulate. I.e. simulate this many from the reference burst in either direction. Don't forget to account for missed bursts!
        self.numburstsobs = numburstsobs # number of observed bursts in your dataset
        self.ref_ind = ref_ind # Define index of reference burst (should be middle of predicted burst train). This burst will not be simulated but will be used as a reference to predict the other bursts.
        self.gti_checking = gti_checking #Option to turn on gti time checking (1 for on, 0 for off):
        self.obsname = obsname #set location of your observation data files
        self.burstname = burstname #set location of your burst data files
        self.gtiname = gtiname #set location of your gti data files
        self.bc = bc #bolometric correction to apply to your persistent flux (1.0 if they are already bolometric fluxes):
        self.restart = restart #if your run crashed and you would like to restart from a previous run, with run_id above, set this to True

        self.x, self.y, self.yerr, self.tref, self.bstart, self.pflux, self.pfluxe, self.tobs, self.fluen = init(ndim, nwalkers, theta, run_id, threads, numburstssim, numburstsobs, ref_ind, gti_checking, obsname, burstname, gtiname,bc,restart)


        # # -------------------------------------------------------------------------#
        # # TEST THE MODEL WORKS
        # # -------------------------------------------------------------------------#
        test, valid = runmodel(self.theta, self.y, self.tref, self.bstart, self.pflux, self.pfluxe, self.tobs, self.numburstssim, self.ref_ind, self.gti_checking)
        print("# -------------------------------------------------------------------------#")
        print("Doing Initialisation..")
        print("Testing the model works..")
        print("result: ", test, valid)


    def lnlike(self, theta_in, x, y, yerr):

        # define y = "data" parameters

        for x,i in zip([ x for x in range(0, len(self.bstart)-1) if x != self.ref_ind], [i for i in range(0, len(self.bstart)-1) if i != self.ref_ind]):
            globals()['t%s' % i] = self.y[x]
        for x,i in zip(range(len(self.bstart)-1, len(self.fluen)+len(self.bstart)-1),range(0,len(self.bstart))):
            globals()['Eb%s' % i] = self.y[x]
        for x,i in zip(range(len(self.fluen)+len(self.bstart)-1, len(self.y)),range(0, len(self.bstart-1))):
            globals()['a%s' % i] = self.y[x]

    # define yerr as variance terms (errors) for our data parameters (listed in same order as for y)
    # *note that we have to enter three time errors for the code to work however in reality the error should be the same for all of them (st0, st2 and st3 are really dummy parameters)

        for x,i in zip([ x for x in range(0, len(self.bstart)-1) if x != self.ref_ind], [i for i in range(0, len(self.bstart)-1) if i != self.ref_ind]):
            globals()['st%s' % i] = self.yerr[x]
        for x,i in zip(range(len(self.bstart)-1, len(self.fluen)+len(self.bstart)-1),range(0,len(self.bstart))):
            globals()['sEb%s' % i] = self.yerr[x]
        for x,i in zip(range(len(self.fluen)+len(self.bstart)-1, len(self.y)),range(0, len(self.bstart-1))):
            globals()['sa%s' % i] = self.yerr[x]


        # define theta = model parameters, which we define priors for

        X, Z, Q_b, f_a, f_E, r1, r2, r3, mass, radius = theta_in

        # Instead of treating s_t as a parameter, we just hardwire it here

        s_t = 10.0 / 1440.0

        # call model from IDL code defined as modeldata(base, z, x, r1, r2 ,r3)
        model, valid = runmodel(
            theta_in, y, self.tref, self.bstart, self.pflux, self.pfluxe, self.tobs,self. numburstssim, self.ref_ind, self.gti_checking
        )

        if not valid:
            return -np.inf, model

        # multiplying by scaling factors to match with the data
        model[len(self.bstart)-1:len(self.fluen)+len(self.bstart)-1] *= r3
        model[len(self.fluen)+len(self.bstart)-1:len(self.y)] *= r2

    # To simplify final likelihood expression we define inv_sigma2 for each data parameter that describe the error.
    # The variance (eg sEb0) is underestimated by some fractional amount, f, for each set of parameters.

        sEb = yerr[len(self.bstart)-1:len(self.fluen)+len(self.bstart)-1]
        sa = yerr[len(self.fluen)+len(self.bstart)-1:len(self.yerr)]

        inv_sigma2 = []
        for i in range (0,len(self.bstart)-1):
            inv_sigma2.append(1.0/(s_t**2))
        for i in range(0,len(self.bstart)):
            inv_sigma2.append(1.0/((sEb[i]*f_E)**2))
        for i in range(0,len(self.bstart)-1):
            inv_sigma2.append(1.0/((sa[i]*f_a)**2))

        # Final likelihood expression
        cpts = (self.y - (model)) ** 2 * inv_sigma2 - (np.log(inv_sigma2))

        # Test if the result string is defined here. It is, so we return the selected elements of result
        # instead of the downselection in model

        base = Q_b
        z = Z
        x = X
        r1 = r1
        r2 = r2
        r3 = r3
        mass = mass
        radius = radius

        model2 = generate_burst_train(
            base,
            z,
            x,
            r1,
            r2,
            r3,
            mass,
            radius,
            self.bstart,
            self.pflux,
            self.pfluxe,
            self.tobs,
            self.numburstssim,
            self.ref_ind
        )

        #model2 =  np.string_(model2, dtype='S1000')
        model2 = str(model2).encode('ASCII')

        # Now also return the model
        return -0.5 * np.sum(cpts), model2


    # -------------------------------------------------------------------------#
    # This is the expression for the prior on the model parameters.

    # Define priors for theta. mr prior function is located in mrprior.py


    def lnZprior(self, z):
        # This beta function for the metallicity prior is from Andy Casey and is an approximation of the metallicity of a mock galaxy
        # at 2.5-4.5 kpc for the location of 1808. Assuming ZCNO = 0.01 is average value.
        from scipy import stats
        import numpy as np

        beta = stats.beta
        ZCNO = 0.01

        return np.log(
            beta(10, 3).pdf((np.log10(z / ZCNO) + 3) / 3.75) / (3.75 * np.log(10) * z)
        )


    def lnprior(self, theta_in):
        import numpy as np

        X, Z, Q_b, f_a, f_E, r1, r2, r3, mass, radius = theta_in

        if 0.00001 < X < 0.76 and 0.00001 < Z < 0.056 and 0.000001 <= Q_b < 5.0 \
            and 0 < f_a < 100 and 0 < f_E < 100 and 0.005 < r1 < 1.0 and 0.005 < r2 < 3.0 and 0 < r3 * 1e3 < 1000 and 1.15 < mass < 2.5 and 9 < radius < 17:  # upper bound and lower bounds of each parameter defined here. Bounds were found by considering an estimated value for each parameter then giving reasonable limits.
            return 0.0 + self.lnZprior(Z) + mr_prior(mass, radius)
        else:
            return -np.inf


    # -------------------------------------------------------------------------#
    # Finally we combine the likelihood and prior into the overall lnprob function, called by emcee

    # define lnprob, which is the full log-probability function
    def lnprob(self, theta_in, x, y, yerr):
        import numpy as np

        lp = self.lnprior(theta_in)

        # Now also returns the model, to accumulate along with the likelihoods

        like, model = self.lnlike(theta_in, x, y, yerr)

        if (not np.isfinite(lp)) or (not np.isfinite(like)):
            return -np.inf, -np.inf, model

        # we return the logprobability as well as the theta parameters at this point so we can extract results later
        return lp + like, lp, model



    # -------------------------------------------------------------- #

    def do_run(self):
        ## Running the chain
        # we use multiprocessing to speed things up. Emcee parameters are defined in runemcee module.

        print("# -------------------------------------------------------------------------#")
        # Testing the various functions. Each of these will display the likelihood value, followed by the model-results "blob"
        print("Testing the prior and likelihood functions..")
        print("lnprior:", self.lnprior(self.theta))
        print("lnlike:", self.lnlike(self.theta, self.x, self.y, self.yerr))
        print("lnprob:", self.lnprob(self.theta, self.x, self.y, self.yerr))
        print("# -------------------------------------------------------------------------#")
        print(f"The theta parameters will begin at: {self.theta}")
        print("# -------------------------------------------------------------------------#")
        print("plotting the initial guess.. (you want the predicted bursts to match approximately the observed bursts here)")
        # make plot of observed burst comparison with predicted bursts:
        # get the observed bursts for comparison:
        X, Z, Q_b, f_a, f_E, r1, r2, r3, mass, radius = self.theta
        base = Q_b
        z = Z
        x = X
        r1 = r1
        r2 = r2
        r3 = r3
        mass = mass
        radius = radius

        model = generate_burst_train(
            base,
            z,
            x,
            r1,
            r2,
            r3,
            mass,
            radius,
            self.bstart,
            self.pflux,
            self.pfluxe,
            self.tobs,
            self.numburstssim,
            self.ref_ind
        )
        timepred = model["time"]
        ebpred = np.array(model["e_b"])*np.array(model["r3"])

        tobs = self.bstart
        ebobs = self.fluen
        plt.figure(figsize=(10,7))
        plt.scatter(tobs,ebobs, color = 'black', marker = '.', label='Observed', s =200)
        plt.scatter(timepred[1:], ebpred, marker = '*',color='darkgrey',s = 100, label = 'Predicted')
        #plt.errorbar(timepred[1:], ebpred, yerr=[ebpred_errup, ebpred_errlow], xerr=[timepred_errup[1:],timepred_errlow[1:]], fmt='.', color='darkgrey')
        #plt.errorbar(tobs, ebobs, fmt='.',color='black')
        plt.xlabel("Time (days after start of outburst)")
        plt.ylabel("Fluence (1e-9 erg/cm$^2$)")
        plt.title("Initial guess of parameters")
        plt.legend(loc=2)
        plt.show()

        print("# -------------------------------------------------------------------------#")
        print("Beginning sampling..")

        sampler = runemcee(self.nwalkers, self.nsteps, self.ndim, self.theta, self.lnprob, self.x, self.y, self.yerr, self.run_id, self.restart) # this will run the chains and save the output as a h5 file
        print(f"Sampling complete!")
        self.do_analysis()
        #if self.restart == False:
            #sampler.reset()

# -------------------------------------------------------------------------#
# Analyse and display the results:

    def do_analysis(self):
       # run_id = "chains_1808/test1"

        # constants:
        c = 2.9979e10
        G = 6.67428e-8

    # -------------------------------------------------------------------------#
        # load in sampler:
        reader = emcee.backends.HDFBackend(filename=self.run_id+".h5")
        #sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, self.lnprob, args=(self.x, self.y, self.yerr), backend=reader)
        #tau = 20
        tau = reader.get_autocorr_time(tol=0) #using tol=0 means we'll always get an estimate even if it isn't trustworthy.
        burnin = int(2 * np.max(tau))
        thin = int(0.5 * np.min(tau))
        samples=reader.get_chain(flat=True, discard=burnin)
        sampler=reader.get_chain(flat=False)
        blobs = reader.get_blobs(flat=True)
        # samples = reader.get_chain(discard=burnin, flat=True, thin=thin)
        # log_prob_samples = reader.get_log_prob(discard=burnin, flat=True, thin=thin)
        # blobs = reader.get_blobs(discard=burnin, flat=True, thin=thin)

        data = []
        for i in range(len(blobs["model"])):
            data.append(eval(blobs["model"][i].decode('ASCII', 'replace')))

    # -------------------------------------------------------------------------#
        # get the acceptance fraction:
        #accept = reader.acceptance_fraction/nsteps #this will be an array with the acceptance fraction for each walker
        #print(f"The average acceptance fraction of the walkers is: {np.mean(accept)}")

        # get the autocorrelation times:
        # print("burn-in: {0}".format(burnin))
        # print("thin: {0}".format(thin))
        # print("flat chain shape: {0}".format(samples.shape))
        # print("flat log prob shape: {0}".format(log_prob_samples.shape))
        # print("flat log prior shape: {0}".format(log_prior_samples.shape))

        # alternate method of checking if the chains are converged:
        # This code is from https://dfm.io/posts/autocorr/ 

        # get autocorrelation time:

        def next_pow_two(n):
            i = 1
            while i < n:
                i = i << 1
            return i

        def autocorr_func_1d(x, norm=True):
            x = np.atleast_1d(x)
            if len(x.shape) != 1:
                raise ValueError("invalid dimensions for 1D autocorrelation function")
            n = next_pow_two(len(x))

            # Compute the FFT and then (from that) the auto-correlation function
            f = np.fft.fft(x - np.mean(x), n=2*n)
            acf = np.fft.ifft(f * np.conjugate(f))[:len(x)].real
            acf /= 4*n
            
            # Optionally normalize
            if norm:
                acf /= acf[0]

            return acf


        # Automated windowing procedure following Sokal (1989)
        def auto_window(taus, c):
            m = np.arange(len(taus)) < c * taus
            if np.any(m):
                return np.argmin(m)
            return len(taus) - 1

        # Following the suggestion from Goodman & Weare (2010)
        def autocorr_gw2010(y, c=5.0):
            f = autocorr_func_1d(np.mean(y, axis=0))
            taus = 2.0*np.cumsum(f)-1.0
            window = auto_window(taus, c)
            return taus[window]

        def autocorr_new(y, c=5.0):
            f = np.zeros(y.shape[1])
            for yy in y:
                f += autocorr_func_1d(yy)
            f /= len(y)
            taus = 2.0*np.cumsum(f)-1.0
            window = auto_window(taus, c)
            return taus[window]

        # Compute the estimators for a few different chain lengths

        #loop through 10 parameters:
        f = plt.figure(figsize=(8,5))

        param = ["$X$", "$Z$", "$Q_{\mathrm{b}}$", "$f_{\mathrm{a}}$", "$f_{\mathrm{E}}$", "$r{\mathrm{1}}$",\
                "$r{\mathrm{2}}$", "$r{\mathrm{3}}$", "$M$", "$R$"]
        for j in range(10):
            chain = sampler[:, :, j].T
            print(np.shape(sampler))

            N = np.exp(np.linspace(np.log(100), np.log(chain.shape[1]), 10)).astype(int)
            print(N)
            gw2010 = np.empty(len(N))
            new = np.empty(len(N))
            for i, n in enumerate(N):
                gw2010[i] = autocorr_gw2010(chain[:, :n])
                new[i] = autocorr_new(chain[:, :n])

            # Plot the comparisons
            #plt.loglog(N, gw2010, "o-", label="G\&W 2010")
            plt.loglog(N, new, "o-", label=f"{param[j]}")
            plt.loglog(N, gw2010, "o-", label=None, color='grey')
            ylim = plt.gca().get_ylim()
            
            #plt.ylim(ylim)
        plt.xlabel("Number of samples, $N$", fontsize='xx-large')
        plt.ylabel(r"$\tau$ estimates",fontsize='xx-large')

            
        plt.plot(N, np.array(N)/50.0, "--k")# label=r"$\tau = N/50$")
        plt.legend(fontsize='large',loc='best',ncol=2) #bbox_to_anchor=(0.99, 1.02)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.savefig('{}_autocorrelationtimes.pdf'.format(self.run_id))
        plt.show()


        print(f"The autocorrelation time for each parameter as calculated by emcee is: {tau}")
    # -------------------------------------------------------------------------#
        # Get parameters for each model run from the blobs structure:

        # get each individual parameter:
        time = [data[i]['time'] for i in range(len(data))]
        e_b = [data[i]['e_b'] for i in range(len(data))]
        alpha = [data[i]['alpha'] for i in range(len(data))]
        X = [data[i]['x_0'] for i in range(len(data))]
        Z = [data[i]['z'] for i in range(len(data))]
        base = [data[i]['base'] for i in range(len(data))]
        mdot = [data[i]['mdot'] for i in range(len(data))]
        r1 = np.array([data[i]['r1'] for i in range(len(data))])
        r2 = np.array([data[i]['r2'] for i in range(len(data))])
        r3 = np.array([data[i]['r3'] for i in range(len(data))])
        mass = np.array([data[i]['mass'] for i in range(len(data))])
        radius = np.array([data[i]['radius'] for i in range(len(data))])

        # calculate redshift and gravity from mass and radius:
        R = np.array(radius)*1e5 #cgs
        M = np.array(mass)*1.989e33 #cgs
        redshift = np.power((1 - (2*G*M/(R*c**2))), -0.5)
        gravity = M*redshift*G/R**2 #cgs

        # calculate distance and inclincation from scaling factors:
        r1 = np.array(r1)
        r2 = np.array(r2)
        r3 = np.array(r3)
        print(np.min(r1))
        print(np.min(r2))
        print(np.min(r3))
        print(np.min(mass))
        print(np.min(X))
        
        sqrt = (r1*r2*r3*1e3)/(63.23*0.74816)
        xip = np.power(sqrt, 0.5)
        xib = (0.74816*xip)/r2
        distance = 10*np.power((r1/xip), 0.5) #kpc
        cosi_2 = 1/(2*xip)
        cosi = 0.5/(2*(xip/xib)-1)


        # to get the parameter middle values and uncertainty use the functions get_param_uncert_obs and get_param_uncert_pred, e.g.

        #t1, t2, t3, t4, t5, t6, t7 = get_param_uncert_obs1(time, self.numburstssim+1)
        #times = [list(t1), list(t2), list(t3), list(t4), list(t5), list(t6), list(t7)]
        times = get_param_uncert_obs(time, self.numburstssim*2+1)
        timepred = [x[0] for x in times]
        timepred_errup = [x[1] for x in times]
        timepred_errlow = [x[2] for x in times]

        ebs = get_param_uncert_obs(e_b, self.numburstssim*2)
        ebpred = [x[0] for x in ebs]
        ebpred_errup = [x[1] for x in ebs]
        ebpred_errlow = [x[2] for x in ebs]

        alphas = get_param_uncert_obs(alpha, self.numburstssim*2)

        Xpred = np.array(list(get_param_uncert(X))[0])
        Zpred = np.array(list(get_param_uncert(Z))[0])
        basepred = np.array(list(get_param_uncert(base))[0])
        dpred = np.array(list(get_param_uncert(distance))[0])
        cosipred = np.array(list(get_param_uncert(cosi))[0])
        xippred = np.array(list(get_param_uncert(xip))[0])
        xibpred = np.array(list(get_param_uncert(xib))[0])
        masspred = np.array(list(get_param_uncert(mass))[0])
        radiuspred = np.array(list(get_param_uncert(radius))[0])
        gravitypred = np.array(list(get_param_uncert(gravity))[0])
        redshiftpred = np.array(list(get_param_uncert(redshift))[0])
        r1pred = np.array(list(get_param_uncert(r1))[0])
        r2pred = np.array(list(get_param_uncert(r2))[0])
        r3pred = np.array(list(get_param_uncert(r3))[0])

        # scale fluences by scaling factor:
        ebpred = np.array(ebpred)*np.array(r3pred[0])
        ebpred_errup = np.array(ebpred_errup)*np.array(r3pred[0])
        ebpred_errlow = np.array(ebpred_errlow)*np.array(r3pred[0])

        # save to text file with columns: paramname, value, upper uncertainty, lower uncertainty

        np.savetxt(f'{self.run_id}_parameterconstraints_pred.txt', (Xpred, Zpred, basepred, dpred, cosipred, xippred, xibpred, masspred, radiuspred,gravitypred, redshiftpred, r1pred, r2pred, r3pred) , header='Xpred, Zpred, basepred, dpred, cosipred, xippred, xibpred, masspred, radiuspred,gravitypred, redshiftpred, r1pred, r2pred, r3pred \n value, upper uncertainty, lower uncertainty')

    # -------------------------------------------------------------------------#
    # PLOTS
    # -------------------------------------------------------------------------#

        # make plot of posterior distributions of your parameters:
        c = ChainConsumer()
        c.add_chain(samples, parameters=["X", "Z", "Qb", "fa", "fE", "r1", "r2", "r3", "M", "R"])
        c.plotter.plot(filename=self.run_id+"_posteriors.pdf", figsize="column")

        # make plot of posterior distributions of the mass, radius, surface gravity, and redshift:
        # stack data for input to chainconsumer:
        mass = mass.ravel()
        radius = radius.ravel()
        gravity = gravity.ravel()
        redshift = redshift.ravel()
        mrgr = np.column_stack((mass, radius, gravity, redshift))

        # plot with chainconsumer:
        c = ChainConsumer()
        c.add_chain(mrgr, parameters=["M", "R", "g", "1+z"])
        c.plotter.plot(filename=self.run_id+"_massradius.pdf",figsize="column")

        # make plot of observed burst comparison with predicted bursts:
        # get the observed bursts for comparison:
        tobs = self.bstart
        ebobs = self.fluen

        plt.figure(figsize=(10,7))

        plt.scatter(tobs,ebobs, color = 'black', marker = '.', label='Observed', s =200)
        #plt.scatter(time_pred_35, e_b_pred_35, marker = '*',color='cyan',s = 200, label = '2 M$_{\odot}$, R = 11.2 km')
        plt.scatter(timepred[1:], ebpred, marker = '*',color='darkgrey',s = 100, label = 'Predicted')
        #plt.scatter(time_pred_18, e_b_pred_18, marker = '*',color='orange',s = 200, label = '1.4 M$_{\odot}$, R = 10 km')

        plt.errorbar(timepred[1:], ebpred, yerr=[ebpred_errup, ebpred_errlow], xerr=[timepred_errup[1:],timepred_errlow[1:]], fmt='.', color='darkgrey')
        plt.errorbar(tobs, ebobs, fmt='.',color='black')

        plt.xlabel("Time (days after start of outburst)")
        plt.ylabel("Fluence (1e-9 erg/cm$^2$)")
        plt.legend(loc=2)

        plt.savefig(f'{self.run_id}_predictedburstscomparison.pdf')
        plt.show()


        # plot the chains:
        ndim = 10

        labels = ["$X$","$Z$","$Q_b$","$f_a$","$f_E$","$r1$","$r2$","$r3$", "$M$", "$R$"]
        plt.clf()
        fig, axes = plt.subplots(ndim, 1, sharex=True, figsize=(8, 9))

        for i in range(ndim):
            axes[i].plot(sampler[:,:,i].T, color="k", alpha=0.4)
            axes[i].yaxis.set_major_locator(MaxNLocator(5))
            axes[i].set_ylabel(labels[i])

        axes[ndim-1].set_xlabel("step number")
        plt.tight_layout(h_pad=0.0)
        plt.savefig(self.run_id+'chain-plot.pdf')
        plt.show()
           
