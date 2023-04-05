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
import astropy.units as u
import astropy.constants as const
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
from configparser import ConfigParser

import pkg_resources  # part of setuptools
try:
    # this will fail if the package is not pip-installed
    __version__ = pkg_resources.require("beans")[0].version
except:
    # in which case just record the path
    __version__ = os.getcwd()

try:
    # Required for the distance_limit method
    import concord as cd
except:
    pass

# -------------------------------------------------------------------------#
## load local  modules
from .settle import settle
from .burstrain import generate_burst_train, next_burst, get_a_b, mean_flux, burstensemble
from .run_model import runmodel
from .get_data import get_obs
from .mrprior import mr_prior
from .get_data import get_obs
from .run_emcee import runemcee
from .analyse import get_param_uncert_obs, get_param_uncert
from .initialise import init

# -------------------------------------------------------------------------#
# Some example prior functions, or you can write your own for input to the code.

# Define priors for theta. mr prior function is located in mrprior.py


def lnZprior(z):
    """
    This beta function for the metallicity prior is from Andy Casey and is an
    approximation of the metallicity of a mock galaxy at 2.5-4.5 kpc for the
    location of 1808. Assuming ZCNO = 0.01 is average value.
    """

    from scipy import stats
    import numpy as np

    beta = stats.beta
    ZCNO = 0.01

    return np.log(
        beta(10, 3).pdf((np.log10(z / ZCNO) + 3) / 3.75) / (3.75 * np.log(10) * z)
    )


def prior_func(theta_in):
    """
    This function implements a simple box prior for all the parameters 
    excluding mass and radius, which comes instead from a separate mr_prior
    function

    :param theta_in: parameter vector
    """

    import numpy as np

    X, Z, Q_b, f_a, f_E, r1, r2, r3, mass, radius = theta_in

    # upper bound and lower bounds of each parameter defined here. Bounds were
    # found by considering an estimated value for each parameter then giving
    # reasonable limits.
    if (0.00001 < X < 0.76) and (0.00001 < Z < 0.056) and \
        (0.000001 <= Q_b < 5.0) and (1 <= f_a < 100) and (1 <= f_E < 100) and \
        (0.005 < r1 < 1.0) and (0.005 < r2 < 3.0) and \
        (0 < r3 * 1e3 < 1000) \
        and (1.15 < mass < 2.5) and (9 < radius < 17):
        #return 0.0 + lnZprior(Z) + mr_prior(mass, radius) #use this option for 1808 prior
        return 0.0 + mr_prior(mass, radius)
    else:
        return -np.inf


class Beans:
    """
    The main object class that includes the basic functionality required for
    beans. The code will read in burst (and observation) data and attempt to
    simulate bursts to match the observed burst properties. There are two
    principle modes; the original function generates a "train" of individual
    bursts, observed (for example) during a transient outburst, as for the 
    original application to the 2002 outburst of SAX J1808.4-3658, observed
    with RXTE/PCA. The alternative is to match to a set of non-contiguous
    bursts ("ensemble" mode)
    """

    def __init__(self, config_file=None, nwalkers=200, nsteps=100,
                 run_id="test", obsname=None, burstname=None, gtiname=None,
                 theta= (0.44, 0.01, 0.18, 2.1, 3.5, 0.108, 0.90, 0.5, 1.4, 11.2),
                 numburstssim=3, bc=2.21, ref_ind=1, prior=prior_func,
                 threads = 4, test_model=True, restart=False, **kwargs):
        """
        Initialise a Beans object

        :param config_file: file to read in configuration parameters (in
          which case the keyword params below are ignored)
        :param nwalkers: number of walkers for the emcee run
        :param nsteps: number of MCMC steps to run
        :param run_id: string identifier for the run, used to label all the
          result files, and where you want output to be saved
        :param obsname: name of observation file, which includes the flux
          history, from which the mdot is estimated to generate to generate
          the burst train (set obsname=None for a non-contiguous, or "ensemble"
          mode run)
        :param burstname: name of the burst data file, listing the bursts
        :param gtiname: name of the GTI file, set to None to turn off checking
        :param theta: initial centroid values for walker model parameters, with
          X, Z, Q_b, f_a, f_E, r1, r2, r3, mass & radius
        :param numburstssim: number of bursts to simulate, for the "train" mode,
          both earlier and later than the reference burst; i.e. set to half
          the total number of bursts you want to simulate. Don't forget to
          account for missed bursts!
        :param bc: bolometric correction to adopt for the flux history (set
          to 1.0 if the fluxes are already bolometric):
        :param ref_ind: rank of "index" burst, against which the other burst
          times are relative to. For the "train" mode, should be around the
          middle of the predicted burst train. This burst will not be
          simulated but will be used as a reference to predict the other bursts.
        :param prior: prior function to use
        :param threads: number of threads for emcee to use (e.g. number of
          cores your computer has). Set to None to use all available
        :param test_model: flag to test the model during the setup process
        :param restart: set to True to continue a previously interrupted run
        :result: Beans object including all the required data
        """

        # Some housekeeping

        if 'ndim' in kwargs.keys():
            print ('** WARNING ** parameter ndim is redundant (ignored), setting from len of param array')
        if 'numburstsobs' in kwargs.keys():
            print ('** WARNING ** parameter numburstsobs is redundant (ignored), setting from len of burst data')
        if 'gti_checking' in kwargs.keys():
            print ('** WARNING ** parameter gti_checking is redundant (ignored), setting from value of gtiname param')

        # MCU: this is a solution to load packaged data compatible
        # with python 3.6 - 3.10
        # however for python 3.10 and later this might become deprecated
        # alternative method is "from importlib.resources import files"
        # see https://setuptools.pypa.io/en/latest/userguide/datafiles.html

        data_path = os.path.join(os.path.dirname(__file__), 'data')
        print("data_path = " + data_path)

        # Only want to set the default values if both obsname and burstname
        # are not set (indicating a default run). This because setting 
        # obsname=None is also how we indicate an "ensemble" mode run
        if (obsname is None) & (burstname is None):
            obsname = os.path.join(data_path, '1808_obs.txt')

        if run_id is None:
            run_id = os.path.join(data_path, '1808/test1')

        if burstname is None:
            burstname = os.path.join(data_path, '1808_bursts.txt')

        # Set up initial conditions:

        if config_file is not None:
            if not os.path.exists(config_file):
                print ('** ERROR ** config file not found, applying keywords')
            print ('Reading run params from {} ...'.format(config_file))
            self.read_config(config_file)
            print ("...done")
        else:
            # apply the keyword values or defaults
            self.nwalkers = nwalkers
            self.nsteps = nsteps
            self.run_id = run_id
            self.theta = theta
            self.threads = threads
            self.numburstssim = numburstssim
            # number of bursts observed (redundant; set below after reading the data)
            # self.numburstsobs = numburstsobs
            self.ref_ind = ref_ind
            self.obsname = obsname
            self.burstname = burstname
            self.gtiname = gtiname
            self.bc = bc

        self.lnprior = prior
        self.restart = restart

        # number of dimensions for the parameter array
        # self.ndim = ndim
        self.ndim = len(theta)

        # self.gti_checking = gti_checking
        self.gti_checking = gtiname is not None

	# determines whether will run as a train of bursts or non-contiguous
	# bursts ("ensemble" mode); previously numerical, default is 1 (True),
	# which means a burst train will be generated; set obsname=None for
	# non-contiguous (ensemble mode) run

        self.train = (obsname is not None)

        # Read in all the measurements and set up all the parameters

        self.x, self.y, self.yerr, self.tref, self.bstart, self.pflux, \
            self.pfluxe, self.tobs, self.fluen, self.st, self.et = init(
            self.ref_ind, self.gti_checking, self.obsname, self.burstname,
            self.gtiname, self.bc)
        self.numburstsobs = len(self.fluen)
        print(self.st, self.et)


        # # -------------------------------------------------------------------------#
        # # TEST THE MODEL WORKS
        # # -------------------------------------------------------------------------#
        print("# -------------------------------------------------------------------------#")
        print("Doing Initialisation..")

        if test_model:

            print("Testing the model works..")


            test, valid, test2 = runmodel(self.theta, self.y, self.tref, self.bstart,
                                   self.pflux, self.pfluxe, self.tobs, self.numburstssim,self.numburstsobs, self.ref_ind,
                                   self.gti_checking, self.train,self.st, self.et,
                                   debug=False) # set debug to True for testing
            print("result: ", test, valid)

            # MCU Note: commented out - no interactive windows for automated testing
            # self.plot_model(test2)


    def __str__(self):
        """
        Show the parameters that the code has been intialised with
        For restart runs could include the number of steps that has
        already been done
        """

        mode = ('ensemble', 'train')
        restart = ('', ', resuming')
        return """== beans dataset =============================================================
See https://beans-7.readthedocs.io

Run ID: {}
Observation data file: {}
  bolometric correction: {}
GTI data file: {}
Burst data file: {}
  comprising {} observed bursts, ref. to #{}
No. of bursts to simulate: {} ({} mode)
  with {} walkers, {} steps, {} threads{}
Initial parameters:
{} 
==============================================================================""".format(self.run_id, self.obsname, self.bc, self.gtiname, self.burstname,
            self.numburstsobs, self.ref_ind,
            self.train+self.numburstssim*(1+self.train), mode[self.train],
            self.nwalkers, self.nsteps,
            self.threads, restart[int(self.restart)], 
            self.theta_table(self.theta, indent=2) )


    def theta_table(self, theta, indent=0):
        """
        Format the run parameter vector as a table
        Could include the errors for a neatly formatted way to present
        results

        :param theta: the model parameter tuple
        :param indent: number of characters to indent the string from the left
        """

        X, Z, Q_b, f_a, f_E, r1, r2, r3, mass, radius = theta

        return """#X = {} \ hydrogen mass fraction
#Z = {} \ CNO mass fraction
#Q_b = {} \ base flux [MeV/nucleon]
#M_NS = {} M_sun \ neutron star mass
#R_NS = {} km \ neutron star radius
#f_a, f_E = {}, {} \ systematic error terms for alpha, fluence
#r_1, r_2, r_3 = {}, {}, {} \ scaling factors to convert predictions""".format(
            X, Z, Q_b, mass, radius, f_a, f_E, r1, r2, r3).replace(
            '#',' '*indent)


    def save_config(self, file=None, clobber=True):
        """
        Routine to write all the configuration parameters to a file, as a
        record of the run; but also to more easily replicate or continue
        a run

        :param file: name of file to save the config as. If None, the run_id
          will be used as a prefix
        :param clobber: set to True to overwrite any existing file
        """

        if file is None:
            file = self.run_id+'.ini'

        if (not clobber) and (os.path.isfile(file)):
            print ('** ERROR ** config file already exists, skipping write')
            return
        else:
           cfgfile = open(file, "w")

           Config = ConfigParser()
           Config.add_section("beans")
           Config.set("beans", "run_id", self.run_id)
           Config.set("beans", "version", __version__)

           Config.add_section("data")
           Config.set("data", "obsname", str(self.obsname))
           Config.set("data", "burstname", self.burstname)
           Config.set("data", "ref_ind", str(self.ref_ind))
           Config.set("data", "gtiname", str(self.gtiname))
           Config.set("data", "bc", str(self.bc))

           Config.add_section("emcee")
           Config.set("emcee", "theta", str(self.theta))
           Config.set("emcee", "numburstssim", str(self.numburstssim))
           # Config.set("emcee", "prior", str(self.lnprior))
           Config.set("emcee", "nwalkers", str(self.nwalkers))
           Config.set("emcee", "nsteps", str(self.nsteps))
           Config.set("emcee", "threads", str(self.threads))

           Config.write(cfgfile)
           cfgfile.close()


    def read_config(self, file=None):
        """
        Routine to read all the configuration parameters from a file, to
        more easily replicate or continue a run

        :param file: name of file to read the config from.
        """

        if file is None:
            data_path = os.path.join(os.path.dirname(__file__), 'data')
            run_id = os.path.join(data_path, 'beans.ini')

        int_params = ('ref_ind','numburstssim','nwalkers','nsteps','threads')

        if not os.path.isfile(file):
            print ('** ERROR ** config file not found')
            return

        config = ConfigParser(allow_no_value=True)
        config.read(file)
        
        # Loop over sections, attributes

        for section in config.sections():
            # print("Section: %s" % section)
            for option in config.options(section):
                # print(
                #     "x %s:::%s:::%s"
                #     % (option, config.get(section, option), str(type(option))))
                if option == 'theta':
                    setattr(self, option, 
                        tuple(map(float, config.get(section, option)[1:-1].split(', '))))
                elif option == 'bc':
                    setattr(self, option, config.getfloat(section, option))
                elif option in int_params:
                    setattr(self, option, config.getint(section, option))
                else:
                    # string options (including "None")
                    _value = config.get(section, option)
                    if _value == 'None':
                        setattr(self, option, None)
                    else:
                        setattr(self, option, _value)


    def lnlike(self, theta_in, x, y, yerr):
        """
        Calculate the "model" likelihood for the current walker position
        Calls runmodel which actually runs the model, either generating a 
        burst train, or a set of runs for "ensemble" mode. Then extracts the
        relevant model outputs and calculates the likelihood.
        Includes an *additional* call to generate_burst_train/burstensemble
        to get the model, so as to be able to return it to the calling function;
        this step is totally redundant

        :param theta_in: model parameter tuple, with X, Z, Q_b, f_a, f_E, r1,
          r2, r3, mass & radius
        :param x: the "independent" variable, passed to lnlike
        :param y: the "dependent" variable (i.e. measurements), passed to lnlike
        :param yerr: erorr estimates on y
        :return: likelihood, model result array
        """

        # define y = "data" parameters
        # I think these "globals" are not used; the only other reference I 
        # can see is in run_model, which is commented out. So I think
        # TODO these siz for loops can be deleted - dkg

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

	      # call model (function runmodel, in run_model.py) to generate the burst
	      # train, or the set of bursts (for "ensemble" mode. In earlier versions
	      # the corresponding IDL function was defined as
        # modeldata(base, z, x, r1, r2 ,r3)

        model, valid, model2 = runmodel(
            theta_in, y, self.tref, self.bstart, self.pflux, self.pfluxe, self.tobs,self.numburstssim,self.numburstsobs, self.ref_ind, self.gti_checking,self.train,
             self.st, self.et
        )
        if not valid:
            return -np.inf, model

        # multiplying by scaling factors to match with the data
        # indices for the fluence and the alphas are shifted by one
        # if we're doing a "train" run, hence the int(self.train) (=1 if
        # self.train == True

        ato = int(self.train) # array "train" offset
        # special to trap "unhashable type" error
        # print (model, ato, len(self.bstart), len(self.fluen))
        model[len(self.bstart)-ato:len(self.fluen)+len(self.bstart)-ato] *= r3
        model[len(self.fluen)+len(self.bstart)-ato:] *= r2

	# To simplify final likelihood expression we define inv_sigma2 for each
	# data parameter that describe the error.  The variance (eg sEb0) is
	# underestimated by some fractional amount, f, for each set of
	# parameters.
        # TODO: assembling inv_sigma2 can probably all be done in one line

        sEb = yerr[len(self.bstart)-ato:len(self.fluen)+len(self.bstart)-ato]
        sa = yerr[len(self.fluen)+len(self.bstart)-ato:]

        inv_sigma2 = []
        if self.train:
            for i in range (0,len(self.bstart)-1):
                inv_sigma2.append(1.0/(s_t**2))
        else:
            for i in range (0,len(self.bstart)):
                inv_sigma2.append(1.0/(yerr[i]**2))
        for i in range(0,len(self.bstart)):
            inv_sigma2.append(1.0/((sEb[i]*f_E)**2))
        for i in range(0,len(self.bstart)-ato):
            inv_sigma2.append(1.0/((sa[i]*f_a)**2))

        # Final likelihood expression
        cpts = (self.y - (model)) ** 2 * inv_sigma2 - (np.log(inv_sigma2))

        # Test if the result string is defined here. It is, so we return the selected elements of result
        # instead of the downselection in model

        model2 = str(model2).encode('ASCII')


        # Now also return the model
        return -0.5 * np.sum(cpts), model2


    def lnprob(self, theta_in, x, y, yerr):
        """
        The full log-probability function incorporating the priors (via 
        lnprior), and and model likelihood (via lnlike), that is passed to
        runemcee when creating the sampler (in the do_run method).

        :param theta_in:
        :param x: the "independent" variable, passed to lnlike
        :param y: the "dependent" variable (i.e. measurements), passed to lnlike
        :param yerr: erorr estimates on y
        :return: total (prior+model) likelihood, prior likelihood, model array
          (from lnlike)
        """
   
        lp = self.lnprior(theta_in)
        # Check if the parameters are consistent with the prior, and skip
        # the model run it if not
        if (not np.isfinite(lp)):
            return -np.inf, -np.inf, None

        # Now also returns the model, to accumulate along with the likelihoods

        like, model = self.lnlike(theta_in, x, y, yerr)

        if (not np.isfinite(like)):
            return -np.inf, -np.inf, model

        # we return the logprobability as well as the theta parameters at this point so we can extract results later
        return lp + like, lp, model


    def plot_model(self, model=None, mdot=True):
        """
        Display a plot of the model results, for a burst train calculated with generate_burst_train
        Adapted from the example at https://matplotlib.org/gallery/api/two_scales.html

        :param model: array of packed model prediction, OR dict giving full
          model results
        :param mdot: flag to show mdot rather than flux (only possible if
          you're passing the full model)
        """
        tobs = self.bstart
        ebobs = self.fluen

        if model is None:
            test, valid, model = runmodel(self.theta, self.y, self.tref,
                self.bstart, self.pflux, self.pfluxe, self.tobs,
                self.numburstssim, self.numburstsobs, self.ref_ind,
                self.gti_checking, self.train,self.st, self.et,
                debug=False)

        full_model = False  # Flag to remember whether we're plotting the full model output of
                            # generate burst train or the packed output array
        # if hasattr(model, "time"):
        if type(model) == dict:
            full_model = True
            timepred = model["time"]
            if len(timepred) == 0:
                print ('** ERROR ** no predicted times to show')
                return
            ebpred = np.array(model["e_b"])*np.array(model["r3"])
        else:
            # The other way to return the model is as an array with the burst times, fluences
            # and alphas all together. So unpack those here
            timepred = model[:self.numburstssim+1]
            # Don't have access to the r3 value to scale, as we did for ebpred above
            ebpred = np.array(model[self.numburstssim+int(self.train):self.numburstssim*2+int(self.train)])*self.theta[6]

        fig, ax1 = plt.subplots()
        # fig.figure(figsize=(10,7))

        flux_color = 'tab:red'
        bursts_color = 'tab:blue'
        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        if mdot and full_model:
            ax1.set_ylabel('Accretion rate (fraction of Eddington)', color=flux_color)
            if self.train:
                ax1.errorbar(self.tobs, self.pflux*model['r1'], self.pfluxe*model['r1'], marker='.',
                         color=flux_color, label='mdot')
                # show the time of the "reference" burst
                ax2.axvline(timepred[model["iref"]], c='k')
            else:
                ax1.errorbar(self.bstart, self.pflux*model['r1'], self.pfluxe*model['r1'], fmt='.',
                         color=flux_color, label='mdot')
        else:
            ax1.set_ylabel('Persistent flux', color=flux_color)
            if self.train:
                ax1.errorbar(self.tobs, self.pflux, self.pfluxe, marker = '.',color=flux_color,
                         label = 'pflux')
                ax2.scatter(timepred[0], ebpred[0], marker = '*',color=bursts_color,s = 100)
                ax1.set_xlabel("Time (days after start of outburst)")
            else:
                # "ensemble" mode plot vs. epoch, rather than observation time
                ax1.errorbar(self.bstart, self.pflux, self.pfluxe, fmt='.', color=flux_color,
                         label='pflux')
                ax1.set_xlabel("Epoch (MJD)")

        ax1.tick_params(axis='y', labelcolor=flux_color)

        # Plot the bursts here
        ax2.set_ylabel("Fluence (1e-9 erg/cm$^2$)", color=bursts_color)
        if self.train:
            if tobs is not None:
                # Plot the observed bursts, if available
                ax2.scatter(tobs,ebobs, color = 'darkgrey', marker = '.', label='observed', s =200)
            ax2.scatter(timepred[1:], ebpred, marker = '*',color=bursts_color,s = 100, label = 'predicted')
        else:
            ax2.scatter(self.bstart,ebobs, color = 'darkgrey', marker = '.', label='observed', s =200)
            ax2.scatter(self.bstart, ebpred, marker = '*',color=bursts_color,s = 100, label = 'predicted')

        ax2.tick_params(axis='y', labelcolor=bursts_color)

        fig.tight_layout()  # otherwise the right y-label is slightly clipped

        # fig.title("Initial guess of parameters")
        fig.legend(loc=2)
        fig.show()

    # -------------------------------------------------------------- #

    def distance_limit(self, distrange=[1,10], Xmin= 0.01, gridsize=10, numtrains=100, skiplo=True):
        """
        Generate burst trains and test how many of them are consistent with the GTIs
        """

        assert len(distrange) == 2

        X_0, Z, Q_b, f_a, f_E, r1_0, r2, r3, mass, radius = self.theta

        # Loop over the distance values

        distance = np.linspace(distrange[0], distrange[1], gridsize)
        X = np.linspace(Xmin, 0.7, gridsize)
        idist = cd.iso_dist(nsamp=gridsize)
        xi_b, xi_p = cd.anisotropy(idist)

        # initialise the array to keep to grid results

        result = np.zeros((gridsize, gridsize))

        for i, _distance in enumerate(distance):

            # For each distance, loop over the h-fraction range
            for j, _X in enumerate(X):

                # initialise result array
                valid_array = []

                print('Drawing one of {} values for xi_p:'.format(len(xi_p.distribution)))
                # For each distance and h-fraction, draw a bunch of inclination values and xi_p's:
                for _xi_p in xi_p.distribution:

                    r1 = _xi_p.value*(_distance/10)**2
                    theta_1 = (_X, Z, Q_b, f_a, f_E, r1, r2, r3, mass, radius)

                    # Run for just one burst to get the initial interval
                    # Set ref_ind to be zero, will subsequently distribute the start burst times
                    # between up to the simulated interval
                    test, valid, test2 = runmodel(theta_1, self.y, 0.0, self.bstart,
                                           self.pflux, self.pfluxe, self.tobs, 1,1, 0.0,
                                           0, self.train, debug=False)
                    print("result: ", test, valid)
                    # self.plot_model(test)

                    # find the initial range of start times to distribute the bursts over

                    if (test is not None) & (len(test['time']) > 1):
                        # iref = np.argmin(abs(np.array(test['time'])-self.tref))
                        iref = np.argmin(abs(np.array(test['time'])))
                        intvl = test['time'][iref+1]
                    else:
                        # If we can't get an initial burst train, distribute the start times
                        # over the entire outburst
                        intvl = max(self.tobs)
                    # This piece of code is probably redundant now that next_burst checks for
                    # predicted burst times beyond tobs
                    if (intvl > 2.*max(self.tobs)) and skiplo:
                        print ('Skipping trials with d={:.4f} and X={:.4f} due to excessive tdel {} >> {}'.format(
                            _distance, _X, intvl, max(self.tobs) ))
                        valid_array.append(True)
                    else:
                        # intvl = min([test['time'][iref+1], max(self.tobs)])

                        trials = np.linspace(0., intvl, numtrains)
                        print ('Simulating {} trains with reference burst distributed within (0, {:.3f})'.format(numtrains, intvl))

                        for trial in trials:

                            print ('trial {:.4f}:'.format(trial))
                            # Set nburstssim to 100 below, just need to make sure it's sufficient to cover
                            # the whole outburst. Replace ref_ind with trial, as the starting burst time
                            # (ref_ind is meaningless if there's no bursts)
                            test, valid, test2 = runmodel(theta_1, self.y, 0.0, self.bstart,
                                                   self.pflux, self.pfluxe, self.tobs, 100,100, trial,
                                                   1,self.train, gti_start=self.st, gti_end=self.et, debug=False)
                            # for debugging
                            # self.plot_model(test)
                            # breakpoint()

                            valid_array.append(valid)
                            print ('  valid={}'.format(valid))


                result[i,j] = len(np.where(valid_array)[0])/len(valid_array)
                print ('End of loop for d={}, X={}, % valid is {}'.format(
                    _distance, _X, 100*result[i,j] ))

        breakpoint()

    # -------------------------------------------------------------- #

    def do_run(self, plot=False, analyse=True, burnin=2000):
        """
        This routine runs the chain for as many steps as is specified in
        the init call.  Emcee parameters are defined in runemcee module.
        Have previously used multiprocessing to speed things up, but that's
        not currently active

        :param plot: set to True to do the plot of the initial guess. This
          seems redundant since it's also plotted at the __init__ stage
        :param analyse: set to True to call do_analysis automatically once the
          chains finish
        :param burnin: number of steps to ignore at the start of the run, 
          passed to do_analysis

        :return:
        """

        # Want to avoid overwriting existing log & config files

        if (self.restart is False) and (os.path.exists(self.run_id+'.h5')):
            print ('** ERROR ** run will overwrite existing log file {}, set restart=True to extend'.format(self.run_id+'.h5'))
            return

        if (os.path.exists(self.run_id+'.ini')):
            print ('** WARNING ** run will overwrite existing config file {}'.format(self.run_id+'.ini'))
            value = input('              enter Y[RETURN] to continue: ')
            if (value != 'y') and (value != 'Y'):
                print ('do_run terminated')
                return
        self.save_config(clobber=True)


        print("# -------------------------------------------------------------------------#")
        print (self)
        print("# -------------------------------------------------------------------------#")
        # Testing the various functions. Each of these will display the likelihood value, followed by the model-results "blob"
        print("Testing the prior and likelihood functions..")
        print("lnprior:", self.lnprior(self.theta))
        print("lnlike:", self.lnlike(self.theta, self.x, self.y, self.yerr))
        print("lnprob:", self.lnprob(self.theta, self.x, self.y, self.yerr))
        print("# -------------------------------------------------------------------------#")
        # print(f"The theta parameters will begin at: {self.theta}")
        # print("# -------------------------------------------------------------------------#")
        if plot:
            print("plotting the initial guess.. (you want the predicted bursts to match approximately the observed bursts here)")
            # make plot of observed burst comparison with predicted bursts:
            # TODO: this section can presumably be replaced by the plot_model
            # method, which also produces a plot at this point
            # get the observed bursts for comparison:
            X, Z, Q_b, f_a, f_E, r1, r2, r3, mass, radius = self.theta
            print(self.train)
            if self.train:
                model = generate_burst_train(
                    Q_b,Z,X,r1,r2,r3,mass,radius,self.bstart,self.pflux,self.pfluxe,self.tobs,self.numburstssim,self.ref_ind
                )
            else:
                model = burstensemble(Q_b,X,Z,r1,r2,r3,mass,radius,self.bstart,self.pflux,self.numburstsobs)
            timepred = model["time"]
            ebpred = np.array(model["e_b"])*np.array(model["r3"])
    
            # Display initial model
            tobs = self.bstart
            ebobs = self.fluen
            plt.figure(figsize=(10,7))
            plt.scatter(tobs,ebobs, color = 'black', marker = '.', label='Observed', s =200)
            if self.train:
                plt.scatter(timepred[1:], ebpred, marker = '*',color='darkgrey',s = 100, label = 'Predicted')
            else:
                # No predicted time in "ensemble" mode so we just plot the
                # fluences, predicted and observed, as a function of epoch
                plt.scatter(tobs, ebpred, marker='*', color='darkgrey', s=100, label='Predicted')
            #plt.errorbar(timepred[1:], ebpred, yerr=[ebpred_errup, ebpred_errlow], xerr=[timepred_errup[1:],timepred_errlow[1:]], fmt='.', color='darkgrey')
            #plt.errorbar(tobs, ebobs, fmt='.',color='black')
            plt.xlabel("Time (days after start of outburst)")
            plt.ylabel("Fluence (1e-9 erg/cm$^2$)")
            plt.title("Initial guess of parameters")
            plt.legend(loc=2)
            plt.show()

        print("# -------------------------------------------------------------------------#")
        print("Beginning sampling...")
        _start = time.time()

        # run the chains and save the output as a h5 file
        sampler = runemcee(self.nwalkers, self.nsteps, self.theta, self.lnprob, self.x, self.y, self.yerr, self.run_id, self.restart, self.threads)
        print(f"...sampling complete!")

        _end = time.time()
        print ("Sampling took {0:.1f} seconds".format(_end-_start))

        if analyse:
            self.do_analysis(burnin=burnin)
        #if self.restart == False:
            #sampler.reset()

# -------------------------------------------------------------------------#

    def plot_autocorr(self, reader=None, savefile=None, figsize=(8,5) ):
        """
        This method shows the estimated autocorrelation time for the run
        Extracted from do_analysis

        :param savefile: name of file to save as (skip if None)
        :param figsize: size for the figure, (tuple, in inches)

        :return:
        """

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


        def auto_window(taus, c):
            """
            Automated windowing procedure following Sokal (1989)
            """

            m = np.arange(len(taus)) < c * taus
            if np.any(m):
                return np.argmin(m)
            return len(taus) - 1

        def autocorr_gw2010(y, c=5.0):
            """
            Following the suggestion from Goodman & Weare (2010)
            """

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

        if reader is None:
            # load in sampler:
            reader = emcee.backends.HDFBackend(filename=self.run_id+".h5")

        #tau = 20
        tau = reader.get_autocorr_time(tol=0) #using tol=0 means we'll always get an estimate even if it isn't trustworthy.
        thin = int(0.5 * np.min(tau)) # seems to be ignored - dkg
        print(f"The autocorrelation time for each parameter as calculated by emcee is: {tau}")
        print ("  mean {:.1f}, min {:.1f}, max {:.1f}".format(np.mean(tau), 
          min(tau), max(tau)))

        # alternate method of checking if the chains are converged:
        # This code is from https://dfm.io/posts/autocorr/

        # get autocorrelation time:

        sampler=reader.get_chain(flat=False)

        # loop through 10 parameters and plot the evolution of the
        # autocorrelation time estimate for each

        f = plt.figure(figsize=figsize)

        param = ["$X$", "$Z$", "$Q_{\mathrm{b}}$", "$f_{\mathrm{a}}$", "$f_{\mathrm{E}}$", "$r{\mathrm{1}}$",\
                "$r{\mathrm{2}}$", "$r{\mathrm{3}}$", "$M$", "$R$"]
        for j in range(10):
            chain = sampler[:, :, j].T
            # print(np.shape(sampler))

            N = np.exp(np.linspace(np.log(100), np.log(chain.shape[1]), 10)).astype(int)
            # print(N)
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
        if savefile is not None:
            print ('plot_autocorr: saving figure as {}'.format(savefile))
            plt.savefig(savefile)
        plt.show()


    def do_analysis(self, burnin=2000, savefig=True):
        """
        This method is for running standard analysis and displaying the
        results.
        Nothing is returned, but by default the method will create several
        files, labeled by the run_id:
          {}_autocorrelationtimes.pdf (via plot_autocorr)
          {}_predictedburstscomparison.pdf
          {}chain-plot.pdf
          {}_massradius.pdf
          {}_posteriors.pdf
          {}_parameterconstraints_pred.txt

        TODO: need to reorganise a bit, and add more options

        :param burnin: number of steps to discard when plotting the posteriors
        :param savefig: set to True to save figures to .pdf files, False to skip

        :return: none
        """

        # constants:
        c = const.c.to('cm s-1')
        G = const.G.to('cm3 g-1 s-2')

    # -------------------------------------------------------------------------#

        # plot autocorrelation times

        print ("Reading in samples to calculate autocorrelation time...")

        # load in sampler:
        reader = emcee.backends.HDFBackend(filename=self.run_id+".h5")

        if savefig:
            self.plot_autocorr(reader, savefile='{}_autocorrelationtimes.pdf'.format(self.run_id))
        else:
            self.plot_autocorr(reader, savefile=None)
            print ('Skipping autocorrelation plot save')
        print ("...done")

        #sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, self.lnprob, args=(self.x, self.y, self.yerr), backend=reader)

        # Read in the full chain to get the number of steps completed
        sampler=reader.get_chain(flat=False)
        nsteps_completed = np.shape(sampler)[0]

        # plot the chains:

        print ("Plotting the chains...")
        labels = ["$X$","$Z$","$Q_b$","$f_a$","$f_E$","$r1$","$r2$","$r3$", "$M$", "$R$"]
        plt.clf()
        fig, axes = plt.subplots(self.ndim, 1, sharex=True, figsize=(8, 9))

        for i in range(self.ndim):
            axes[i].plot(sampler[:,:,i].T, color="k", alpha=0.4)
            axes[i].yaxis.set_major_locator(MaxNLocator(5))
            axes[i].set_ylabel(labels[i])

        axes[self.ndim-1].set_xlabel("step number")
        plt.tight_layout(h_pad=0.0)
        if savefig:
            print ('Saving chain plot to {}chain-plot.pdf'.format(self.run_id))
            plt.savefig(self.run_id+'chain-plot.pdf')
        else:
            print ('Skipping chain plot save')

        plt.show()
        print ("...done")

        # moved burnin to be a parameter, so we can pass that from do_run

        if burnin >= nsteps_completed*0.9:
            print ('** WARNING ** discarding burnin {} will leave too few steps ({} total), ignoring'.format(burnin, nsteps_completed))
            burnin = 0

        # Also read in the "flattened" chain, for the posteriors

        print ("Reading in flattened samples to show posteriors...")
        samples=reader.get_chain(flat=True, discard=burnin)

        # make plot of posterior distributions of your parameters:
        cc = ChainConsumer()
        cc.add_chain(samples, parameters=["X", "Z", "Qb", "fa", "fE", "r1", "r2", "r3", "M", "R"])
        cc.plotter.plot(filename=self.run_id+"_posteriors.pdf", figsize="column")
        print ("...done")

        # and finally read in the model realisations
        # This loop can take a LOOOOOONG time for long runs

        print ("Reading in and processing blobs...")
        blobs = reader.get_blobs(flat=True)
        # samples = reader.get_chain(discard=burnin, flat=True, thin=thin)
        # log_prob_samples = reader.get_log_prob(discard=burnin, flat=True, thin=thin)
        # blobs = reader.get_blobs(discard=burnin, flat=True, thin=thin)

        data = []
        for i in range(len(blobs["model"])):
            data.append(eval(blobs["model"][i].decode('ASCII', 'replace')))
        print ("...done")

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
        # keep the parameters that we're going to calculate limits on below, 
        # dimensionless

        R = np.array(radius)*1e5*u.cm #cgs
        M = np.array(mass)*const.M_sun.to('g') #cgs

	# ChainConsumer's plot method can't handle Quantity objects, so we need
	# to convert gravity and redshift back to numpy arrays here
        redshift = np.power((1 - (2*G*M/(R*c**2))), -0.5).value
        gravity = (M*redshift*G/R**2 / (u.cm/u.s**2)).value #cgs

        # calculate distance and inclincation from scaling factors:
        r1 = np.array(r1)
        r2 = np.array(r2)
        r3 = np.array(r3)
        print(np.min(r1))
        print(np.min(r2))
        print(np.min(r3))
        print(np.min(mass))
        print(np.min(X))

        xip = np.power( (r1*r2*r3*1e3)/(63.23*0.74816), 0.5)
        xib = (0.74816*xip)/r2
        distance = 10*np.power((r1/xip), 0.5) #kpc
        cosi_2 = 1/(2*xip)
        cosi = 0.5/(2*(xip/xib)-1)


        # to get the parameter middle values and uncertainty use the functions get_param_uncert_obs and get_param_uncert_pred, e.g.

        #t1, t2, t3, t4, t5, t6, t7 = get_param_uncert_obs1(time, self.numburstssim+1)
        #times = [list(t1), list(t2), list(t3), list(t4), list(t5), list(t6), list(t7)]
        if self.train:
            times = get_param_uncert_obs(time, self.numburstssim*2+1)
        else:
            times = get_param_uncert_obs(time, self.numburstsobs)
        timepred = [x[0] for x in times]
        timepred_errup = [x[1] for x in times]
        timepred_errlow = [x[2] for x in times]

        if self.train:
            ebs = get_param_uncert_obs(e_b, self.numburstssim*2)
        else:
            ebs = get_param_uncert_obs(e_b, self.numburstsobs)
        ebpred = [x[0] for x in ebs]
        ebpred_errup = [x[1] for x in ebs]
        ebpred_errlow = [x[2] for x in ebs]
        if self.train:
            alphas = get_param_uncert_obs(alpha, self.numburstssim*2)
        else:
            alphas = get_param_uncert_obs(alpha, self.numburstssim)
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

        # make plot of posterior distributions of the mass, radius, surface gravity, and redshift:
        # stack data for input to chainconsumer:
        mass = mass.ravel()
        radius = radius.ravel()
        gravity = np.array(gravity).ravel()
        redshift = redshift.ravel()
        mrgr = np.column_stack((mass, radius, gravity, redshift))

        # plot with chainconsumer:
        cc = ChainConsumer()
        cc.add_chain(mrgr, parameters=["M", "R", "g", "1+z"])
        cc.plotter.plot(filename=self.run_id+"_massradius.pdf",figsize="column")

        # make plot of observed burst comparison with predicted bursts:
        # get the observed bursts for comparison:
        tobs = self.bstart
        ebobs = self.fluen

        plt.figure(figsize=(10,7))

        plt.scatter(tobs,ebobs, color = 'black', marker = '.', label='Observed', s =200)
        #plt.scatter(time_pred_35, e_b_pred_35, marker = '*',color='cyan',s = 200, label = '2 M$_{\odot}$, R = 11.2 km')
        if self.train:
            plt.scatter(timepred[1:], ebpred, marker='*', color='darkgrey', s=100, label='Predicted')
            plt.errorbar(timepred[1:], ebpred, yerr=[ebpred_errup, ebpred_errlow],xerr=[timepred_errup[1:], timepred_errlow[1:]], fmt='.', color='darkgrey')
            plt.errorbar(tobs, ebobs, fmt='.', color='black')
        else:
            plt.scatter(timepred, ebpred, marker='*', color='darkgrey', s=100, label='Predicted')
            plt.errorbar(timepred, ebpred, yerr=[ebpred_errup, ebpred_errlow],xerr=[timepred_errup, timepred_errlow], fmt='.', color='darkgrey')
            plt.errorbar(tobs, ebobs, fmt='.', color='black')

        plt.xlabel("Time (days after start of outburst)")
        plt.ylabel("Fluence (1e-9 erg/cm$^2$)")
        plt.legend(loc=2)

        if savefig:
            print ('Saving burst comparison plot to {}_predictedburstscomparison.pdf'.format(self.run_id))
            plt.savefig(f'{self.run_id}_predictedburstscomparison.pdf')
        else:
            print ('Skipping burst comparison plot save')
        plt.show()


