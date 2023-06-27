"""Main module. This has functions that do the sampling, save the chains, and analyse the results."""
## Python packages required:
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import emcee
from astropy.io import ascii
import astropy.units as u
import astropy.constants as const
from matplotlib.ticker import MaxNLocator
from scipy.stats import gaussian_kde
# import scipy.stats as stats
from scipy.interpolate import splrep, BSpline, splint
from chainconsumer import ChainConsumer
from multiprocessing import Pool
import os
import time
from configparser import ConfigParser

import pkg_resources  # part of setuptools
try:
    # this will fail if the package is not pip-installed
    __version__ = pkg_resources.require("beansp")[0].version
except:
    # in which case just record the path
    __version__ = os.getcwd()

# Some constants

BSTART_ERR = 10*u.min # Default uncertainty for burst start times

# -------------------------------------------------------------------------#
## load local  modules
from .settle import settle
from .burstrain import generate_burst_train, next_burst, burstensemble
from .run_model import runmodel
from .get_data import get_obs
from .mrprior import mr_prior
from .run_emcee import runemcee
from .analyse import get_param_uncert_obs, get_param_uncert

# -------------------------------------------------------------------------#


__all__ = (
    "Beans"
)

# Some example prior functions, or you can write your own for input to the code.

# Define priors for theta. mr prior function is located in mrprior.py


def lnZprior(z):
    """
    This beta function for the metallicity prior is from Andy Casey and is an
    approximation of the metallicity of a mock galaxy at 2.5-4.5 kpc for the
    location of 1808. Assuming ZCNO = 0.01 is average value.
    """

    from scipy import stats

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


def prior_1808(theta_in):
    """
    This function implements a simple box prior for all the parameters 
    excluding mass, radius, and the metallicity, which come instead from
    separate functions

    This prior is explicitly intended for use with SAX J1808.4-3653, and
    should only be used with extreme caution in other cases!

    :param theta_in: parameter vector
    """

    X, Z, Q_b, f_a, f_E, r1, r2, r3, mass, radius = theta_in

    # upper bound and lower bounds of each parameter defined here. Bounds were
    # found by considering an estimated value for each parameter then giving
    # reasonable limits.
    if (0.00001 < X < 0.76) and (Z > 0.00001) and\
        (0.000001 <= Q_b < 5.0) and (1 <= f_a < 100) and (1 <= f_E < 100) and \
        (0.005 < r1 < 1.0) and (0.005 < r2 < 3.0) and \
        (0 < r3 * 1e3 < 1000) \
        and (1.15 < mass < 2.5) and (9 < radius < 17):

        return 0.0 + lnZprior(Z) + mr_prior(mass, radius) 
    else:
        return -np.inf


def calc_dist_anisotropy(r1, r2, r3):
    '''
    Calculate distance and inclincation from scaling factors.  These
    functions are taken (more or less) from Goodwin et al. 2019, eq.  18-20, 
    but there seem to be some slight differences

    TODO also infer the inclination, which will depend on the model
    adopted for the anisotropy; because the burst and persistent emission
    anisotropy are free to vary independently, there is no guarantee there
    would be an inclination value consistent with every pair of values for
    a given model

    :param r1: ratio of mdot to bolometric (isotropic) fluence (eq. 12)
    :param r2: ratio of observed to predicted alpha (eq. 13)
    :param r3: ratio of fluence to isotropic burst energy (eq. 14)

    :return: distance, xi_b, xi_p
    '''

    xi_p = np.power( (r1*r2*r3*1e3)/(63.23*0.74816), 0.5)
    xi_b = (0.74816*xi_p)/r2
    distance = 10*np.power((r1/xi_p), 0.5) #kpc

    return distance, xi_b, xi_p


def model_str(model):
    """
    Prints a compressed string representation of the model dict, with
    reduced precision to reduce the size of the record saved to the .h5
    file

    :param model: model dictionary as returned by runmodel

    :return: string representation of the model dict
    """

    return ("{'time': ["+','.join(['{:.4f}'.format(x) for x in model['time']]) \
        +"], 'mdot': ["+','.join(['{:.5f}'.format(x) for x in model['mdot']]) \
        +"], 'alpha': ["+','.join(['{:.3f}'.format(x) for x in model['alpha']]) \
        +"], 'e_b': ["+','.join(['{:.4f}'.format(x) for x in model['e_b']]) \
        +"]}").replace(' ','')


def mean_flux_spline(t1, t2, bean):
    """
    Calculates the mean flux between t1 and t2 from the spline set up at
    the __init__ phase

    :param t1: start time for averaging
    :param t2: end time for averaging
    :param bean: Beans object, from which the remaining parameters are drawn:
      tck_s

    :result: mean flux
    """

    return splint(t1,t2,bean.tck_s)/(t2-t1)


def mean_flux_linear(t1, t2, bean):
    """
    Calculates the mean flux between t1 and t2 from the piecewise linear
    interpolation of tobs,a,b

    :param t1: start time for averaging
    :param t2: end time for averaging
    :param bean: Beans object, from which the remaining parameters are drawn:
      tobs, a, b

    :result: mean flux
    """

    tobs, a, b = bean.tobs, bean.a, bean.b
    na = len(a)

    if len(tobs) != na + 1 or len(b) != na:
        print("** ERROR ** some problem with tobs, a, b arrays")
        return -1
    # i1 is max of where t1 > tobs

    i1 = max(np.where(t1 > tobs))
    if len(i1) == 0:
        i1 = -1
    else:
        i1 = max(i1)

    i2 = max(np.where(t2 > tobs))
    if len(i2) == 0:
        i2 = -1
    else:
        i2 = max(i2)

    sum = 0.0
    if (i1 < 0) and (i2 < 0):

        # Modified this section to just report the flux from the first measurement

        sum += (t2 - t1) * (a[0] + b[0] * tobs[0])
    else:

        # Add the contribution from the start time through to tobs[0], using the
        # gradient between the first pair of observations

        if i1 < 0:
            sum = sum + (tobs[0] - t1) * np.mean(
                [a[0] + b[0] * t1, a[0] + b[0] * tobs[0]]
            )

        # Add the contributions between each pair of observations

        for i in range(max([0, i1]), min([i2, na - 1]) + 1):
            sum = sum + (min([t2, tobs[i + 1]]) - max([t1, tobs[i]])) * np.mean(
                [a[i] + b[i] * max([t1, tobs[i]]), a[i] + b[i] * min([t2, tobs[i + 1]])]
            )

        # Add the contribution for the last section which overlaps with the
        # observation times

        if i1 < na and i2 == na:
            sum = sum + (t2 - tobs[i2]) * np.mean(
                [a[i2 - 1] + b[i2 - 1] * tobs[i2], t2]
            )

        # Now add any contribution *completely* beyond the observations. Previously
        # this gave an exception

        if i1 == na and i2 == na:
            sum = sum + (t2 - t1) * np.mean(
                [a[na - 1] + b[na - 1] * t1, a[na - 1] + b[na - 1] * t2]
            )

    return sum / (t2 - t1)


def get_a_b(pflux, pfluxe, tobs):
    """
    Do piecewise continuous fits to the flux evolution, to determine the
    appropriate parameters for each interval for use by mean_flux_linear:

    :param pflux: persistent flux measurements to interpolate
    :param pfluxe: uncertainty on persistent flux (not used)
    :param tobs: time (midpoint of observation extent) for flux measurement
    :result: a, b arrays for use with mean_flux
    """

    # Now actually calculate the coefficients for the flux fit

    # Linear fit
    ng = len(tobs)

    b0 = np.zeros(ng - 1)
    a0 = np.zeros(ng - 1)

    for i in range(1, ng):
        b0[i - 1] = (pflux[i] - pflux[i - 1]) / (tobs[i] - tobs[i - 1])
        a0[i - 1] = pflux[i - 1] - b0[i - 1] * tobs[i - 1]

    return a0, b0


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

    HAS_CONCORD = False
    try:
        # Required for the distance_limit method
        import concord as cd
        HAS_CONCORD = True
    except:
        pass

    def __init__(self, config_file=None, nwalkers=200, nsteps=100,
                 run_id="test", obsname=None, burstname=None, gtiname=None,
                 theta= (0.44, 0.01, 0.18, 2.1, 3.5, 0.108, 0.90, 0.5, 1.4, 11.2),
                 numburstssim=3, bc=2.21, ref_ind=1, prior=prior_func,
                 threads = 4, test_model=True, restart=False, 
                 interp='linear', smooth=0.02, **kwargs):
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
        :param interp: interpolation mode for the flux; possible values are
          'linear', or 'spline'
        :param smooth: smoothing factor for spline interpolation
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

        self.lnprior = prior

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
            self.obsname = obsname
            self.burstname = burstname
            self.gtiname = gtiname
            self.theta = theta
            self.bc = bc
            self.ref_ind = ref_ind
            self.threads = threads
            self.numburstssim = numburstssim
            # number of bursts observed (redundant; set below after reading the data)
            # self.numburstsobs = numburstsobs

        if self.lnprior(self.theta) == -np.inf:
            print ('** ERROR ** supplied parameter vector is excluded by the prior')
            return

        # below set the parameters which are not part of the config
        # file

        self.restart = restart

        # number of dimensions for the parameter array

        self.ndim = len(self.theta)

        self.gti_checking = self.gtiname is not None

	# determines whether will run as a train of bursts or non-contiguous
	# bursts ("ensemble" mode); previously numerical, default is 1 (True),
	# which means a burst train will be generated; set obsname=None for
	# non-contiguous (ensemble mode) run

        self.train = (self.obsname is not None)

        self.bstart_err = BSTART_ERR.to('d').value

        # Read in all the measurements and set up all the parameters
        # This function now operates on the Beans object directly, and the
        # required attributes:
        # x, y, yerr, tref, bstart, pflux, pfluxe, tobs, fluen, fluen_err, 
        #     st, et
        # are set in that routine
        # bypasses the earlier init function, and instead calls get_obs
        # directly

        get_obs(self)

        self.numburstsobs = len(self.fluen)

        # Set interpolation mode, and define averaging function

        if not hasattr(self, 'interp'):
            self.interp = interp
        assert self.interp in ('linear','spline')
        if self.interp == 'linear':
            self.a, self.b = get_a_b(self.pflux, self.pfluxe, self.tobs)
            self.mean_flux = mean_flux_linear
        else:
            if not hasattr(self, 'smooth'):
                self.smooth = smooth
            self.tck_s = splrep(self.tobs, self.pflux, s=self.smooth)
            self.mean_flux = mean_flux_spline


        # # -------------------------------------------------------------------------#
        # # TEST THE MODEL WORKS
        # # -------------------------------------------------------------------------#
        print("# -------------------------------------------------------------------------#")
        print("Doing Initialisation..")

        if test_model:

            print("Testing the model works..")


            test, valid, test2 = runmodel(self.theta, self, debug=False) # set debug to True for testing
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

        dist, xi_b, xi_p = calc_dist_anisotropy(r1, r2, r3)

        return """#X = {} \\ hydrogen mass fraction
#Z = {} \\ CNO mass fraction
#Q_b = {} \\ base flux [MeV/nucleon]
#M_NS = {} M_sun \\ neutron star mass
#R_NS = {} km \\ neutron star radius
#f_a, f_E = {}, {} \\ systematic error terms for alpha, fluence
#r_1, r_2, r_3 = {}, {}, {} \\ scaling factors to convert predictions
#  (equivalent to d = {:.2f} kpc, xi_b = {:.3f}, xi_p = {:.3f})""".format(
    X, Z, Q_b, mass, radius, f_a, f_E, r1, r2, r3, dist, xi_b, xi_p).replace(
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
           Config.set("data", "interp", str(self.interp))
           if self.interp == 'spline':
               Config.set("data", "smooth", str(self.smooth))

           Config.add_section("emcee")
           Config.set("emcee", "theta", str(self.theta))
           Config.set("emcee", "numburstssim", str(self.numburstssim))
           Config.set("emcee", "prior", str(self.lnprior))
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
        float_params = ('bc', 'smooth')

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
                elif option in float_params:
                    setattr(self, option, config.getfloat(section, option))
                elif option =='prior':
                    function_name = config.get(section, option).split(' ')[1]
                    if function_name != str(self.lnprior).split(' ')[1]:
                        print ('''
** WARNING ** config file lists prior function as {}, but supplied prior is {}
              To fully replicate the previous run you need to specify the same prior using the prior flag on init
'''.format(function_name, self.lnprior))
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

        # define theta = model parameters, which we define priors for

        X, Z, Q_b, f_a, f_E, r1, r2, r3, mass, radius = theta_in

	# call model (function runmodel, in run_model.py) to generate the burst
	# train, or the set of bursts (for "ensemble" mode. In earlier versions
	# the corresponding IDL function was defined as
        # modeldata(base, z, x, r1, r2 ,r3)

        assert np.allclose(y, self.y)
        model, valid, model2 = runmodel( theta_in, self)
        if not valid:
            return -np.inf, model

        # multiplying by scaling factors to match with the data
        # indices for the fluence and the alphas are shifted by one
        # if we're doing a "train" run, hence the int(self.train) (=1 if
        # self.train == True

        ato = int(self.train) # array "train" offset
        # special to trap "unhashable type" error
        # print (model, ato, len(self.bstart), len(self.fluen))
        model[self.numburstsobs-ato:len(self.fluen)+self.numburstsobs-ato] *= r3
        model[len(self.fluen)+self.numburstsobs-ato:] *= r2

	# To simplify final likelihood expression we define inv_sigma2 for each
	# data parameter that describe the error.  The variance (eg sEb0) is
	# underestimated by some fractional amount, f, for each set of
	# parameters.

        err_fac = np.concatenate(( np.full(self.numburstsobs-ato,1.),
            np.full(self.numburstsobs,f_E), np.full(self.numburstsobs-ato,f_a)))
        inv_sigma2 = 1./(yerr*err_fac)**2

        # Final likelihood expression
        cpts = (self.y - (model)) ** 2 * inv_sigma2 - (np.log(inv_sigma2))

        # Test if the result string is defined here. It is, so we return the selected elements of result
        # instead of the downselection in model

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

        # encoding below is so we have a suitable object for including in
        # the blobs (see the dtype specification in runemcee)

        return lp + like, lp, model_str(model).encode('ASCII')



    def plot_model(self, model=None, mdot=True, title=None):
        """
	Display a plot of the model results, for a burst train calculated
	with generate_burst_train Adapted from the example at
        https://matplotlib.org/gallery/api/two_scales.html

        :param model: array of packed model prediction, OR dict giving full
          model results
        :param mdot: flag to show mdot rather than flux (only possible if
          you're passing the full model)
        :param title: add a title, if required
        """

        tobs = self.bstart
        ebobs = self.fluen

        # for the default linear interpolation connect the flux
        # measurements by lines

        ls = '-'
        if self.interp == 'spline':
            ls = ''

        # the ratios are no longer part of the model, so extract them from the param array here

        r1, r2, r3 = self.theta[5:8]

        if model is None:
            test, valid, model = runmodel(self.theta, self, debug=False)

        full_model = False  # Flag to remember whether we're plotting the full model output of
                            # generate burst train or the packed output array
        # if hasattr(model, "time"):
        if type(model) == dict:
            full_model = True
            timepred = model["time"]
            if len(timepred) == 0:
                print ('** ERROR ** no predicted times to show')
                return
            # ebpred = np.array(model["e_b"])*np.array(model["r3"])
            ebpred = np.array(model["e_b"])*np.array(r3)
        else:
            # The other way to return the model is as an array with the burst times, fluences
            # and alphas all together. So unpack those here
            timepred = model[:self.numburstssim+1]
            ebpred = np.array(model[self.numburstssim+int(self.train):self.numburstssim*2+int(self.train)])*r3

        fig, ax1 = plt.subplots()
        # fig.figure(figsize=(10,7))

        flux_color = 'tab:red'
        bursts_color = 'tab:blue'
        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

        if mdot and full_model:
            ax1.set_ylabel('Accretion rate (fraction of Eddington)', color=flux_color)
            if self.train:
                ax1.errorbar(self.tobs, self.pflux*r1, self.pfluxe*r1, 
                    marker='.', ls=ls, color=flux_color, label='mdot')
                if self.interp == 'spline':
                    t = np.arange(min(self.tobs), max(self.tobs), 0.1)
                    ax1.plot(t, BSpline(*self.tck_s)(t)*r1, color=flux_color)
                # show the time of the "reference" burst
                # ax2.axvline(timepred[self.iref], c='k')
                ax2.axvline(self.bstart[self.ref_ind], c='k')
            else:
                ax1.errorbar(self.bstart, self.pflux*r1, self.pfluxe*r1, fmt='.',
                         color=flux_color, label='mdot')
        else:
            ax1.set_ylabel('Persistent flux', color=flux_color)
            if self.train:
                ax1.errorbar(self.tobs, self.pflux, self.pfluxe,  
                    marker = '.', ls=ls, color=flux_color, label = 'pflux')
                if self.interp == 'spline':
                    t = np.arange(min(self.tobs), max(self.tobs), 0.1)
                    ax1.plot(t, BSpline(*self.tck_s)(t)*r1, color=flux_color)
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
            # and the averaged mdot over the burst interval (predicted)
            av_mdot = []
            for i in range(len(timepred)-1):
                av_mdot.append(self.mean_flux(timepred[i], timepred[i+1], self)*r1)
            av_mdot.insert(0, av_mdot[0])
            ax1.step(timepred, av_mdot, where='pre', color=flux_color)
        else:
            ax2.scatter(self.bstart,ebobs, color = 'darkgrey', marker = '.', label='observed', s =200)
            ax2.scatter(self.bstart, ebpred, marker = '*',color=bursts_color,s = 100, label = 'predicted')

        ax2.tick_params(axis='y', labelcolor=bursts_color)

        ax1.set_xlabel("Time (days after outburst start)")
        if title is not None:
            plt.title(title)

        fig.tight_layout()  # otherwise the right y-label is slightly clipped

        fig.legend(loc=1)
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
        print("lnlike:", self.lnlike(self.theta, None, self.y, self.yerr))
        print("lnprob:", self.lnprob(self.theta, None, self.y, self.yerr))
        print("# -------------------------------------------------------------------------#")
        # print(f"The theta parameters will begin at: {self.theta}")
        # print("# -------------------------------------------------------------------------#")
        if plot:
            print("plotting the initial guess.. (you want the predicted bursts to match approximately the observed bursts here)")
            # make plot of observed burst comparison with predicted bursts:
            self.plot_model(title='Initial guess of parameters')
            value = input('Press [RETURN] to continue: ')

        _start = time.time()

        # run the chains and save the output as a h5 file
        sampler = runemcee(self.nwalkers, self.nsteps, self.theta, self.lnprob, self.lnprior, 
            None, self.y, self.yerr, self.run_id, self.restart, self.threads)
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


    def do_analysis(self, options=['autocor','posteriors'], 
                          truths=None, burnin=2000, savefig=True):
        """
        This method is for running standard analysis and displaying the
        results.
        Nothing is returned, but by default the method will create several
        files, labeled by the run_id; drawn from
          {}_autocorrelationtimes.pdf (via plot_autocorr)
          {}_predictedburstscomparison.pdf
          {}chain-plot.pdf
          {}_massradius.pdf
          {}_posteriors.pdf
          {}_parameterconstraints_pred.txt

        TODO: need to reorganise a bit, and add more options

        :param options: array of strings corresponding to various analysis 
            options, listed in the analyses dict below
        :param truths: parameter vector to overplot on (one of the) corner
          plots (TODO: need to check if >1 corner plot are selected
	  +truths, which will likely result in an error due to
          incompatible number of parameters)
        :param burnin: number of steps to discard when plotting the posteriors
        :param savefig: set to True to save figures to .pdf files, False to skip

        :return: none
        """

        # constants:
        c = const.c.to('cm s-1')
        G = const.G.to('cm3 g-1 s-2')

        # list of available analyses

        analyses = {'autocor': 'autocorrelation times as a function of timestep',
                    'chain': 'first 300 iterations of the chains',
                    'posteriors': 'raw posteriors and the input values',
                    'mrcorner': 'corner plot with M, R, g and 1+z',
                    'fig6': 'corner plot with xi_b, xi_p, d, Q_b, Z',
                    'fig8': 'xi_b vs. xi_p and models for comparison',
                    'comparison': 'observed and predicted burst times, fluences' }

        # check the chosen option is one of those implemented

        for option in options:
            if option not in analyses.keys():
                print ('** ERROR ** {} is not an available analysis option; choose from'.format(option))
                for key in analyses.keys():
                    print ('  {}: {}'.format(key, analyses[key]))
                return

        # ---------------------------------------------------------------------#
        # PLOTS
        # ---------------------------------------------------------------------#

        if not hasattr(self, 'reader'):

            print ("Reading in samples...")# to calculate autocorrelation time...")

            # load in sampler:
            self.reader = emcee.backends.HDFBackend(filename=self.run_id+".h5")

            # Read in the full chain to get the number of steps completed
            self.sampler = self.reader.get_chain(flat=False)
            self.nsteps_completed = np.shape(self.sampler)[0]

            print ("... done. Got {} steps completed".format(self.nsteps_completed))

        # moved burnin to be a parameter, so we can pass that from do_run

        if burnin >= self.nsteps_completed*0.9:
            print ('** WARNING ** discarding burnin {} will leave too few steps ({} total), ignoring'.format(burnin, self.nsteps_completed))
            burnin = 0

        # print ("Reading in flattened samples to show posteriors...")
        # samples = self.reader.get_chain(flat=True, discard=burnin)
        self.samples = self.sampler[burnin:,:,:].reshape((-1,10))

        # ---------------------------------------------------------------------#
        if 'autocor' in options:

            # plot autocorrelation times

            if savefig:
                self.plot_autocorr(self.reader, savefile='{}_autocorrelationtimes.pdf'.format(self.run_id))
            else:
                self.plot_autocorr(self.reader, savefile=None)
                print ('Skipping autocorrelation plot save')
            print ("...done")

        #sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, self.lnprob, args=(self.x, self.y, self.yerr), backend=reader)

        # ---------------------------------------------------------------------#
        if 'chain' in options:

            # plot the chains:

            print ("Plotting the chains...")
            labels = ["$X$","$Z$","$Q_b$","$f_a$","$f_E$","$r1$","$r2$","$r3$", "$M$", "$R$"]
            # plt.clf()
            fig, axes = plt.subplots(self.ndim, 1, sharex=True, figsize=(8, 9))

            for i in range(self.ndim):
                axes[i].plot(self.sampler[:,:,i].T, color="k", alpha=0.4)
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

        # ---------------------------------------------------------------------#
        if 'posteriors' in options:

            # Also read in the "flattened" chain, for the posteriors

            # make plot of posterior distributions of your parameters:
            cc = ChainConsumer()
            cc.add_chain(self.samples, parameters=["X", "Z", "Qb", "fa", "fE", "r1", "r2", "r3", "M", "R"])

            if truths is None:
                truths = list(self.theta)

            if savefig:
                cc.plotter.plot(filename=self.run_id+"_posteriors.pdf",
                    figsize="page", truth=truths)
            else:
                fig = cc.plotter.plot(figsize="page", truth=truths)
                fig.show()
            print ("...done")

        # ---------------------------------------------------------------------#
        if (('mrcorner' in options) or ('comparison' in options) \
            or ('fig6' in options) or ('fig8' in options)):
            # maybe can try to record that we've already read these
            # parameters in, like so
            # and (self.burnin_excluded != burnin):

            # and finally read in the model realisations
            # This loop can take a LOOOOOONG time for long runs
            # TODO save this to the Beans object so we only need to read
            # it in once; need to rationalise what arrays are kept etc.

            print ("Reading in and processing blobs, ignoring first {}...".format(burnin))
            blobs = self.reader.get_blobs(flat=True)

            # Get predictions for each model run from the blobs structure:

            time, e_b, alpha, mdot = [], [], [], []
            # X = [] # just for testing; won't be included in future blobs
            for i in range(burnin*self.nwalkers, len(blobs["model"])):
                _model = eval(blobs["model"][i].decode('ASCII', 'replace'))
                time.append(_model['time'])
                e_b.append(_model['e_b'])
                alpha.append(_model['alpha'])
                mdot.append(_model['mdot'])
                # X.append(_model['x_0'][0])
            print ("...done")

            # we don't include these in the blobs anymore, to save space
            # (and they're redundant, as also included in the samples)
            # X = [data[i]['x_0'][0] for i in range(len(data))]
            # Z = [data[i]['z'][0] for i in range(len(data))]
            # base = [data[i]['base'][0] for i in range(len(data))]
            # r1 = np.array([data[i]['r1'][0] for i in range(len(data))])
            # r2 = np.array([data[i]['r2'][0] for i in range(len(data))])
            # r3 = np.array([data[i]['r3'][0] for i in range(len(data))])
            # mass = np.array([data[i]['mass'][0] for i in range(len(data))])
            # radius = np.array([data[i]['radius'][0] for i in range(len(data))])
            # assert np.allclose(X, self.samples[:,0])
            X = self.samples[:,0]
            Z = self.samples[:,1]
            base = self.samples[:,2]
            r1 = self.samples[:,5]
            r2 = self.samples[:,6]
            r3 = self.samples[:,7]
            mass = self.samples[:,8]
            radius = self.samples[:,9]

            # calculate redshift and gravity from mass and radius:
	    # keep the parameters that we're going to calculate limits on
	    # below, dimensionless

            R = np.array(radius)*1e5*u.cm #cgs
            M = np.array(mass)*const.M_sun.to('g') #cgs

	    # ChainConsumer's plot method can't handle Quantity objects,
	    # so we need to convert gravity and redshift back to numpy
	    # arrays here
            redshift = np.power((1 - (2*G*M/(R*c**2))), -0.5).value
            gravity = (M*redshift*G/R**2 / (u.cm/u.s**2)).value #cgs

            # calculate distance and inclincation from scaling factors:

            distance, xib, xip = calc_dist_anisotropy(r1, r2, r3)

            cosi_2 = 1/(2*xip)
            cosi = 0.5/(2*(xip/xib)-1)

	    # to get the parameter middle values and uncertainty use the
	    # functions get_param_uncert_obs and get_param_uncert_pred,
	    # e.g.

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

            Xpred = get_param_uncert(X)
            Zpred = get_param_uncert(Z)
            basepred = get_param_uncert(base)
            dpred = get_param_uncert(distance)
            cosipred = get_param_uncert(cosi)
            xippred = get_param_uncert(xip)
            xibpred = get_param_uncert(xib)
            masspred = get_param_uncert(mass)
            radiuspred = get_param_uncert(radius)
            gravitypred = get_param_uncert(gravity)
            redshiftpred = get_param_uncert(redshift)
            r1pred = get_param_uncert(r1)
            r2pred = get_param_uncert(r2)
            r3pred = get_param_uncert(r3)

            # scale fluences by scaling factor:
            ebpred = np.array(ebpred)*np.array(r3pred[0])
            ebpred_errup = np.array(ebpred_errup)*np.array(r3pred[0])
            ebpred_errlow = np.array(ebpred_errlow)*np.array(r3pred[0])

            # save to text file with columns: paramname, value, upper uncertainty, lower uncertainty

            np.savetxt(f'{self.run_id}_parameterconstraints_pred.txt', (Xpred, Zpred, basepred, dpred, cosipred, xippred, xibpred, masspred, radiuspred,gravitypred, redshiftpred, r1pred, r2pred, r3pred) , header='Xpred, Zpred, basepred, dpred, cosipred, xippred, xibpred, masspred, radiuspred,gravitypred, redshiftpred, r1pred, r2pred, r3pred \n value, upper uncertainty, lower uncertainty')

	    # make plot of posterior distributions of the mass, radius,
	    # surface gravity, and redshift: stack data for input to
	    # chainconsumer:
            mass = mass.ravel()
            radius = radius.ravel()
            gravity = np.array(gravity).ravel()
            redshift = redshift.ravel()

        # ---------------------------------------------------------------------#
        if 'mrcorner' in options:

            mrgr = np.column_stack((mass, radius, gravity, redshift))

            # plot with chainconsumer:
            cc = ChainConsumer()
            cc.add_chain(mrgr, parameters=["M", "R", "g", "1+z"])
            if savefig:
                cc.plotter.plot(filename=self.run_id+"_massradius.pdf",
                    truth=truths, figsize="page")
            else:
                fig = cc.plotter.plot(figsize="page", truth=truths)
                fig.show()

        # ---------------------------------------------------------------------#
        if 'fig6' in options:

            # fig6data = np.column_stack((xip, xib, distance, base, Z, X))
            fig6data = np.column_stack((X, Z, base, distance, xib, xip))

            # plot with chainconsumer:

            cc = ChainConsumer()

            # configure params below copied from Adelle's jupyter notebook
            cc.add_chain(fig6data, parameters=["X", "$Z$", "$Q_b$ (MeV)",
                "$d$ (kpc)", "$\\xi_b$", "$\\xi_p$"])\
                .configure(flip=False, bins=0.7, summary=False, \
                diagonal_tick_labels=False, max_ticks=3, shade=True, \
                shade_alpha=1.0 ,bar_shade=True, tick_font_size='xx-large', \
                label_font_size='xx-large',smooth=True, \
                sigmas=np.linspace(0, 3, 4))
            if savefig:
                cc.plotter.plot(filename=self.run_id+"_fig6.pdf",
                    truth=truths, figsize="page")
            else:
                fig = cc.plotter.plot(figsize="page", truth=truths)
                fig.show()

        # ---------------------------------------------------------------------#
        if 'fig8' in options:

            # here we read in data from the anisotropy models. There's
            # probably a better way to do this, via concord (if it's
            # available)

            counts, ybins, xbins, image = plt.hist2d(np.array(xip), 
                np.array(xib), bins=500, norm=LogNorm(), cmap='OrRd')

            xi_p_model2 = np.arange(0, 2.5, 0.01)
            xi_b_model2 = np.empty(len(xi_p_model2))

            for i in range(0,250):
    
                xi_b_model2[i] = 1./((1./(2*xi_p_model2[i])) + 0.5)
    
            # overplot the various models

            plt.plot(xi_p_model2, xi_b_model2, color = 'black',ls='-', label = 'Fujimoto (1988)')  

            if Beans.HAS_CONCORD:
                # setup dict with list of models, legend labels and linestyles
                he16_models = {'he16_a': ('He & Keek (2016) model A', '--'),
                               'he16_b': ('model B', '-.'),
                               'he16_c': ('model C', (0, (3, 5, 1, 5))), 
                               'he16_d': ('model D', (0, (1, 5))) }

                for model in he16_models.keys():

                    _model = self.cd.diskmodel.load_he16(model)

                    model_theta = _model['col1']
                    model_xid = _model['col2']
                    model_xir = _model['col3']
                    model_xip1 = _model['col4']
                    model_xib1 = model_xid + model_xir

                    model_xib = 1./model_xib1
                    model_xip = 1./model_xip1

                    #modela:
                    plt.plot(model_xip, model_xib, color='darkblue', 
                        ls=he16_models[model][1], label=he16_models[model][0])

            else:
                print ('''
** WARNING ** install concord if you want to overplot model curves
              See https://github.com/outs1der/concord''')

            plt.xlabel(r'$\xi_{\mathrm{p}}$',fontsize='xx-large')
            plt.ylabel(r'$\xi_{\mathrm{b}}$',fontsize='xx-large')

            plt.legend(loc='best',fontsize='large')

            plt.axis([0.,2.1,0.,2.1])

            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)

            if savefig:
                plt.savefig('{}_xipvsxib_models_contourlines.pdf'.format(self.run_id))
            else:
                plt.show()

        # ---------------------------------------------------------------------#
        if 'comparison' in options:

            # make plot of observed burst comparison with predicted bursts:

            plt.figure(figsize=(10,7))

            # plt.scatter(self.bstart, self.fluen, color = 'black', marker = '.', label='Observed', s =200)
            plt.errorbar(self.bstart, self.fluen, yerr=self.fluene, 
                color='black', linestyle='', marker='.', ms=13, label='Observed')
            #plt.scatter(time_pred_35, e_b_pred_35, marker = '*',color='cyan',s = 200, label = '2 M$_{\odot}$, R = 11.2 km')
            if self.train:
                # plt.scatter(timepred[1:], ebpred, marker='*', color='darkgrey', s=100, label='Predicted')
                plt.errorbar(timepred[1:], ebpred, 
                    yerr=[ebpred_errup, ebpred_errlow],
                    xerr=[timepred_errup[1:], timepred_errlow[1:]], 
                    marker='*', ms=11, color='darkgrey', linestyle='', 
                    label='Predicted')
            else:
                plt.scatter(timepred, ebpred, marker='*', color='darkgrey', s=100, label='Predicted')
                plt.errorbar(timepred, ebpred, yerr=[ebpred_errup, ebpred_errlow],xerr=[timepred_errup, timepred_errlow], fmt='.', color='darkgrey')
                plt.errorbar(self.bstart,  self.fluen, fmt='.', color='black')

            plt.xlabel("Time (days after start of outburst)")
            plt.ylabel("Fluence (1e-9 erg/cm$^2$)")
            plt.legend(loc=2)

            if savefig:
                print ('Saving burst comparison plot to {}_predictedburstscomparison.pdf'.format(self.run_id))
                plt.savefig(f'{self.run_id}_predictedburstscomparison.pdf')
            else:
                print ('Skipping burst comparison plot save')
            plt.show()


