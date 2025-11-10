# ___.                                    ___          ___
# \_ |__   ____ _____    ____   ______   /  / ______   \  \
#  | __ \_/ __ \\__  \  /    \ /  ___/  /  /  \____ \   \  \
#  | \_\ \  ___/ / __ \|   |  \\___ \  (  (   |  |_> >   )  )
#  |___  /\___  >____  /___|  /____  >  \  \  |   __/   /  /
#      \/     \/     \/     \/     \/    \__\ |__|     /__/
#
"""Main module. This has functions that do the sampling, save the chains, and analyse the results."""

## Python packages required:
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import emcee
import astropy.units as u
import astropy.constants as const
from astropy.table import Table, MaskedColumn
from astropy.time import Time
from matplotlib.ticker import MaxNLocator
from scipy.interpolate import splrep, BSpline, splint
from chainconsumer import ChainConsumer
from multiprocessing import cpu_count
from math import erf
import os, sys
import logging
import gzip
import time
from configparser import ConfigParser
import pickle

import pkg_resources  # part of setuptools
try:
    # this will fail if the package is not pip-installed
    __version__ = pkg_resources.require("beansp")[0].version
except:
    # in which case just record the path
    __version__ = os.getcwd()

# Set the default font to Times; this doesn't seem to affect the
# ChainConsumer plots

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times']
plt.rcParams['text.usetex'] = True

# Some constants & standard units

BSTART_ERR = 10*u.min # Default uncertainty for burst start times
M_NS = 1.4 # canonical NS mass [M_sun]
R_NS = 11.2 # canonical NS mass [km]
FLUX_U = 1e-9*u.erg/u.cm**2/u.s
FLUEN_U = 1e-6*u.erg/u.cm**2

# Adopted value of the Eddington luminosity/peak luminosity of PRE bursts
# (preliminary values from Galloway et al. 2026, in prep)

L_EDD = 3.5e38*u.erg/u.s
L_EDD_ERR = 1.6e38*u.erg/u.s

BILBY_OUTPUT = 'bilby_out'

# switch to remember what I'm using in lnprob

LNPROB_USES_LNLIKE_SYS = True

# -------------------------------------------------------------------------#
## load local  modules
from .settle import settle
from .burstrain import punkt_train, generate_burst_train
from .run_model import runmodel, burst_time_match
from .get_data import get_obs
from .mrprior import mr_prior
from .run_emcee import runemcee
from .analyse import get_param_uncert_part, get_param_uncert

# multiepoch_mcmc is required to use the grid_interp model
has_multiepoch_mcmc = True
try:
    # TODO perhaps pull the path from an environment variable instead
    sys.path.append('/Users/duncan/python/multiepoch_mcmc')
    from multiepoch_mcmc import grid_interpolator

    from .grid_interp import grid_interp
except:
    has_multiepoch_mcmc = False


def create_logger():
    """
    Create a logger instance where messages are sent
    See https://docs.python.org/3/library/logging.html
    """
    logger = logging.getLogger(__name__)
    if not logger.handlers: # Check if created before, otherwise a reload will add handlers
        logger.setLevel(logging.INFO)
        handler = logging.StreamHandler()
        formatter = logging.Formatter('%(levelname)s: %(name)s: %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.propagate = False # prevent duplicate messages

    return logger


logger = create_logger()


# -------------------------------------------------------------------------#


__all__ = (
    "Beans"
)

def strmeas(val, err, err_hi=None, mask_str='--'):
    '''
    Function to return nicely (TeX) formatted measurements, with errors
    Adapted (and simplified) from strmeas.pro

    For simplicity we consider only the case where the error is smaller
    than the value. Since the error and the value need to be
    plotted to the same digit, we have essentially just two cases:
    integer and integer, or float and float

    :param val: central value (50th percentile or whatever)
    :param err: if err_hi is not supplied, this is the symmetric error; otherwise, the lower error
    :param err_hi: upper uncertainty

    :returns: formatted string
    '''

    # check for string values
    if (type(val) == str) | (type(val) == np.str_):
        try:
            val, err = float(val), float(err)
            if err_hi is not None:
                err_hi = float(err_hi)
        except:
            # if we have strings that are not floats, just return
            # them
            return val

    # check for masked values
    if np.ma.is_masked(val):
        return mask_str

    eta=1e-20     # Threshold for non-zero measurements
    sym_templ = '${{{}}}\pm{{{}}}$'
    asym_templ = '${{{}}}_{{{{{{{}}}}}}}^{{{{{{{}}}}}}}$'

    # get the number of significant figures of each of the errors
    fudge=0.0222764
    n_lo, n_hi = -1, -1
    if abs(err) > eta:
        n_lo=np.floor(np.log10(abs(err))+fudge)

    if err_hi is None:
        err_hi = err

    if abs(err_hi) > eta:
        n_hi=np.floor(np.log10(abs(err_hi))+fudge)

    lmin = int(min([n_lo, n_hi]))
    if (err < 2.*10.**n_lo) | (err_hi < 2.*10.**n_hi):
        # need an extra significant figure if the leading digit is 1
        lmin -= 1

    if (lmin < 0):

        # Floating point value

        rtnfmt = ':.{}f'.format(-lmin)
        if np.allclose(round(err, -lmin), round(err_hi, -lmin)):
            # symmetric errors
            res = sym_templ.format(rtnfmt, rtnfmt).format(val, err)
        else:
            # asymmetric errors
            errfmt = ':+.{}f'.format(-lmin)
            res = asym_templ.format(rtnfmt, errfmt, errfmt).format(val, -err, err_hi)
    else:

        # Integer value

        rtnfmt=':d'
        # round the quantities and convert to integers
        ival = int(np.floor(val/10.**lmin+0.5)*10.**lmin)
        ierr_lo = int(np.floor(err/10.**lmin+0.5)*10**lmin)
        ierr_hi = int(np.floor(err_hi/10.**lmin+0.5)*10**lmin)
        if np.allclose(ierr_lo, ierr_hi):
            # symmetric errors
            res = sym_templ.format(rtnfmt, rtnfmt).format(ival, ierr_lo)
        else:
            # asymmetric errors
            errfmt = ':+d'
            res = asym_templ.format(rtnfmt, errfmt, errfmt).format(ival, -ierr_lo, ierr_hi)


    return res


# Some example prior functions, or you can write your own for input to the code.

# Define priors for theta. mr prior function is located in mrprior.py


def lnZprior(z):
    """
    This beta function for the metallicity prior is from Andy Casey and is an
    approximation of the metallicity of a mock galaxy at 2.5-4.5 kpc for the
    location of SAX J1808.4-3658. Assuming ZCNO = 0.01 is average value.

    :param z: CNO metallicity (mass fraction)

    :return: prior probability
    """

    from scipy import stats

    beta = stats.beta
    ZCNO = 0.01

    return np.log(
        beta(10, 3).pdf((np.log10(z / ZCNO) + 3) / 3.75) / (3.75 * np.log(10) * z)
    )


def prior_func(theta_in):
    """
    This function is the default prior and implements a simple box prior
    for all the parameters. Notably the systematic error on the burst
    times is limited to 10, which gives a maximum absolute error on the
    times of 10*BSTART_ERR = 100 min

    :param theta_in: parameter vector, with *X*, *Z*, *Q_b*, *d*, *xi_b*,
      *xi_p*, and (optionally) *mass*, *radius* & *f_t*

    :return: prior probability
    """

    X, Z, Q_b, dist, xi_b, xi_p, *extra = theta_in
    mass, radius, f_t = extra+[M_NS, R_NS, 1.0][len(extra):]

    # upper bound and lower bounds of each parameter defined here. Bounds were
    # found by considering an estimated value for each parameter then giving
    # reasonable limits.
    if (0.00001 < X < 0.76) and (0.00001 < Z < 0.056) and \
        (0.000001 <= Q_b < 5.0) and (1 < dist < 20) and \
        (0.01 < xi_b < 2) and (0.01 < xi_p < 2) and \
        (1 <= f_t < 10) and \
        (1.15 < mass < 2.5) and (9 < radius < 17):
        return 0.0
    else:
        return -np.inf


def prior_kepler(theta_in):
    """
    This function applies a more restrictive prior suitable for use with
    the corr_kepler correction function, which is only defined over the
    ranges X in (0.2,0.8), Z in (0.005,0.1)

    Otherwise identical to :func:`beans.prior_func`;

    :param theta_in: parameter vector, with *X*, *Z*, *Q_b*, *d*, *xi_b*,
      *xi_p*, and (optionally) *mass*, *radius* & *f_t*

    :return: prior probability
    """

    X, Z, Q_b, dist, xi_b, xi_p, *extra = theta_in
    mass, radius, f_t = extra+[M_NS, R_NS, 1.0][len(extra):]

    # upper bound and lower bounds of each parameter defined here. Bounds were
    # found by considering an estimated value for each parameter then giving
    # reasonable limits.
    if (0.2 < X < 0.8) and (0.005 < Z < 0.1) and \
        (0.000001 <= Q_b < 5.0) and (1 < dist < 20) and \
        (0.01 < xi_b < 2) and (0.01 < xi_p < 2) and \
        (1 <= f_t < 10) and \
        (1.15 < mass < 2.5) and (9 < radius < 17):
        return 0.0
    else:
        return -np.inf


def prior_1808(theta_in):
    """
    This function implements a simple box prior for all the parameters
    excluding mass and radius, which comes instead from
    :meth:`Beans.mr_prior` and the metallicity, which comes instead from
    :meth:`Beans.lnZprior`

    This prior is explicitly intended for use with SAX J1808.4-3658, and
    should only be used with extreme caution in other cases!

    :param theta_in: parameter vector, with *X*, *Z*, *Q_b*, *d*, *xi_b*,
      *xi_p*, *mass*, *radius*, and (optionally) *f_t*

    :return: prior probability
    """

    X, Z, Q_b, dist, xi_b, xi_p, *extra = theta_in
    mass, radius, f_t = extra+[M_NS, R_NS, 1.0][len(extra):]

    # upper bound and lower bounds of each parameter defined here. Bounds were
    # found by considering an estimated value for each parameter then giving
    # reasonable limits.
    if (0.00001 < X < 0.76) and (Z > 0.000010000001) and \
        (0.000001 <= Q_b < 5.0) and (1 < dist < 20) and \
        (0.01 < xi_b < 2) and (0.01 < xi_p < 2) and \
        (1 <= f_t < 10):

        return 0.0 + lnZprior(Z) + mr_prior(mass, radius)
    else:
        return -np.inf


def prior_grid(theta_in):
    """
    This function is to accompany grid_interp and implements a simple
    box prior for all the parameters, respecting the grid range

    Otherwise identical to :func:`beans.prior_func`;

    :param theta_in: parameter vector, with *X*, *Z*, *Q_b*, *d*, *xi_b*,
      *xi_p*, and (optionally) *mass*, *radius* & *f_t*

    :return: prior probability
    """

    c = const.c.to('cm s-1')
    G = const.G.to('cm3 g-1 s-2')

    X, Z, Q_b, dist, xi_b, xi_p, *extra = theta_in
    mass, radius, f_t = extra+[M_NS, R_NS, 1.0][len(extra):]

    R = radius*1e5*u.cm #cgs
    M = mass*const.M_sun.to('g') #cgs
    redshift = np.power((1 - (2*G*M/(R*c**2))), -0.5).value
    gravity = (M*redshift*G/R**2 / (u.cm/u.s**2)).value / 1e14 #cgs

    # upper bound and lower bounds of each parameter defined here. Bounds are
    # determined by the limits of the grid
    if (0.64 <= X <= 0.76) and (0.0025 <= Z <= 0.03) and \
        (0.0 <= Q_b <= 0.6) and (1 < dist < 20) and \
        (0.01 < xi_b < 2) and (0.01 < xi_p < 2) and \
        (1 <= f_t < 10) and \
        (1.86 < gravity < 3.45):
        return 0.0
    else:
        return -np.inf


def corr_goodwin19(burst, **kwargs):
    """
    This is an example Settle correction function that applies the correction
    of Goodwin et al. (2019), multiplying the recurrence time and burst
    energy by 0.65. With this correction there is no need to modify the
    alpha value.
    The ``**kwargs`` can be used to incorporate other parameters into the
    correction; from settle we pass the base flux F, mdot M, H-fraction X,
    metallicity Z, and neutron star radius R and mass M (none used for
    this example)

    :param burst: 3-element tuple output from settle, with alpha, trec [hr], and burst energy [1e39 erg], all values in the observer frame

    :return: tuple with any necessary corrections performed
    """

    return (burst[0], burst[1]*0.65, burst[2]*0.65)


def corr_kepler(burst, **kwargs):
    """
    This is an example Settle correction function that applies the correction
    to match kepler burst results.
    The ``**kwargs`` can be used to incorporate other parameters into the
    correction; from settle we pass the base flux F, mdot M, H-fraction X,
    metallicity Z, and neutron star radius R and mass M (none used for
    this example)

    :param burst: 3-element tuple output from settle, with alpha, trec [hr], and burst energy [1e39 erg], all values in the observer frame

    :return: tuple with any necessary corrections performed
    """

    #Scale Fluence
    Fluence_coef= [-0.01494263, 0.61500959,  0.13263262]
    E_error = Fluence_coef[0]*burst[2]**2 + Fluence_coef[1]*burst[2] + Fluence_coef[2]

    #Scale trec
    a = burst[1]; b = kwargs['Z']; c = kwargs['X']
    form = [1, a, b, c, a ** 2, a * b, a * c, b ** 2, b * c, c ** 2]
    tdel_coff = [4.00396086e-24,  1.78459388e+00, 6.37416752e+01, -4.06404901e+00,
                2.94290411e-02,  4.52264987e+00, -2.35649145e+00, -4.64376836e+02,
                -2.73559070e+01,  8.08332362e+00]
    intercept = -0.8976578105076352
    t_error = np.dot(tdel_coff,form) + intercept

    # alpha ~ tdel/E_b, thus scale alpha as well
    t_fac = (burst[1] - t_error) / burst[1]
    e_fac = (burst[2] - E_error) / burst[2]
    alp_scl = burst[0] * (1 - t_fac) / (1 - e_fac)

    return (alp_scl, t_error, E_error)


def model_str(model):
    """
    Prints a compressed string representation of the model dict, with
    selected parameters and reduced precision to reduce the size of the
    record saved to the .h5 file

    :param model: model dictionary as returned by runmodel

    :return: string representation of the model dict
    """

    return ("{'time': ["+','.join(['{:.4f}'.format(x) for x in model['time']])
        +"], 'mdot': ["+','.join(['{:.5f}'.format(x) for x in model['mdot']])
        +"], 'alpha': ["+','.join(['{:.3f}'.format(x) for x in model['alpha']])
        +"], 'e_b': ["+','.join(['{:.4f}'.format(x) for x in model['e_b']])
        +"]}").replace(' ','')


def mean_flux_spline(t1, t2, bean):
    """
    Calculates the mean flux between t1 and t2 from the spline set up at
    the __init__ phase

    :param t1: start time for averaging
    :param t2: end time for averaging
    :param bean: :class:`beansp.Beans` object, from which the remaining parameters are drawn:
      tck_s

    :result: mean flux
    """

    return splint(t1,t2,bean.tck_s)/(t2-t1)


def mean_flux_linear(t1, t2, bean):
    """
    Calculates the mean flux between ``t1`` and ``t2`` from the piecewise
    linear interpolation of ``tobs``, ``a``, ``b``

    :param t1: start time for averaging
    :param t2: end time for averaging
    :param bean: :class:`beansp.Beans` object, from which the remaining parameters are drawn:
      tobs, a, b

    :result: mean flux
    """

    tobs, a, b = bean.tobs, bean.a, bean.b
    na = len(a)

    if len(tobs) != na + 1 or len(b) != na:
        logger.error("some problem with tobs, a, b arrays")
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

    HAS_BILBY = False
    try:
        # required for using dynesty sampler
        import bilby
        from .run_bilby import runbilby
        HAS_BILBY = True
    except:
        pass

    def __init__(self, prior=prior_func, corr=None, config_file=None,
                 run_id="test", obsname=None, burstname=None, gtiname=None,
                 continuous=True, maxgap=2,
                 interp='linear', smooth=0.02, model = settle,
                 theta= (0.58, 0.013, 0.4, 3.5, 1.0, 1.0, 1.5, 11.8),
                 sampler='emcee', fluen=True, alpha=False, pflux=True,
                 numburstssim=3, bc=2.21, ref_ind=1, threads = 4,
                 test_model=True, restart=False, **kwargs):
        """
        Initialise a Beans object. In addition to the parameters below,
        the user can specify additional parameter relevant to the various
        samplers: nwalkers, nsteps, stretch_a for emcee; and nlive, dlogz
        for dynesty. If these are not specified default values will be
        applied (when :meth:`Beans.do_run` is called).

        :param prior: prior function to use
        :param corr: correction function for bursts, or None
        :param config_file: file to read in configuration parameters (in
          which case the keyword params below are ignored)
        :param run_id: string identifier for the run, used to label all the
          result files, and where you want output to be saved
        :param obsname: name of the file including the flux history, from which
          the mdot is estimated to generate to generate the burst train
          (set ``obsname=None`` for a non-contiguous, or "ensemble" mode run)
        :param burstname: name of the burst data file, listing the bursts
        :param gtiname: name of the GTI file, set to ``None`` to turn off
          checking
        :param continuous: for burst "train" modes, set to True to generate a
          continuous train (i.e. using generate_burst_train rather than
          the newer punkt_train)
        :param interp: interpolation mode for the flux; possible values are
          'linear', or 'spline'
        :param smooth: smoothing factor for spline interpolation
        :param model: burst model to use, one of settle or grid_interp
        :param theta: initial centroid values for walker model parameters, with
          *X*, *Z*, *Q_b*, *d*, *xi_b*, *xi_p*, and (optionally) *mass*,
          *radius*, & *f_t*
        :param fluen: set to True (default) to include the fluences in the
          data for comparison, or False to omit
        :param alpha: set to True to include the alphas in the
          data for comparison, or False (default) to omit
        :param pflux: set to True (default) to include the peak fluxes in
          the data for comparison, if they have been supplied
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
        :param threads: number of threads for emcee to use (e.g. number of
          cores your computer has). Set to ``None`` to use all available
        :param test_model: flag to test the model during the setup process
        :param restart: set to ``True`` to continue a previously interrupted run

        :result: Beans object including all the required data
        """

        # Some housekeeping

        if 'ndim' in kwargs.keys():
            logger.warning('parameter ndim is redundant (ignored), setting from len of param array')
        if 'numburstsobs' in kwargs.keys():
            logger.warning('parameter numburstsobs is redundant (ignored), setting from len of burst data')
        if 'gti_checking' in kwargs.keys():
            logger.warning('parameter gti_checking is redundant (ignored), setting from value of gtiname param')

        # Conversion factor between model predicted burst energy and
        # observed fluence. Multiply by this value to convert the burst
        # energy from settle (in 1e39, the observer frame) to fluence at 1
        # kpc in units of 1e-6 erg/cm^2

        self.dist_fac = (1./(4*np.pi*u.kpc**2)).decompose()
        self.fluen_fac = ((1e39*u.erg * self.dist_fac)
            / (1e-6*u.erg/u.cm**2)).decompose()

        # Eddington luminosity and error, already multiplied by the
	# distance factor & divided by flux units (i.e. units of 1e-9
	# erg/cm^2/s kpc^2); for use comparing the peak fluxes in runmodel

        self.L_Edd = (L_EDD*self.dist_fac/FLUX_U).decompose().value
        self.L_Edd_err = (L_EDD_ERR*self.dist_fac/FLUX_U).decompose().value

        # Conversion factor between persistent flux (in units of 1e-9
        # erg/cm^2/s) and the accretion rate

        self.r1 = (1e-9*u.erg/u.cm**2/u.s*u.kpc**2
            / (u.km*const.c)**2).decompose()

        # Conversion factor for redshift, GM/Rc^2

        self.gmrc2 = (2.*const.G*const.M_sun/(const.c**2*u.km)).decompose()

        print("# ---------------------------------------------------------------------------#")
        logger.info("initialising Beans object ...")

        # MCU: this is a solution to load packaged data compatible
        # with python 3.6 - 3.10
        # however for python 3.10 and later this might become deprecated
        # alternative method is "from importlib.resources import files"
        # see https://setuptools.pypa.io/en/latest/userguide/datafiles.html

        self.data_path = os.path.join(os.path.dirname(__file__), 'data')
        # print("data_path = " + data_path)

        # here we want to parse the model function, which should have a
        # string representation something like <function settle at 0x1599b68b0>
        _tmp = str(model).split(' ')
        model_type, self.model_name = _tmp[0], None
        if len(_tmp) > 1:
            self.model_name = _tmp[1]
        if (model_type != '<function') | (self.model_name not in['settle','grid_interp']):
            logger.error('model unknown or invalid, use one of settle or grid_interp')
            return
        elif self.model_name == 'grid_interp':
            if not has_multiepoch_mcmc:
                logger.error('can\'t find multiepoch_mcmc, cannot initialise GridInterpolator')
                return
            # initialise the GridInterpolator object
            self.model_file = '/'.join((self.data_path, 'burst_model_grid.txt'))
            self.gi = grid_interpolator.GridInterpolator(file=self.model_file,
                reconstruct=True)
            self.grid_mdot_min = min(self.gi.grid['mdot'])
            self.grid_mdot_max = max(self.gi.grid['mdot'])
        else:
            self.model_file, self.gi = None, None

        # Only want to set the default values if both obsname and burstname
        # are not set (indicating a default run). This because setting
        # obsname=None is also how we indicate an "ensemble" mode run
        # and we might want to set burstname=None if we're doing some
        # simulations, e.g. using the :meth:`Beans.sim_data` method).
        if (obsname is None) & (burstname is None):
            obsname = os.path.join(self.data_path, '1808_obs.txt')
            burstname = os.path.join(self.data_path, '1808_bursts.txt')

        if run_id is None:
            # run_id = os.path.join(self.data_path, '1808/test1')
            # moving away from including the path in the run_id
            run_id = 'test1'

        # apply the keyword values or defaults

        self.run_id = run_id
        self.version = __version__
        self.obsname = obsname
        self.burstname = burstname
        self.gtiname = gtiname
        self.ref_ind = ref_ind
        self.bc = bc
        self.continuous = continuous
        # the original "train" model generates a continuous train of
        # bursts covering the entire history
        self.train_model = generate_burst_train
        if (not continuous) & (maxgap > 0):
            # the newer version generates a train but limits the
            self.train_model = punkt_train
            self.maxgap = maxgap

        self.theta = theta

        self.numburstssim = numburstssim
        self.lnprior = prior
        self.corr = corr
        self.model = model
        self.threads = threads
        self.sampler = sampler

        # Here we want to be agnostic about the sampler, so as to avoid
        # having to code up every possible attribute for the __init__
        # method. Instead we check the kwargs for these parameters,
        # including:
        #
        # :param nwalkers: number of walkers for the emcee run
        # :param nsteps: number of MCMC steps to run
        # :param stretch_a: the Goodman & Weare (2010) stretch move scale parameter, passed to emcee

        for key in kwargs:
            if key in ['stretch_a','nwalkers','nsteps', # emcee params
                'nlive','dlogz', # dynesty params
                'outdir']: # bilby param
                setattr(self, key, kwargs[key])
            else:
                logger.warning('ignored init option {}={}'.format(key, kwargs[key]))

        # following parameters control what goes into the likelihood
        # calcualtion

        self.cmpr_fluen = fluen
        self.cmpr_alpha = alpha
        self.cmpr_pflux = pflux

        # check for the config_file, and if detected read in (will
        # override the defaults above):

        config_file_exists = None
        if (config_file is not None):
            if (os.path.exists(config_file)):
                config_file_exists = True
                cmpr_alpha, cmpr_fluen, cmpr_pflux = False, True, True
                logger.info('Reading run params from {} ...'.format(config_file))
                self.read_config(config_file)
                logger.info('... done')
                # special here for the alpha and fluen parameters, which
		# are replaced by the actual data values (if the option is
		# True); (pre-v2.10 config files don't list alpha or fluen)
                # pflux is even more confusing, because this is the name
                # for the persistent fluxes (while the peak burst fluxes
                # are bpflux)
                if hasattr(self, 'alpha'):
                    self.cmpr_alpha = self.alpha
                    alpha = self.alpha
                if hasattr(self, 'fluen'):
                    self.cmpr_fluen = self.fluen
                    fluen = self.fluen
                if hasattr(self, 'pflux'):
                    self.cmpr_pflux = self.pflux
                    pflux = self.pflux

            else:
                logger.error('config file not found, applying keywords')
                config_file_exists = False

        # Some checks here

        if (sampler != 'emcee') & (not self.HAS_BILBY):
            logger.error('using samplers other than emcee requires bilby')
            return

        if alpha and (not fluen):
            logger.error('need to include fluences as well as alphas')
            return

        if alpha and (model == 'grid_interp'):
            logger.error('grid_interp doesn\'t currently provide alphas, set alpha=False')
            return

        if self.lnprior(self.theta) == -np.inf:
            logger.error('supplied parameter vector is excluded by the prior')
            return

        # for some older runs we need to add parameters here, that
        # are not listed in the .ini file

        if (self.sampler == 'emcee') & (not hasattr(self, 'stretch_a')):
            logger.warning('stretch_a not set or provided by .ini file, assuming default')
            self.stretch_a = 2.0

        # below set the parameters which are not part of the config
        # file

        self.restart = restart

        # number of dimensions for the parameter array
        # we want to record whether we're including systematic errors
        # here; previously this was via the (boolean) has_systematic

        self.ndim = len(self.theta)
        if (self.ndim < 6) | (self.ndim > 9):
            logger.error('number of dimensions of input parameter vector should be 6-9')
            return
        self.num_systematic = self.ndim-8
        if (self.num_systematic > 0) & (not LNPROB_USES_LNLIKE_SYS):
            logger.warning('''likelihood calculation does not currently support systematic errors;
                you need to swap out the lnlike method for lnlike_sys''')
        if ((self.ndim == 9) & (not self.cmpr_fluen)) | \
            ((self.ndim == 10) & (not (self.cmpr_alpha and self.cmpr_fluen))):
            logger.warning("systematic errors are provided for ignored quantities!")

        self.gti_checking = self.gtiname is not None

        # determines whether will run as a train of bursts or non-contiguous
        # bursts ("ensemble" mode); previously numerical, default is 1 (True),
        # which means a burst train will be generated; set obsname=None for
        # non-contiguous (ensemble mode) run

        self.train = (self.obsname is not None)

        self.bstart_err = BSTART_ERR.to('d').value
        self.M_NS = M_NS
        self.R_NS = R_NS

        self.matches = None

        # Read in all the measurements and set up all the parameters
        # This function now operates on the Beans object directly, and the
        # required attributes:
        # x, y, yerr, tref, bstart, pflux, pfluxe, tobs, fluen, fluen_err,
        #     st, et
        # are set in that routine
        # bypasses the earlier init function, and instead calls get_obs
        # directly

        get_obs(self, logger, alpha, fluen, pflux)

        # pre-calculate the sigmas and other parameters, for use in lnlike

        self.inv_sigma2 = 1. / self.yerr[:self.ly] ** 2
        # In versons 2.56.1 and earlier, there was a missing factor of
        # 2\pi in this expression!
        self.log_2pi_sigma2 = np.log(2.*np.pi/self.inv_sigma2)
        if self.train:
            # self.log_inv_sigma2[self.ref_ind] = 0.  # no contribution from the reference burst
            self.log_2pi_sigma2[self.ref_ind] = 0.  # no contribution from the reference burst

        # Set interpolation mode, and define averaging function

        if not hasattr(self, 'interp'):
            self.interp = interp
        assert self.interp in ('linear','spline')
        if self.interp == 'linear':
            self.a, self.b = get_a_b(self.pflux, self.pfluxe, self.tobs)
            self.mean_flux = mean_flux_linear
            self.smooth = None
        else:
            if not hasattr(self, 'smooth'):
                self.smooth = smooth
            self.tck_s = splrep(self.tobs, self.pflux, s=self.smooth)
            self.mean_flux = mean_flux_spline

        # some final checks for default sampler

        if (self.sampler == 'emcee'):
            if not hasattr(self, 'nwalkers'):
                self.nwalkers=100
            if not hasattr(self, 'nsteps'):
                self.nsteps=100

        print("# ---------------------------------------------------------------------------#")

        # # -------------------------------------------------------------------------#
        # # TEST THE MODEL WORKS
        # # -------------------------------------------------------------------------#
        if test_model:

            logger.info("testing the model works & is valid ...")


            test, valid, test2 = runmodel(self.theta, self, debug=False) # set debug to True for testing
            logger.info("\nresult: {} ({})".format( test, "OK" if valid else
                "can't match observations & model predictions"))

            # MCU Note: commented out - no interactive windows for automated testing
            # self.plot_model(test2)

            if (not valid) & (self.numburstsobs > 0):
                logger.warning ('''
the model is not valid. You need to adjust the model
                parameters to better suit the data.  ''')


    def __str__(self):
        """
        Show the parameters that the code has been intialised with
        For restart runs could include the number of steps that has
        already been done
        """

        return """== beans dataset =============================================================
See https://beans-7.readthedocs.io

Run ID: {}
{}{}Burst data file: {}
  comprising {} observed bursts{}
No. of bursts to simulate: {} ({} mode{})
  model: {}
  sampler: {}{}
  likelihood includes: {}times{}{}{}
  {}/{} threads
Initial parameters:
{}
==============================================================================""".format(self.run_id,
            'Observation data file: {}\n  bolometric correction: {}\n'.format( self.obsname, self.bc) if self.train else '',
            '' if self.gtiname is None else 'GTI data file: {}\n'.format(self.gtiname),
            self.burstname, self.numburstsobs,
            '' if (self.obsname is None) | ~self.continuous else ', ref. to #{}'.format(self.ref_ind),
            self.train+self.numburstssim*(1+self.train) if self.continuous else self.numburstsobs,
            'train' if self.train else 'ensemble',
            ', continuous' if (self.continuous & self.train) else '',
            self.model_name,
            # sampler-specific part
            self.sampler, (' with {} walkers, {} steps{}, a={}'.format(
                self.nwalkers, self.nsteps, ', resuming' if self.restart else '', self.stretch_a)
                if (self.sampler=='emcee') | (self.sampler=='bilby') else
                ' with nlive={}, dlogz={}'.format(self.nlive, self.dlogz)),
            'recurrence ' if not self.train else '',
            ', fluences' if self.cmpr_fluen else '',
            ', alphas' if self.cmpr_alpha else '',
            ', peak fluxes' if self.cmpr_pflux else '',
            self.threads, cpu_count(),
        self.theta_table(self.theta, indent=2) )


    def theta_table(self, theta, indent=0):
        """
        Format the run parameter vector as a table
        Could include the errors for a neatly formatted way to present
        results

        :param theta: the model parameter tuple, with *X*, *Z*, *Q_b*, *d*,
          *xi_b*, *xi_p*, *mass*, *radius*, and (optionally) *f_t*
        :param indent: number of characters to indent the string from the left
        """

        X, Z, Q_b, dist, xi_b, xi_p, *extra = theta
        mass, radius, f_t = extra+[self.M_NS, self.R_NS, 1.0][len(extra):]

        result = """#X = {:.3f} \\ hydrogen mass fraction
#Z = {:.5f} \\ CNO mass fraction
#Q_b = {:.3f} \\ base flux [MeV/nucleon]
#d = {:.2f} kpc \\ source distance
#xi_b = {:.3f} \\ anisotropy factor for burst emission
#xi_p = {:.3f} \\ anisotropy factor for persistent emission""".format(
    X, Z, Q_b, dist, xi_b, xi_p)
        result = result+"\n#M_NS = {:.3f} M_sun \\ neutron star mass".format(mass)
        if len(theta) <= 6:
            result += ' (fixed)'
        result = result+"\n#R_NS = {:.3f} km \\ neutron star radius".format(radius)
        if len(theta) <= 7:
            result += ' (fixed)'
        if self.num_systematic == 1:
            return (result+"""
#f_t = {:.3f} \\ systematic error term for burst times""".format(
    f_t)).replace('#',' '*indent)

        return result.replace('#', ' '*indent)


    def mdot_Edd(self, X, radius):
        """
        Calculate the Eddington accretion rate (per unit area) as used
        by Settle, for converting to physical units. The Eddington rate is
        (1.75*(1.7/(1+G.X))*(1e-8)*(5.01837638e24))/(G.R*G.R)
        (settle.cc, line 131) where G.X is the hydrogen mass fraction, and
        G.R is the radius in cm; the constant is  M_sun/(4*365.25×86400)
        (i.e. the conversion of  M_sun /yr to g/s) to better than 1 part
        in 1000) in the NS frame

        :param X: accreted H-fraction
        :param radius: NS radius (km)

        :return: accretion rate in g/cm^2/s
        """

        return (1.75e-8*1.7/(1+X)*5.01837638e24)/(radius*1e5)**2 * u.g/u.cm**2/u.s
        # alternative value here for the grid interpolation
        # Johnston et al. 2020: The Eddington-limited accretion rate,
        # return 8.775E4*u.g/u.cm**2/u.s


    def save_config(self, file=None, clobber=True):
        """
        Routine to write all the configuration parameters to a file, as a
        record of the run; but also to more easily replicate or continue
        a run. (Incomplete) list of all the parameters saved to the
        configuration file:
        ``run_id``, ``obsname``, ``burstname``, ``gtiname``, ``alpha``,
        ``fluen``, ``ref_ind``, ``bc``, ``interp``, ``smooth``, ``theta``,
        ``numburstssim``, ``prior``, ``nwalkers``, ``nsteps``, ``threads``

        :param file: name of file to save the config as. If None, the run_id
          will be used as a prefix
        :param clobber: set to True to overwrite any existing file
        """

        if file is None:
            file = self.run_id+'.ini'

        if (not clobber) and (os.path.isfile(file)):
            logger.error('config file already exists, skipping write')
            return
        else:
           cfgfile = open(file, "w")

           Config = ConfigParser()
           Config.add_section("beans")
           Config.set("beans", "run_id", self.run_id)
           Config.set("beans", "version", __version__)
           Config.set("beans", "model", str(self.model))

           Config.add_section("data")
           Config.set("data", "obsname", str(self.obsname))
           Config.set("data", "burstname", self.burstname)
           Config.set("data", "gtiname", str(self.gtiname))
           Config.set("data", "alpha", str(self.cmpr_alpha))
           Config.set("data", "fluen", str(self.cmpr_fluen))
           Config.set("data", "pflux", str(self.cmpr_pflux))

           if self.obsname is not None:
               # These parameters only used for "train" mode
               Config.set("data", "tref", str(self.tref))
               Config.set("data", "ref_ind", str(self.ref_ind))
               Config.set("data", "bc", str(self.bc))
               Config.set("data", "interp", str(self.interp))
               if self.interp == 'spline':
                   Config.set("data", "smooth", str(self.smooth))
               Config.set("beans", "continuous", str(self.continuous))
               if not self.continuous:
                   Config.set("beans", "maxgap", str(self.maxgap))

           # Config.add_section("emcee")
           Config.add_section("sampler")
           Config.set("sampler", "theta", str(self.theta))
           Config.set("sampler", "numburstssim", str(self.numburstssim))
           Config.set("sampler", "prior", str(self.lnprior))
           Config.set("sampler", "corr", str(self.corr))
           Config.set("sampler", "threads", str(self.threads))
           Config.set("sampler", "sampler", str(self.sampler))
           if (self.sampler == 'emcee') | (self.sampler == 'bilby'):
               Config.set("sampler", "nwalkers", str(self.nwalkers))
               Config.set("sampler", "nsteps", str(self.nsteps))
           elif (self.sampler == 'dynesty'):
               Config.set("sampler", "nlive", str(self.nlive))
               Config.set("sampler", "dlogz", str(self.dlogz))
           if (self.sampler == 'emcee'):
               Config.set("sampler", "stretch_a", str(self.stretch_a))
           if (self.sampler != 'emcee'):
               Config.set("sampler", "outdir", self.outdir)

           Config.write(cfgfile)
           cfgfile.close()


    def read_config(self, file=None):
        """
        Routine to read all the configuration parameters from a file, to
        more easily replicate or continue a run

        :param file: name of file to read the config from.
        """

        if file is None:
            run_id = os.path.join(self.data_path, 'beans.ini')

        int_params = ('ref_ind', 'numburstssim', 'nwalkers', 'nsteps', 'threads', 'maxgap', 'nlive')
        float_params = ('bc', 'smooth', 'tref', 'stretch_a', 'dlogz')

        if not os.path.isfile(file):
            logger.error('config file not found')
            return

        config = ConfigParser(allow_no_value=True)
        config.read(file)

        # Loop over sections, attributes

        for section in config.sections():
            # print("Section: %s" % section)
            if section == 'emcee':
                # read the older formats prior to the different sampling options
                setattr(self, 'sampler', 'emcee')
            for option in config.options(section):
                # print(
                #     "x %s:::%s:::%s"
                #     % (option, config.get(section, option), str(type(option))))
                if option == 'theta':
                    setattr(self, option,
                        tuple(map(float, config.get(section, option)[1:-1].split(', '))))
                elif option in float_params:
                    setattr(self, option, config.getfloat(section, option))
                elif (option == 'prior'):
                    function_name = config.get(section, option).split(' ')[1]
                    if (option == 'prior') & (function_name != str(self.lnprior).split(' ')[1]):
                        logger.warning ('''config file lists the prior function as {},
                but supplied prior is {}
              To fully replicate the previous run you need to specify the
                same prior using the prior=beans.{} flag on init
'''.format(function_name, str(self.lnprior).split(' ')[1], function_name))
                elif (option == 'model'):
                    function_name = config.get(section, option).split(' ')[1]
                    if (function_name != str(self.model).split(' ')[1]):
                        logger.warning('''config file lists the model as {},
                but supplied value is {}
              To fully replicate the previous run you need to specify the
                same model using the model={} specification on init
'''.format(function_name, str(self.model).split(' ')[1], function_name))
                elif (option == 'corr'):
                    function_name = config.get(section, option)
                    if function_name != 'None':
                        function_name = function_name.split(' ')[1]
                    _scorr = str(self.corr)
                    if _scorr != 'None':
                        _scorr = _scorr.split(' ')[1]
                    if (function_name != _scorr):
                        logger.warning('''config file lists {} for the correction,
                but supplied value is {}
              To fully replicate the previous run you need to specify the
                same option/condition using the corr=beans.{} flag on init
'''.format(function_name, _scorr, function_name))
                elif option in int_params:
                    setattr(self, option, config.getint(section, option))
                elif (option == 'continuous') & (config.get(section, option) == 'False'):
                    self.continuous = False
                    self.train_model = punkt_train
                else:
                    # string options (including "None")
                    _value = config.get(section, option)
                    if _value == 'None':
                        setattr(self, option, None)
                    elif _value == 'True':
                        setattr(self, option, True)
                    elif _value == 'False':
                        setattr(self, option, False)
                    else:
                        setattr(self, option, _value)

        if (self.sampler != 'emcee') & (not hasattr(self, 'outdir')):
            setattr(self, 'outdir', BILBY_OUTPUT)


    def flux_to_mdot(self, X, dist, xi_p, mass, radius, flux=None):
        """
        Function to convert fluxes to accretion rate, in units of the
        Eddington rate, calculated using the mdot_Edd function

        This routine uses the precise calculation of Q_grav, i.e.  c^2z/(1+z)

        :param X: accreted H-fraction
        :param dist: source distance (kpc)
        :param xi_p: anisotropy of persistent emission
        :param mass: NS mass (M_sun)
        :param radius: NS radius (km)
        :param flux: flux value or array to convert. If None, then we just
          return the observed persistent flux and error

        :return: mdot, mdot_err OR individual flux values
        """

        if flux is None:
            return self.flux_to_mdot(X, dist, xi_p, mass, radius, self.pflux), \
                self.flux_to_mdot(X, dist, xi_p, mass, radius, self.pfluxe)

        # Johnston et al. 2020: The Eddington-limited accretion rate,
        # scaled independent of composition or NS radius.
        # for use with the grid_interp model

        _mdot_Edd = 8.775E4*u.g/u.cm**2/u.s
        if self.model_name == 'settle':
            # different prescription for settle, via the mdot_Edd function
            _mdot_Edd = self.mdot_Edd(X, radius)

        opz = 1./(np.sqrt(1.-self.gmrc2*mass/radius))

        return (self.r1*flux*self.bc*dist**2*xi_p*opz**2
            / (radius**2*(opz-1)) / _mdot_Edd ).decompose().value


    def lnlike(self, theta_in, x, y, yerr, components=False):
        """
        Calculate the "model" likelihood for the current walker position
        Calls runmodel which actually runs the model, either generating a
        burst train, or a set of runs for "ensemble" mode. Then extracts the
        relevant model outputs and calculates the likelihood.

        :param theta_in: parameter vector, with *X*, *Z*, *Q_b*, *d*, *xi_b*,
          *xi_p* and (optionally) *mass* & *radius*,
        :param x: the "independent" variable, passed to lnlike
        :param y: the "dependent" variable (i.e. measurements), passed to lnlike
        :param yerr: erorr estimates on y

        :return: likelihood, model result array
        """

        # define theta_in = model parameters, which we define priors for

        X, Z, Q_b, dist, xi_b, xi_p, *extra = theta_in
        mass, radius = extra + [self.M_NS, self.R_NS][len(extra):]

        # call model (function runmodel, in run_model.py) to generate the burst
        # train, or the set of bursts (for "ensemble" mode. In earlier versions
        # the corresponding IDL function was defined as
        # modeldata(base, z, x, r1, r2 ,r3)

        assert np.allclose(y, self.y)
        model, valid, model2 = runmodel(theta_in, self, debug=False)
        if not valid:
            return -np.inf, model

        # Final likelihood expression
        # Because the y (observed value) vector may or may not include the
        # alpha values, we need to truncate the other vectors here

        # cpts = (self.y - (model)) ** 2 * inv_sigma2 - (np.log(inv_sigma2))
        # cpts = (self.y - model[:self.ly]) ** 2 * self.inv_sigma2 - self.log_inv_sigma2
        cpts = -0.5 * (self.y - model[:self.ly]) ** 2 * self.inv_sigma2 + self.log_2pi_sigma2

        # if the peak fluxes are defined, add the components for the
        # likelihood here. First the lower limits for the non-PRE bursts:

        for _i in self.non_pre:
            _log_arg = .5-.5*erf((self.bpflux[_i]-model2['F_Edd'][0])/model2['F_Edd'][1]) 
            if _log_arg == 0.0:
                return -np.inf, model
            cpts = np.append(cpts, np.log(_log_arg) )

        # and then for the PRE bursts

        for _i in self.pre:
            cpts = np.append(cpts, 
          -0.5 * (((self.bpflux[_i]-model2['F_Edd'][0])/model2['F_Edd'][1])**2
          + np.log(2.*np.pi*model2['F_Edd'][1]**2) ) )

        # if the components flag is set, also add those to the model dict

        if components:
            model2['cpts'] = cpts

        # Test if the result string is defined here. It is, so we return the selected elements of result
        # instead of the downselection in model

        # Now also return the model
        return np.sum(cpts), model2


    def lnlike_sys(self, theta_in, x, y, yerr, components=False):
        """
        As for :meth:`Beans.lnlike`, but this version also allows the "systematic" error
        factor for the burst times. If you want to run with that
        you need to swap this function out for :meth:`Beans.lnlike` in :meth:`Beans.lnprob`

        :param theta_in: parameter vector, with *X*, *Z*, *Q_b*, *d*, *xi_b*,
          *xi_p*, and (optionally) *mass*, *radius* & *f_t*
        :param x: the "independent" variable, passed to lnlike
        :param y: the "dependent" variable (i.e. measurements), passed to lnlike
        :param yerr: erorr estimates on y

        :return: likelihood, model result array
        """

        # define theta_in = model parameters, which we define priors for

        X, Z, Q_b, dist, xi_b, xi_p, *extra = theta_in
        mass, radius, f_t = extra + [self.M_NS, self.R_NS, 1.0][len(extra):]

        # call model (function runmodel, in run_model.py) to generate the burst
        # train, or the set of bursts (for "ensemble" mode. In earlier versions
        # the corresponding IDL function was defined as
        # modeldata(base, z, x, r1, r2 ,r3)

        assert np.allclose(y, self.y)
        model, valid, model2 = runmodel(theta_in, self)
        if not valid:
            return -np.inf, model

        # To simplify final likelihood expression we define inv_sigma2 for each
	# data parameter that describe the error.  The variance for the
	# burst times is (we hypothesize) underestimated by some
	# fractional amount, f_t; in early versions we also had systematic
	# contributions to fluence and alpha, but these are no longer used

        # Forgot to include an extra err_fac for the extra time, now that
        # we keep the reference burst time in the data array (now fixed)
        # so the length of the err_fac times cpt is self.numburstsobs, not
        # self.numburstsobs-ato

        ato = int(self.train) # array "train" offset
        err_fac = np.concatenate(( np.full(self.numburstsobs, 1./f_t**2),
            np.full(self.numburstsobs, 1.0), np.full(self.numburstsobs-ato, 1.0)))
        inv_sigma2 = self.inv_sigma2[:self.ly] * err_fac[:self.ly]

        # Final likelihood expression
        # Because the y (observed value) vector may or may not include the
        # alpha values, we need to truncate the other vectors here
        # In versons 2.56.1 and earlier, there was a missing factor of
        # 2\pi in this expression!

        # cpts = (self.y - (model)) ** 2 * inv_sigma2 - (np.log(inv_sigma2))
        cpts = -0.5 * (self.y - model[:self.ly]) ** 2 * inv_sigma2 \
            + np.log(2.*np.pi/inv_sigma2)

        # if the peak fluxes are defined, add the components for the
        # likelihood here; also want to trap underflows
        # First the lower limits for the non-PRE bursts:

        for _i in self.non_pre:
            _log_arg = .5-.5*erf((self.bpflux[_i]-model2['F_Edd'][0])/model2['F_Edd'][1]) 
            if _log_arg == 0.0:
                return -np.inf, model

            cpts = np.append(cpts, np.log(_log_arg) )

        # and then for the PRE bursts

        for _i in self.pre:
            cpts = np.append(cpts, 
          -0.5 * (((self.bpflux[_i]-model2['F_Edd'][0])/model2['F_Edd'][1])**2
          + np.log(2.*np.pi*model2['F_Edd'][1]**2) ) )

        # if the components flag is set, also add those to the model dict

        if components:
            model2['cpts'] = cpts

        # Test if the result string is defined here. It is, so we return the selected elements of result
        # instead of the downselection in model

        # Now also return the model
        return np.sum(cpts), model2

    # -------------------------------------------------------------------------#
    # Finally we combine the likelihood and prior into the overall lnprob function, called by emcee

    # define lnprob, which is the full log-probability function
    def lnprob(self, theta_in, x, y, yerr):
        """
        The full log-probability function incorporating the priors (via
        the ``Beans.lnprior`` attribute), and model likelihood (via
        :meth:`Beans.lnlike`), that is passed to ``runemcee`` when creating
        the sampler (in the :meth:`Beans.do_run` method).

        :param theta_in: parameter vector, with *X*, *Z*, *Q_b*, *d*, *xi_b*,
          *xi_p*, *mass*, *radius*, and (optionally) *f_t*
        :param x: the "independent" variable, passed to :meth:`Beans.lnlike`
        :param y: the "dependent" variable (i.e. measurements), passed to
          :meth:`Beans.lnlike`
        :param yerr: erorr estimates on y

        :return: total (prior+model) likelihood, prior likelihood, model array
          (from :meth:`Beans.lnlike`)
        """

        lp = self.lnprior(theta_in)
        # Check if the parameters are consistent with the prior, and skip
        # the model run it if not
        if (not np.isfinite(lp)):
            return -np.inf, -np.inf, None

        # Now also returns the model, to accumulate along with the likelihoods
        # if you swap lnlike for lnlike_sys below (or vice versa) make
        # sure you update LNPROB_USES_LNLIKE_SYS at the start of this file

        like, model = self.lnlike_sys(theta_in, x, y, yerr)

        if (not np.isfinite(like)):
            return -np.inf, -np.inf, model

        # encoding below is so we have a suitable object for including in
        # the blobs (see the dtype specification in runemcee)

        return lp + like, lp, model_str(model).encode('ASCII')


    def plot(self, show_model=True, model=None, mdot=True, imatch=None,
        title=None, savefig=False):
        """
        Display a plot of the data and model results, for a burst train
        calculated with :func:`burstrain.generate_burst_train`;
        or burst "ensemble" data calculated with
        :func:`burstrain.burstensemble`. Adapted from the
        example at
        https://matplotlib.org/gallery/api/two_scales.html

        :param show_model: set to False to skip the model generation step, in which case only the observed data will be plotted
        :param model: array of packed model prediction, OR dict giving full model results
        :param mdot: flag to show mdot rather than flux (only possible if you're passing the full model)
        :param imatch: can pass a "matching" array for train models, to use instead of the automatically-calculated version
        :param title: add a title, if required
        :param savefig: set to True to save the figure to PDF
        """

        flux_colour = 'tab:red'
        bursts_colour = 'tab:blue'
        obs_colour = 'tab:grey'

        # for the default linear interpolation connect the flux
        # measurements by lines

        ls = '-'
        if self.interp == 'spline':
            ls = ''

        X, Z, Q_b, dist, xi_b, xi_p, *extra = self.theta
        mass, radius, f_t = extra+[self.M_NS, self.R_NS, 1.0][len(extra):]

        itoff = 1 # time offset for models produced with generate_burst_train
        if (model is None) & show_model:
            # no need to do the matching here
            test, valid, model = runmodel(self.theta, self,
                match=(False or not self.continuous),
                debug=False)

            show_model = valid

        elif model is not None:
            # assume a passed model is valid
            valid = show_model

        if (imatch is None) & ('imatch' not in model) & (valid & self.train & (self.numburstsobs > 0)):
            # ... but it's useful to know if it's possible
            # (not relevant for punkt_train, or ensemble mode)

            imatch = burst_time_match(self.ref_ind, self.bstart,
                model['iref'], np.array(model['time']))

            if imatch is None:
                logger.warning ("can't match predicted bursts to observations")

        elif 'imatch' in model:
            imatch = model['imatch']
            itoff = 0

        full_model = False  # Flag to remember whether we're plotting the
                            # full model output of generate burst train or
                            # the packed output array
        if type(model) == dict:
            full_model = True
            timepred = model["time"]
            if (len(timepred) == 1) | (not valid):
                logger.error('no/insufficient predicted times to show')
                show_model = False
            else:
                ebpred = np.array(model["fluen"])
        elif type(model) == list:
	    # The other way to return the model is as an array with the
	    # burst times, fluences and alphas all together. So unpack
	    # those here
            # This mode is deprecated hence the breakpoint
            breakpoint()
            timepred = model[:self.numburstssim+1]
            ebpred = np.array(model[self.numburstssim+int(self.train):self.numburstssim*2+int(self.train)]) * self.fluen_fac/xi_b/dist**2

        # Now set up the figure

        if self.train & show_model & (imatch is not None):
            fig, axs = plt.subplot_mosaic([['main'],['main'],['resid']])
            ax1 = axs['main']
        else:
            fig, ax1 = plt.subplots()
        # fig.figure(figsize=(10,7))

        ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
        if self.train & self.gti_checking:
            #plot satellite gtis
            for i in range(1,len(self.st)):
                ax1.axvspan(self.st[i],self.et[i],facecolor='0.5', alpha=0.2)
            ax1.axvspan(self.st[0],self.et[0], facecolor='0.5',alpha=0.2,label='Satellite GTIs')

        # first show the flux/mdot
        if mdot and full_model:
            # this is the usual (only?) plot style, showing mdot vs. time,
            # along with the burst model and observation comparison
            ax1.set_ylabel('Accretion rate (fraction of Eddington)', color=flux_colour)
            _mdot, _mdot_err = self.flux_to_mdot(X, dist, xi_p, mass, radius)
            if self.train:
                # show the mdot values
                ax1.errorbar(self.tobs, _mdot, _mdot_err,
                    marker='.', ls=ls, color=flux_colour, label='mdot')
                if self.interp == 'spline':
                    t = np.arange(min(self.tobs), max(self.tobs), 0.1)
                    ax1.plot(t, self.flux_to_mdot(X, dist, xi_p, mass, radius,
                        BSpline(*self.tck_s)(t)), color=flux_colour)
                # show the time of the "reference" burst
                # ax2.axvline(timepred[self.iref], c='k')
                if self.numburstsobs > 0:
                    ax2.axvline(self.bstart[self.ref_ind], c='k', ls='--')
                ax1.set_xlabel("Time (days after MJD {})".format(self.tref))
            else:
                # show the ensemble comparison, which is much simpler
                ax1.errorbar(self.bstart, _mdot, _mdot_err, fmt='.',
                         color=flux_colour, label='mdot')
                ax1.set_xlabel("Epoch (MJD)")
        else:
            # showing here the persistent flux rather than mdot
            ax1.set_ylabel('Persistent flux ($10^{-9}\\,{\\rm erg\\,cm^{-2}\\,s^{-1}}$)', color=flux_colour)
            if self.train:
                ax1.errorbar(self.tobs, self.pflux, self.pfluxe,
                    marker = '.', ls=ls, color=flux_colour, label = 'persistent flux')
                if self.interp == 'spline':
                    t = np.arange(min(self.tobs), max(self.tobs), 0.1)
                    ax1.plot(t, BSpline(*self.tck_s)(t), color=flux_colour)
                if show_model:
                    ax2.scatter(timepred[0], ebpred[0], marker = '*',color=bursts_colour,s = 100)
                ax1.set_xlabel("Time (days after MJD {})".format(self.tref))
            else:
                # "ensemble" mode plot vs. epoch, rather than observation time
                ax1.errorbar(self.bstart, self.pflux, self.pfluxe, fmt='.', color=flux_colour,
                         label='persistent flux')
                ax1.set_xlabel("Epoch (MJD)")

        ax1.tick_params(axis='y', labelcolor=flux_colour)

        # Plot the bursts here
        ax2.set_ylabel("Fluence ($10^{-6}\\,{\\rm erg\\,cm^{-2}}$)", color=obs_colour)
        if self.train:
            if self.bstart is not None:
                # Plot the observed bursts, if available
                # ax2.scatter(tobs,ebobs, color = 'darkgrey', marker = '.', label='observed', s =200)
                ax2.errorbar(self.bstart[self.ifluen], self.fluen[self.ifluen],
                    yerr=self.fluene[self.ifluen],
                    color=obs_colour, linestyle='', marker='.', ms=13,
                    label='observed bursts')
                for i in range(self.numburstsobs):
                    if (i not in self.ifluen) & (i != self.ref_ind):
                        ax2.axvline(self.bstart[i], color=obs_colour, ls='--')

            if show_model:
                ax2.scatter(timepred[itoff:], ebpred, marker = '*',color=bursts_colour,s = 100, label = 'predicted bursts')
            # we have time but not fluence for the first burst
                if itoff == 1:
                    ax2.axvline(timepred[0], color=bursts_colour, ls='--')
                # and the averaged mdot over the burst interval (predicted)
                av_mdot = []
                for i in range(len(timepred)-1):
                    av_mdot.append(self.flux_to_mdot(X, dist, xi_p, mass, radius,
                        self.mean_flux(timepred[i], timepred[i+1], self)))
                av_mdot.insert(0, av_mdot[0])
                ax1.step(timepred, av_mdot, where='pre', color=flux_colour)

                if imatch is not None:
                    # show the burst time comparison
                    resid = -(self.bstart-np.array(timepred)[imatch])*24.
                    axs['resid'].plot(imatch, resid,
                        linestyle='', marker='.', ms=13, color=obs_colour)
                    for i in range(self.numburstsobs):
                        axs['resid'].annotate(' {}'.format(i+1),
                            (imatch[i], resid[i]) )
                    axs['resid'].axhline(0.0, color=obs_colour, ls='--')
                    axs['resid'].set_xlabel('Burst number (predicted)')
                    axs['resid'].set_ylabel('Time offset (hr)')
                    logger.info ('RMS obs-model offset = {:.4f} hr'.format(
                        np.sqrt(np.mean(resid**2))))

        else:
            # ensemble mode plot here
            ax2.scatter(self.bstart, self.fluen, color = 'darkgrey', marker = '.', label='observed', s =200)
            if show_model:
                ax2.scatter(self.bstart, ebpred, marker = '*',color=bursts_colour,s = 100, label = 'predicted')

        ax2.tick_params(axis='y', labelcolor=obs_colour)

        if title is not None:
            plt.title(title)

        fig.tight_layout()  # otherwise the right y-label is slightly clipped

        fig.legend(loc=1)

        if savefig:
            file = f'{self.run_id}_plot.pdf'
            logger.info ('Saving figure plot to {}'.format(file))
            fig.savefig(file)

        fig.show()

    def sim_data(self, file=None, ensemble=False):
        """
        This method generates a table of simulated data based on the
        current theta vector, that can be used for tests of the MCMC
        method

        Errors are taken as the mean of the observational errors; bursts
        matching observed bursts (according to the matching algorithm)
        are flagged as ``True`` in the 6th column

        In lieu of an output file, you can just paste the results into a
        text file and use as the input (specify as ``burstname``)

        :param file: file to save the results to
        :param ensemble: set to True to output "ensemble" data, i.e.  including persistent fluxes, and with recurrence times rather than burst times

        :return:
        """

        logger.info ("generating simulated data based on parameter vector\n theta={}...".format(self.theta))

        test, valid, full_model = runmodel(self.theta, self,
            match=(self.numburstsobs > 0))
        _n = len(full_model['time'])

        if _n <= 1:
            logger.error('no bursts generated by model, adjust your model parameters')
            return

        if file is None:
            f = sys.stdout
            _comm = ''
        else:
            f = open(file, 'w')
            print ('''# File {}
# Created {}
#
# Simulated burst data generated by beansp v{} (see
#   https://github.com/adellej/beans for the code, and
#   https://beans-7.readthedocs.io/en/latest/index.html for documentation)
#
# Parameter vector theta={}'''.format(file, time.strftime("%Y-%m-%d %H:%M:%S"),
            __version__, self.theta), file=f)
            _comm = '# '

        X, Z, Q_b, dist, xi_b, xi_p, *extra = self.theta
        mass, radius, f_t = extra+[self.M_NS, self.R_NS, 1.0][len(extra):]

        opz = 1./(np.sqrt(1.-self.gmrc2*mass/radius))
        print ('{}corresponding to source at\n{}  d={:.2f} kpc, xi_b={:.4f}, xi_p={:.4f}, X={:.2f}, Z={:.3f}, M_NS={:.2f}, R_NS={:.2f}, 1+z={:.3f}, Q_b={:.2f}'.format(
            _comm, _comm, dist,  xi_b, xi_p, X, Z, mass, radius, opz, Q_b), file=f)
        if self.obsname is not None:
            print ('{}Persistent fluxes from {}, interp={}{}, bc={}'.format(
                _comm, self.obsname, self.interp,
                ('' if self.smooth is None else ', smooth={}'.format(self.smooth)), self.bc), file=f)

        if self.corr is not None:
            print ('{}Correction function {}'.format(_comm, self.corr), file=f)

        # assign fluence and alpha errors here, and match string, if
        # we have observed bursts to match to
        if self.numburstsobs > 0:
            _fluene, _alphae = np.mean(self.fluene), np.mean(self.alphae[1:])
            _matched = ['	'+str(i in full_model['imatch']) for i in range(_n)]
            if hasattr(self, 'tdele'):
                # only for ensemble mode input data, I think
                _tdele = mean(self.tdele)
            else:
                _tdele = 0.1 # hr
        else:
            _fluene, _alphae, _tdele = 0.01, 1.0, 0.1 # dummy values
            _matched = [''] * _n

        if self.train & (not ensemble):
            print ('{}\n{}"burst start time [d]"	"bolometric fluence [1e-6 erg/cm^2]"	"error on bolometric fluence"	"alpha"	"error on alpha"'.format(_comm, _comm), file=f)
            _line_fmt = '{:.5f}	{:.4f}	{:.4f}	{:.1f}	{:.1f}{}'
            for i in range(_n):
                if i == 0:
                    print(_line_fmt.format(full_model['time'][0]+self.tref,
                        full_model['fluen'][0],_fluene,0.,0., _matched[i]), file=f)
                else:
                    print(_line_fmt.format(full_model['time'][i]+self.tref,
                        full_model['fluen'][i-1], _fluene,
                        full_model['alpha'][i-1], _alphae, _matched[i]), file=f)

        else:
            # originally this was for the ensemble mode data, but I don't
            # think I ever tested it (and it looks like the code was just
            # a placeholder that replicated the input data).
            # The version implemented below is where we have the "train"
            # data but we're going to save it as ensemble mode, taking
            # each burst interval as a component of the "ensemble"

            print ('{}\n{}"burst epoch [MJD]"	"bolometric fluence [1e-6 erg/cm^2]"	"error on bolometric fluence"	"alpha value"	"error on alpha value"	"bolometric persistant flux [1e-9 erg/cm^2/s^-1]"	"error on bol.  per flux"	"inferred recurrence time [hr]"	"error on recurrence time"'.format(_comm, _comm), file=f)
            _line_fmt = '{:.5f}	{:.3f}	{:.3f}	{:.1f}	{:.1f}	{:.3f}	{:.3f}	{:.3f}	{:.3f}'
            # for i in range(len(self.bstart)):
            for i in range(1, _n):
                print(_line_fmt.format(
                    full_model['time'][i], full_model['fluen'][i-1], _fluene,
                    full_model['alpha'][i-1], _alphae,
                    self.mean_flux(full_model['time'][i-1],full_model['time'][i],self), 0.0,
                    (full_model['time'][i]-full_model['time'][i-1])*24., _tdele
                    # self.bstart[i], self.fluen[i], self.fluene[i],
                    # self.alpha[i],self.alphae[i],self.pflux[i],
                    # self.pfluxe[i],self.tdel[i],self.tdele[i]
                    ), file=f)

        if file is not None:
            f.close()
            logger.info('Simulated data saved in file {}.'.format(file))

    # -------------------------------------------------------------- #

    def distance_limit(self, distrange=[1,10], Xmin= 0.01, gridsize=10, numtrains=100, skiplo=True):
        """
        Generate burst trains and test how many of them are consistent with the GTIs
        """

        assert len(distrange) == 2

        X_0, Z, Q_b, dist, xi_b, xi_p, *extra = self.theta
        mass, radius, f_t = extra+[self.M_NS, self.R_NS, 1.0][len(extra):]

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
                    test, valid, result = runmodel(theta_1, self.y, 0.0, self.bstart,
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
                            test, valid, result = runmodel(theta_1, self.y, 0.0, self.bstart,
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

    def do_run(self, sampler=None, plot=False, analyse=True, burnin=2000, **kwargs):
        """
        This routine runs the chain for as many steps as is specified in
        the init call.  Emcee parameters are defined in runemcee module.
        Have previously used multiprocessing to speed things up, but that's
        not currently active

        :param sampler: define the sampler to use; default is emcee
        :param plot: set to True to do the plot of the initial guess. This
          seems redundant since it's also plotted at the __init__ stage
        :param analyse: set to True to call do_analysis automatically once the
          chains finish
        :param burnin: number of steps to ignore at the start of the run,
          passed to do_analysis

        :return:
        """

        # sampler attribute passed to do_run overrides the existing value
        if hasattr(self, 'sampler'):
            if sampler is not None:
                if sampler != self.sampler:
                    logger.warning ('sampler={} passed to do_run, overriding Beans attribute {}'.format(
                        sampler, self.sampler))
            else:
                sampler = self.sampler
        else:
            sampler = 'emcee'

        if sampler not in ('emcee','bilby','dynesty'):
            logger.error('sampler {} not yet implemented'.format(sampler))

        # store the most recent value of the sampler
        self.sampler = sampler
        if self.sampler == 'emcee':

            # Check here if we've already run the analysis
            if hasattr(self, 'reader'):
                logger.warning ('''re-running the sampler after calling do_analysis can lead
                    to memory issues; proceed with caution!''')
                value = input('              Press [RETURN] to continue: ')

            # Want to avoid overwriting existing log & config files

            if (self.restart is False) and (os.path.exists(self.run_id+'.h5')):
                logger.error ('run will overwrite existing log file {}, set restart=True to extend'.format(self.run_id+'.h5'))
                return

        else:
            if not hasattr(self, 'outdir'):
                self.outdir = BILBY_OUT

        if (self.sampler == 'emcee') | (self.sampler == 'bilby'):

            # check for start parameters and set defaults

            if not hasattr(self, 'nwalkers'):
                self.nwalkers = 200
            if not hasattr(self, 'nsteps'):
                self.nsteps = 100
            if not hasattr(self, 'stretch_a'):
                self.stretch_a = 2.0

        elif (self.sampler == 'dynesty'):

            # check for start parameters and set defaults

            if not hasattr(self, 'nlive'):
                self.nlive = 1000
            if not hasattr(self, 'dlogz'):
                self.nlive = 1.
            # params are passed via kwargs, so need to make sure they're
            # set here
            kwargs['nlive'] = self.nlive
            kwargs['dlogz'] = self.dlogz

        # Check the number of threads here

        ncpu = cpu_count()
        if self.threads > ncpu:
            logger.error ('number of threads is greater than number of CPUs ({} > {}); this seems ill-advised'.format(self.threads, ncpu))
            return

        if (os.path.exists(self.run_id+'.ini')):
            logger.warning ('run will overwrite existing config file {}'.format(self.run_id+'.ini'))
            value = input('              enter Y[RETURN] to continue: ')
            if (value != 'y') and (value != 'Y'):
                logger.info ('do_run terminated')
                return

        self.save_config(clobber=True)

        print("# ---------------------------------------------------------------------------#")
        print (self)
        print("# ---------------------------------------------------------------------------#")
        # Testing the various functions. Each of these will display the likelihood value, followed by the model-results "blob"
        logger.info("testing the prior and likelihood functions..")
        if LNPROB_USES_LNLIKE_SYS:
            print("lnlike_sys:", self.lnlike_sys(self.theta, None, self.y, self.yerr))
        else:
            print("lnlike:", self.lnlike(self.theta, None, self.y, self.yerr))
        if self.sampler == 'emcee':
            print("lnprior:", self.lnprior(self.theta))
            print("lnprob:", self.lnprob(self.theta, None, self.y, self.yerr))
        print("# ---------------------------------------------------------------------------#")

        # for supplied positions, you want to check that they are all
        # valid
        # I think this is only necessary to define for emcee and bilby/emcee

        if 'pos' in kwargs:

            # now also has the option to read in the positions from a
            # pickle file

            if type(kwargs['pos']) == str:
                logger.info ('reading in positions from file {}...'.format(kwargs['pos']))
                new_pos = pickle.load(open(kwargs['pos'],'rb'))
                logger.info ('... done')
            else:
                new_pos = kwargs['pos'].copy()

            if np.shape(new_pos) != (self.nwalkers, self.ndim):
                logger.error ('supplied positions array has wrong dimensions; {} != {}'.format(np.shape(new_pos), (self.nwalkers, self.ndim)))
                return None

            logger.warning('walkers will start at provided position vector.\n              checking supplied positions...')
            # assert len(kwargs['pos']) == self.nwalkers
            # we need to make sure the supplied positions will give a
            # valid model BUT also that the likelihoods will not all be
            # -inf
            bad = np.full(len(new_pos), False)
            for i, pos in enumerate(new_pos):
                # _test, _valid, _model = runmodel(pos, self)
                # print (pos, _test, _valid, _model)
                _test, _prior, _model = self.lnprob(pos, None, self.y, self.yerr )
                # print (pos, _test, _prior, _model)
                bad[i] = not np.isfinite(_test)
            nbad = len(np.where(bad)[0])

            if nbad == 0:
                logger.info ('... done. all OK, continuing')
            # elif nbad == self.nwalkers:
            elif nbad > self.nwalkers*2./3.:
                # if you have too many bad positions, you'll get an error
                # running with too many duplicates; this is a rule of
                # thumb which probably needs testing
                logger.error ('... done. Too many (>2/3) of the supplied positions are invalid (lnprob = -np.inf), run will likely fail')
                return
            else:
                logger.info ('... done. {}/{} positions invalid; redistributing...'.format(nbad, self.nwalkers))
                for i in (np.where(bad)[0]):
                    new_pos[i] = new_pos[np.random.choice(np.where(~bad)[0])] #+ scale*np.random.randn(ndim)
                logger.info ('... done')

            kwargs['pos'] = new_pos

        if plot:
            logger.info("plotting the initial guess.. (you want the predicted bursts to match approximately the observed bursts here)")
            # make plot of observed burst comparison with predicted bursts:
            self.plot(title='Initial guess of parameters')
            value = input('Press [RETURN] to continue: ')

        _start = time.time()

        if sampler == 'emcee':
            # This is the original MCMC implementation for beansp, which uses vanilla emcee
            #   https://github.com/dfm/emcee
            # run the chains and save the output as a h5 file
            # TODO to simplify the subsequent analysis I think this object should be added to the Beans object
            result = runemcee(self.nwalkers, self.nsteps,
                self.theta, self.lnprob, self.lnprior, None, self.y, self.yerr,
                self.run_id, self.restart, self.threads, self.stretch_a, **kwargs)

        elif sampler == 'bilby':
            # This is (will eventually?) be the native bilby sampler, or perhaps bilby's implementation of emcee,
            # for comparison; see https://bilby-dev.github.io/bilby/api/bilby.core.sampler.emcee.Emcee.html#bilby.core.sampler.emcee.Emcee
            logger.error ("'bilby' sampler not yet implemented; might call bilby/emcee, or bilby's native sampler")

        elif sampler == 'dynesty':
            # bilby's implementation of dynesty, see
            #
            result = runbilby(self, sampler=sampler, **kwargs)

        _end = time.time()
        logger.info ("  run_id {} took {:.1f} seconds for {} steps".format(
            self.run_id, _end-_start, self.nsteps))

        if analyse & (sampler == 'emcee'):
            logger.info ("running analses with do_analysis...")
            self.do_analysis(burnin=burnin)
        #if self.restart == False:
            #sampler.reset()

# -------------------------------------------------------------------------#

    def plot_autocorr(self, samples=None, reader=None, title=None, savefile=None, figsize=(8,5) ):
        """
        This method shows the estimated autocorrelation time for the run
        Extracted from do_analysis, and originally based on the analysis
        described in the example for emcee:
        https://emcee.readthedocs.io/en/stable/tutorials/autocorr

        :param samples: numpy array with samples, to calculate the
          autocorrelation from
        :param reader: object to get the chains from (if not supplied via
	  the samples parameter), via the get_chain method. TODO just pass
          the samples to make this more general
        :param title: title for the figure; set to None for some generic info, or False to skip
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


        if (samples is None) and (reader is None):
            # load in sampler:
            logger.info ('loading reader from {}.h5...'.format(self.run_id))
            reader = emcee.backends.HDFBackend(filename=self.run_id+".h5")

        if reader is not None:
            samples = reader.get_chain(flat=False)
            #tau = 20
            tau = reader.get_autocorr_time(tol=0) #using tol=0 means we'll always get an estimate even if it isn't trustworthy.
            thin = int(0.5 * np.min(tau)) # seems to be ignored - dkg
            print(f"The autocorrelation time for each parameter as calculated by emcee is: {tau}")
            print ("  mean {:.1f}, min {:.1f}, max {:.1f}".format(np.mean(tau),
                min(tau), max(tau)))

        # alternate method of checking if the chains are converged:
        # This code is from https://dfm.io/posts/autocorr/

        # get autocorrelation time:

        # loop through 10 parameters and plot the evolution of the
        # autocorrelation time estimate for each

        f = plt.figure(figsize=figsize)

        param = ["$X$", "$Z$", "$Q_{\mathrm{b}}$", "$d$", "$\\xi_b$",
            "$\\xi_p$", "$M$", "$R$", "$f_{\mathrm{E}}$", "$f_{\mathrm{a}}$"]

        for j in range(self.ndim):
            chain = samples[:, :, j].T
            # print(np.shape(sampler))

            N = np.exp(np.linspace(np.log(100), np.log(chain.shape[1]), self.ndim)).astype(int)
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

        if title is not False:
            plt.title('beansp v{} run {}'.format(
                self.version, self.run_id) if title is None else title,loc='right')

        plt.plot(N, np.array(N)/50.0, "--k")# label=r"$\tau = N/50$")
        plt.legend(fontsize='large',loc='best',ncol=2) #bbox_to_anchor=(0.99, 1.02)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        if savefile is not None:
            logger.info ('saving autocorrelation plot to {}'.format(savefile))
            plt.savefig(savefile)
        else:
            logger.info ('skipping autocorrelation plot save')

        plt.show()


    def burst_table(self, show=True, predicted=False, key=None):
        """
	This method creates and, optionally, displays a table from the
        model and observed bursts.

        The table is returned as an astropy Table object that can be saved
        as an MRT file, like so: ``tab.write('test.dat',format='mrt')``
	although you will have to then edit the file afterwards to
        complete the metadata

        :param show: set to True to display each row (in LaTeX format)
        :param predicted: set to True to show predicted bursts also (not yet implemented)
        :param key: for the case where we have more than one set of predictions (e.g. for different numbers of predicted bursts), the key will identify which to return

        :return: astropy table of burst properties
        """

        if not self.HAS_CONCORD:
            logger.error ('alpha calculation requires concord')
            return None

        if not hasattr(self, 'model_pred'):
            logger.error ('no model predictions available, run the comparison first')
            return None

        if (key is None) & (len(self.model_pred['time_stats'].keys()) > 1):
            logger.error ('multiple solutions are available, please specify one with the key keyword')
            print (self.model_pred['time_stats'].keys())
            return None

        elif key is None:
            key = list(self.model_pred['time_stats'].keys())[0]

        timepred = [x[0] for x in self.model_pred['time_stats'][key]]
        ref_tpred = np.argmin(np.abs(self.bstart[self.ref_ind]-timepred))

        # have to do different matching procedures if not continuous
        # cf. with the comparison option of do_analysis

        if self.continuous:
            imatch = burst_time_match(self.ref_ind, self.bstart, ref_tpred, np.array(timepred))
        else:
            imatch = burst_time_match(0, self.bstart, 0, np.array(timepred))

        bursts = Table()
        bursts['num'] = np.arange(self.numburstsobs)+1
        bursts['num'].info.description = 'Burst number'
        bursts['minbar_id'] = np.full(self.numburstsobs, 9999) # filler
        bursts['minbar_id'].info.description = 'MINBAR DR1 ID'
        bursts['time'] = Time(self.bstart+self.tref, format='mjd')
        # we can set this, but it doesn't seem to be written when output
        # as MRT
        bursts['time'].info.description = 'Burst start time'
        bursts['bfluen'] = MaskedColumn(self.fluen, mask=self.fluen <= 0.,
            unit=FLUEN_U, description='Integrated burst fluence')
        bursts['e_bfluen'] = MaskedColumn(self.fluene, mask=self.fluen <= 0.,
            unit=FLUEN_U, description='Error on burst fluence')
        # following list comprehension mask expression is to trap the
        # blanks that come from MINBAR
        # bursts['alpha_obs'] = MaskedColumn(self.alpha, mask=self.alpha <=0.,
        bursts['alpha_obs'] = MaskedColumn(self.alpha,
            mask=[float(x) <= 0.0 if x != '--' else True for x in self.alpha],
            description='Burst alpha-value')
        # bursts['e_alpha_obs'] = MaskedColumn(self.alphae, mask=self.alpha <=0.,
        bursts['e_alpha_obs'] = MaskedColumn(self.alphae,
            mask=[float(x) <= 0.0 if x != '--' else True for x in self.alpha],
            description='Error on alpha-value')

        if show:
            # show the headers for the LaTeX table
            print ('''
\\begin{tabular}{ccccccc} 
  \hline
        & MINBAR & Start  & $\Delta t$ & Fluence & \multicolumn{2}{c}{$\\alpha$-value} \\\\
  Burst & ID     & (MJD)  & (hr)       & $10^{-6}\\ {\\rm erg\,cm^{-2}}$ & Measured & Inferred \\\\
  \hline ''')

        # Now loop over the bursts and calculate the derived quantities
        dt, e_dt, alpha, e_alpha, E_alpha = [], [], [], [], []
        _sel = np.where(np.array(self.model_pred['partition']) == key)[0]
        for i in np.arange(self.numburstsobs):
            if imatch[i] > 0:
                if imatch[i]-imatch[i-1] == 1:
                    # no missed bursts
                    _dt = (self.bstart[i]-self.bstart[i-1])*24.
                    perflx = self.mean_flux(self.bstart[i-1], self.bstart[i], self)
                else:
		    # one or more missed bursts, so use the model
		    # predictions for the previous time
                    # NOTE have to preserve the key selection here
                    _dt = [(self.bstart[i] - self.model_pred['times'][x][imatch[i]-1])*24.0 for x in _sel]
                    perflx = [self.mean_flux(self.model_pred['times'][x][imatch[i]-1], self.bstart[i], self) for x in _sel]

                # calculate the alpha based on the model-predicted time, but only
                # if the measured fluence is nonzero
                if self.fluen[i] > 0.0:
                    _alpha = self.cd.alpha(_dt,(self.fluen[i], self.fluene[i]), perflx).value
                else:
                    _alpha = (0.0, 0.0, 0.0)

                if np.shape(_dt) == ():
                    # scalar _dt, exact recurrence time
                    dt.append(_dt)
                    e_dt.append(0.)
                    if show:
                        print ("{} & [minbar ID] & {} & {:.2f} & {} & {} & {} \\\\".format(
                            bursts['num'][i], bursts['time'][i], dt[-1],
                            (strmeas(bursts['bfluen'][i], bursts['e_bfluen'][i]) if bursts['bfluen'][i] > 0. else '--'),
                            strmeas(bursts['alpha_obs'][i], bursts['e_alpha_obs'][i]),
                            (strmeas(_alpha[0], _alpha[1], _alpha[2]) if bursts['bfluen'][i] > 0. else '--')))
                else:
                    dt_stats = np.percentile(np.array(_dt), [16,50,84])
                    dt.append(dt_stats[1])
                    e_dt.append((dt_stats[2]-dt_stats[0])*0.5)
                    if show:
                        print ("{} & [minbar ID] & {} & {} & {} & {} & {} \\\\".format(
                            bursts['num'][i], bursts['time'][i],
                            strmeas(dt[-1], e_dt[-1]),
                            (strmeas(bursts['bfluen'][i], bursts['e_bfluen'][i]) if bursts['bfluen'][i] > 0. else '--'),
                            strmeas(bursts['alpha_obs'][i], bursts['e_alpha_obs'][i]),
                            (strmeas(_alpha[0], _alpha[1], _alpha[2]) if bursts['bfluen'][i] > 0. else '--')))
                alpha.append(_alpha[0])
                e_alpha.append(_alpha[1])
                E_alpha.append(_alpha[2])
            else:
                # not sure what this case is for (excluded bursts?)
                # I think first observed burst, where we have no prior event
                dt.append(0.)
                e_dt.append(0.)
                alpha.append(0.)
                e_alpha.append(0.)
                E_alpha.append(0.)
                if show:
                    print ("{} & [minbar ID] & {} & \\nodata & {} & \\nodata & \\nodata \\\\".format(
                        bursts['num'][i], bursts['time'][i],
                        strmeas(bursts['bfluen'][i], bursts['e_bfluen'][i])))

        if show:
            # print the end of the LaTeX table
            print ('  \hline\n\end{tabular}')

        bursts['trec'] = MaskedColumn(dt, mask=np.array(dt) <=0., unit=u.hr,
            description='Burst recurrence time')
        bursts['e_trec'] = MaskedColumn(e_dt, mask=np.array(dt) <=0., unit=u.hr,
            description='Error on recurrence time')
        bursts['alpha_mod'] = MaskedColumn(alpha, mask=np.array(alpha) <= 0.,
            description='Model-informed alpha value')
        bursts['e_alpha_mod'] = MaskedColumn(e_alpha, mask=np.array(alpha) <= 0.,
            description='Lower limit on model alpha')
        bursts['E_alpha_mod'] = MaskedColumn(E_alpha, mask=np.array(alpha) <= 0.,
            description='Upper limit on model alpha')

        # set the formats here. These should be OK for most cases!
        bursts['trec'].info.format='5.2f'
        bursts['e_trec'].info.format='5.2f'
        bursts['alpha_mod'].info.format='6.2f'
        bursts['e_alpha_mod'].info.format='4.2f'
        bursts['E_alpha_mod'].info.format='4.2f'

        return bursts


    def write_param_uncert(self, clobber=False, latex=True):
        """
        Calculate and save the parameter uncertainties; previously part of
        do_analysis

        :param clobber: if True, overwrite file without prompting
        :param latex: if True, also display the information in LaTeX form, suitable for inclusion in a paper
        """

        if not hasattr(self, 'reader'):
            logger.error ('no samples available, run do_analysis first')
            return

        file = f'{self.run_id}_parameterconstraints_pred.txt'

        latex_header = "Parameter"
        latex_rows = ["$X$", "$Z$", "$Q_{\\rm b}$", "$d$", "$\\xi_b$",
            "$\\xi_p$", "$M_{\\rm NS}$", "$R_{\\rm NS}$", "$g$", "$1+z$", "$f_t$"]
        latex_unit = ["", "", "MeV/nucleon", "kpc", "", "", "$M_\\odot$",
            "km", "$10^{14}\\ {\\rm cm\\,s^{-2}}$", "", ""]

        if (clobber is False) and (os.path.exists(file)):
            logger.error ('will overwrite existing parameter file {}, set clobber=True to replace'.format(file))
            return

        # constants:
        c = const.c.to('cm s-1')
        G = const.G.to('cm3 g-1 s-2')

        header='''beansp v{} parameter file

run_id {}, nsteps_completed={}, skipping {} steps for burnin

Each row has the 50th percentile value, upper & lower 68% uncertainties

Rows are:
H mass fraction (X), metallicity (Z), base flux (Q_b), distance (d, kpc),
persistent anisotropy factor (xi_p), burst anisotropy factor (xi_b)
'''.format(__version__, self.run_id, self.nsteps_completed,self.samples_burnin)

        # Now get the remainder of the parameter uncertainties and
        # save to the text file

        n_samples = np.shape(self.samples)[0]
        if hasattr(self, 'model_pred'):
            parts = list(set(self.model_pred['partition']))
            sel = np.array(self.model_pred['partition']) == parts[0]
            if len(parts) > 1:
                # something missing here?
                pass
        else:
            # select all the samples for the parameter ranges
            parts = ['all']
            sel = np.full(n_samples, True)

        # have to duplicate the calculation of redshift and gravity here, as
        # those quantities are not saved to the samples
        if self.ndim >= 8:
            M = np.array(self.samples[:,6])*const.M_sun.to('g') #cgs
            R = np.array(self.samples[:,7])*1e5*u.cm #cgs

            redshift = np.power((1 - (2*G*M/(R*c**2))), -0.5).value
            gravity = (M*redshift*G/R**2 / (u.cm/u.s**2)).value #cgs

        with open(file, 'w') as f:

            for i, _part in enumerate(parts):
                Xpred = get_param_uncert(self.samples[sel,0])
                Zpred = get_param_uncert(self.samples[sel,1])
                basepred = get_param_uncert(self.samples[sel,2])
                dpred = get_param_uncert(self.samples[sel,3])
                # cosipred = get_param_uncert(cosi)
                xibpred = get_param_uncert(self.samples[sel,4])
                xippred = get_param_uncert(self.samples[sel,5])

                if len(parts) > 1:
                    latex_header += ' & {} ({:.2f}%)'.format(_part,
                        100*self.model_pred['part_stats'][_part]/n_samples)
                else:
                    latex_header += ' & Value'
                latex_rows[0] += ' & {}'.format(strmeas(Xpred[0], Xpred[2], Xpred[1]))
                latex_rows[1] += ' & {}'.format(strmeas(Zpred[0], Zpred[2], Zpred[1]))
                latex_rows[2] += ' & {}'.format(strmeas(basepred[0], basepred[2], basepred[1]))
                latex_rows[3] += ' & {}'.format(strmeas(dpred[0], dpred[2], dpred[1]))
                latex_rows[4] += ' & {}'.format(strmeas(xibpred[0], xibpred[2], xibpred[1]))
                latex_rows[5] += ' & {}'.format(strmeas(xippred[0], xippred[2], xippred[1]))

                if self.ndim >= 7:
                    masspred = get_param_uncert(self.samples[sel,6])
                    latex_rows[6] += ' & {}'.format(strmeas(masspred[0], masspred[2], masspred[1]))
                if self.ndim >= 8:
                    radiuspred = get_param_uncert(self.samples[sel,7])
                    redshiftpred = get_param_uncert(redshift[sel])
                    gravitypred = get_param_uncert(gravity[sel]/1e14)
                    latex_rows[7] += ' & {}'.format(strmeas(radiuspred[0], radiuspred[2], radiuspred[1]))
                    latex_rows[8] += ' & {}'.format(strmeas(gravitypred[0], gravitypred[2], gravitypred[1]))
                    latex_rows[9] += ' & {}'.format(strmeas(redshiftpred[0], redshiftpred[2], redshiftpred[1]))

                if self.ndim >=9:
                    ftpred = get_param_uncert(self.samples[sel,8])
                    latex_rows[10] += ' & {}'.format(strmeas(ftpred[0], ftpred[2], ftpred[1]))

                if len(parts) > 1:
                    header += '''
Sample subset {} of {}, label {}, {}%'''.format(i+1,len(parts),_part,
    100*self.model_pred['part_stats'][_part]/n_samples)

	        # save to text file with columns: value, upper uncertainty,
	        # lower uncertainty
                # have to have a few options here depending on whether the
                # mass is included or not
                # TODO also include f_E, f_a if they're included

                if self.ndim == 6:
                    np.savetxt(f, (Xpred, Zpred, basepred, dpred, xippred, xibpred),
                        fmt='%9.6f', header=header)
                elif self.ndim == 7:
                    header = header+',\nNS mass (M_sun)'
                    np.savetxt(f, (Xpred, Zpred, basepred, dpred, xippred, xibpred,
                        masspred), fmt='%9.6f', header=header)
                elif self.ndim == 8:
                    header = header+',\nNS mass (M_sun), NS radius (km), gravity (g, 10^14 cm/s^2), redshift (1+z)'
                    np.savetxt(f, (Xpred, Zpred, basepred, dpred, xippred, xibpred,
                        masspred, radiuspred, gravitypred, redshiftpred),
                        fmt='%9.6f', header=header)
                else:
                    header = header+',\nNS mass (M_sun), NS radius (km), gravity (g, 10^14 cm/s^2), redshift (1+z), f_t'
                    np.savetxt(f, (Xpred, Zpred, basepred, dpred, xippred, xibpred,
                        masspred, radiuspred, gravitypred, redshiftpred, ftpred),
                        fmt='%9.6f', header=header)

                if i < len(parts)-1:
                    sel = np.array(self.model_pred['partition']) == parts[i+1]
                    header = ''

            # output as LaTeX as well, to paste into a paper (for example)

            print ('\\begin{tabular}{c'+'c' * len(parts)+'l}\n\\hline')
            print (latex_header+' & Units \\\\\n\\hline')
            for i, row in enumerate(latex_rows):
                print (row+' & '+latex_unit[i]+' \\\\')
            print ('\\hline')


    def do_analysis(self, options=['autocor','posteriors'],
                          part=None, truths=None, burnin=2000,
                          ilab=None, title=None, savefig=False):
        """
        This method is for running standard analysis and displaying the
        results.
        Nothing is returned, but various options for analysis are
        available, which will (optionally) populate some of the
        :class:`beansp.Beans` attributes, including

        | reader - sampler object read in from HDF file
        | sampler - chain data
        | nsteps_completed - number of steps completed
        | samples - flattened samples, omitting first samples_burnin
        | last - last walker positions
        | probs - likelihoods (total, prior, and contributions from each parameter) for the last walker positions
        | cc - ChainConsumer object with samples and derived cosi, gravity, redshift
        | cc_parameters - plot labels for cc object
        | model_pred - dictionary with model realisations read in from the "blobs"

	By default the method will also create several files, labeled by
        the run_id; drawn from

        | {}_autocorrelationtimes.pdf (via plot_autocorr)
        | {}_posteriors.pdf
        | {}_parameterconstraints_pred.txt

        TODO: need to reorganise a bit, and add more options

        :param options: array of strings corresponding to various analysis
          options, listed in the analyses dict below
	:param part: string or array "partition" dividing the set of
          samples into two or more separate groups for analysis
        :param truths: parameter vector to overplot on (one of the) corner
          plots (TODO: need to check if >1 corner plot are selected
          +truths, which will likely result in an error due to
          incompatible number of parameters)
        :param burnin: number of steps to discard when plotting the
          posteriors. If -ve, discard all but that number
        :param title: set to a string to add a title to the plot, or False to omit; if set to None will print some generic information
        :param ilab: set to a list of (integer) inclination values to label in the anisotropy model plot ('fig8')
        :param savefig: set to True to save figures to .pdf files, False to skip

        :return: none
        """

        # constants:
        c = const.c.to('cm s-1')
        G = const.G.to('cm3 g-1 s-2')

        # list of available analyses

        analyses = {'autocor': 'autocorrelation times as a function of timestep',
                    'chain': 'visualisation of the walker evolution',
                    'posteriors': 'raw posteriors and the input values',
                    'mrcorner': 'corner plot with M, R, g and 1+z',
                    'fig6': 'corner plot with xi_b, xi_p, d, Q_b, Z',
                    'fig8': 'xi_b vs. xi_p and models for comparison',
                    'converge': 'check the convergence via posterior evolution',
                    'comparison': 'observed and predicted burst times, fluences',
                    'last': 'analyse last walker positions' }

        if options == 'all':
            options = analyses.keys()

        # check the chosen option is one of those implemented

        for option in options:
            if option not in analyses.keys():
                logger.error ('{} is not an available analysis option; choose from'.format(option))
                for key in analyses.keys():
                    print ('  {}: {}'.format(key, analyses[key]))
                return

        # blobs are not yet implemented for samplers other than emcee

        if ('comparison' in options) & (self.sampler != 'emcee'):
            logger.error ('blobs are not yet available for bilby runs')
            return

        # chains are not available for dynesty

        if ('chain' in options) & (self.sampler == 'dynesty'):
            logger.error ('chains are not available for dynesty runs')
            return

        if not hasattr(self, 'reader'):

            # if you haven't already done so, read in the walker data

            logger.info ("Reading in samples...")# to calculate autocorrelation time...")

            # load in sampler:
            if self.sampler == 'emcee':
                self.reader = emcee.backends.HDFBackend(filename=self.run_id+".h5")
                self.result = self.reader.get_chain(flat=False)
            else:
                # bilby outputs
                self.reader = bilby.result.read_in_result(outdir=self.outdir, label=self.run_id)
                if self.sampler == 'bilby':
                    self.result = self.reader.walkers.transpose(1,0,2)
                elif self.sampler == 'dynesty':
                    self.result = None
                    self.samples = self.reader.posterior.iloc[:,:self.ndim].to_numpy()
                else:
                    logger.error ('sampler {} not yet implemented'.format(self.sampler))
                    return

            # temporary for bilby output, might work
            # import pickle
            # f = open('bilby_out/emcee_canon_b/sampler.pickle','rb')
            # self.reader = pickle.load(f)
            # f.close()

            # Read in the full chain to get the number of steps completed
            if self.result is not None:
                self.nsteps_completed = np.shape(self.result)[0]
                logger.info ("... done. Got {} steps completed".format(self.nsteps_completed))
            else:
                self.nsteps_completed = 0
                logger.info ("... done.")

            self.samples_burnin = None
            self.models_burnin = None

        # convert -ve burnin here

        if burnin < 0:
            if -burnin > self.nsteps_completed:
                logger.warning ('steps to retain {} is > total ({}), ignoring'.format(-burnin, self.nsteps_completed))
            burnin = max([0, self.nsteps_completed+burnin])

        # want to make sure we're using at least about 1000 samples for
        # our statistics

        # if burnin >= self.nsteps_completed*0.9:
        if (self.sampler == 'emcee') & ((self.nwalkers * (self.nsteps_completed - burnin) < 1000) | \
            (self.nsteps_completed <= burnin)):
            logger.warning('discarding burnin {} will leave too few steps ({} total), ignoring'.format(
                burnin, self.nsteps_completed))
            burnin = 0

        if (burnin != self.samples_burnin) | \
            ((burnin != self.models_burnin)&(self.models_burnin is not None)):

            # moved burnin to be a parameter, so we can pass that from do_run
            # and also keep track of what we've used, so that we don't
            # need to read in the samples and create the ChainConsumer
            # object every time

            # print ("Reading in flattened samples to show posteriors...")
            # samples = self.reader.get_chain(flat=True, discard=burnin)
            if self.result is not None:
                samples = self.result[burnin:,:,:]
                self.last = samples[-1,:,:]
                self.samples = samples.reshape((-1,self.ndim))
            self.samples_burnin = burnin

            cosi_2 = 1/(2*self.samples[:,5])
            cosi = 0.5/(2*(self.samples[:,5]/self.samples[:,4])-1)

            # here now create the overall chainconsumer object that will
            # enable the various posterior plots below

            labels = {"X": "$X$", "Z": "$Z$", "Qb": "$Q_b$ (MeV)",
                "d": "$d$ (kpc)", "xi_b": "$\\xi_b$", "xi_p": "$\\xi_p$"}

            if self.ndim >= 7:
                labels["M"] = "$M$ ($M_\odot$)"
                mass = self.samples[:,6]
                # masspred = get_param_uncert(mass)
            if self.ndim >= 8:
                labels["R"] = "$R$ (km)"
                radius = self.samples[:,7]
                # radiuspred = get_param_uncert(radius)

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

                # redshiftpred = get_param_uncert(redshift)
                # gravitypred = get_param_uncert(gravity)

            if self.ndim >= 9:
                labels["ft"] = "$f_t$"

            # now create the chainconsumer object

            _plot_labels = labels
            _plot_labels['cosi'] = '$\cos i$'
            if self.ndim >= 8:
                _plot_labels['g'] = '$g$ (cm s$^{-2}$)'
                _plot_labels['1+z'] = '$1+z$'
                _samples =np.column_stack((self.samples, cosi, gravity, redshift))
            else:
                _samples =np.column_stack((self.samples, cosi))

            self.cc = ChainConsumer()
            # self.cc.add_chain(_samples, parameters=_labels)
            self.cc.add_chain(_samples, parameters=list(_plot_labels.values()),
                name='beansp v{} run {} last {}/{}'.format(
                self.version, self.run_id, self.nsteps_completed-self.samples_burnin, self.nsteps_completed))
            # configure params below copied from Adelle's jupyter notebook
            # we apply them here for consistency across all the posterior
            # plots
	    # despite the original sigmas setting, only 2 contours are
	    # shown...? (the first two, 1 & 2 sigma)
            self.cc.configure(usetex=True, serif=True,
                flip=False, summary=False,
                bins=0.7, # has the effect of light smoothing of the histograms
                diagonal_tick_labels=False, max_ticks=3, shade=True, \
                shade_alpha=1.0 ,bar_shade=True, tick_font_size='xx-large', \
                label_font_size='xx-large',smooth=True, \
                sigma2d=False, sigmas=[1,2]) #np.linspace(0, 3, 4))
            self.samples = _samples # keep the samples up to date
            self.cc_parameters = _plot_labels
            self.cc_nchain = 1 # initially

        # ---------------------------------------------------------------------#
        if 'last' in options:

            # show likelihood contributions for last step
	    # TODO could also use the reader.get_last_sample method, at
	    # least for p_tot, but this won't give you the likelihood
            # contributions. It does include as the 4th element the prior
            # and model realisations (blobs) though

            logger.info('calculating probabilities for last walker step...')

            probs = pd.DataFrame(columns = ['p_tot','prior','p_time','p_fluen','p_alpha'],
                index = np.arange(self.nwalkers) )
            p_fluen, p_alpha = 0.0, 0.0
            for _i in np.arange(np.shape(self.last)[0]):
                ptot, model = self.lnlike_sys(self.last[_i,:], None, self.y, self.yerr, components=True)
                ato = int(self.train)
                if model is None:
                    # The model will not always be valid
                    p_time, p_fluen, p_alpha = 0., 0., 0.
                else:
                    p_time = -0.5*np.sum(model['cpts'][:self.numburstsobs-ato])
                    if self.cmpr_fluen:
                        p_fluen = -0.5*np.sum(model['cpts'][self.numburstsobs-ato:2*self.numburstsobs-ato])
                    if self.cmpr_alpha:
                        p_alpha = -0.5*np.sum(model['cpts'][2*self.numburstsobs-ato:])
                pprior = self.lnprior(self.last[_i,:])
                probs.loc[_i] = [ptot, pprior, p_time, p_fluen, p_alpha]

            self.probs = probs
            logger.info('... done. Likelihoods and breakdown stored in probs attribute')

        # ---------------------------------------------------------------------#
        # PLOTS
        # ---------------------------------------------------------------------#

        if 'autocor' in options:

            # plot autocorrelation times

            self.plot_autocorr(reader=self.reader, title=title,
                savefile='{}_autocorrelationtimes.pdf'.format(self.run_id)
                if savefig else None)

        #sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, self.lnprob, args=(self.x, self.y, self.yerr), backend=reader)

        # ---------------------------------------------------------------------#
        if 'chain' in options:

            # plot the chains:

            logger.info ("plotting the chains...")
            labels = ["$X$","$Z$","$Q_b$","$d$", "$\\xi_b$", "$\\xi_p$", "$M$", "$R$","$f_t$"]
            # plt.clf()
            fig, axes = plt.subplots(self.ndim, 1, sharex=True, figsize=(8, 9))

            for i in range(self.ndim):
                # Previously the transposed sampler object below meant
                # that we were plotting with walker number on the x-axis,
                # instead of step. Now fixed
                axes[i].plot(self.result[:,:,i], color="k", alpha=0.4)
                axes[i].yaxis.set_major_locator(MaxNLocator(5))
                axes[i].set_ylabel(labels[i])

            axes[self.ndim-1].set_xlabel("step number")
            if title is not False:
                axes[0].set_title('beansp v{} run {}'.format(
                    self.version, self.run_id) if title is None else title,loc='right')
            plt.tight_layout(h_pad=0.0)

            if savefig:
                savefile = '{}_chain-plot.pdf'.format(self.run_id)
                logger.info ('saving chain plot to {}'.format(savefile))
                plt.savefig(savefile)
            else:
                logger.info ('skipping chain plot save')

            plt.show()
            logger.info ("... done")

        # ---------------------------------------------------------------------#
        if ('posteriors' in options) | ('mrcorner' in options) | ('fig6'
in options):

            # common truths/title settings for all corner plots

            if truths is None:
                truths = list(self.theta)
            elif truths is False:
                truths = None

            if title is not None:
                if type(title) == str:
                    logger.warning("can't override title for posterior plots, sorry")
                    title=True
            else:
                title=True

        if 'posteriors' in options:

            # make plot of posterior distributions of your parameters:
            # (using the already-created ChainConsumer object)

            fig = self.cc.plotter.plot(parameters=list(self.cc_parameters.values())[:self.ndim],
                figsize="page", truth=truths, legend=title, display=False)

            fig.show()

            if savefig:
                # save the figure
                savefile = '{}_posteriors.pdf'.format(self.run_id)
                logger.info ('saving posteriors plot to {}'.format(savefile))
                plt.savefig(savefile)
            else:
                logger.info ('skipping posteriors plot save')

        # ---------------------------------------------------------------------#
        if ('mrcorner' in options) & (self.ndim < 8):
            logger.error ('can''t show M-R posteriors as one or both of M & R are fixed')
        elif 'mrcorner' in options:

            # mrgr = np.column_stack((mass, radius, gravity, redshift))

            # plot with chainconsumer:
            # cc = ChainConsumer()
            # cc.add_chain(mrgr, parameters=["M", "R", "g", "1+z"])

            fig = self.cc.plotter.plot(parameters=[self.cc_parameters[x] for x in ["M", "R", "g", "1+z"]],
                truth=truths, figsize="page", legend=title, display=False)

            fig.show()

            if savefig:
                # save the figure
                savefile = '{}_massradius.pdf'.format(self.run_id)
                logger.info ('saving mass-radius posteriors plot to {}'.format(savefile))
                plt.savefig(savefile)
            else:
                logger.info ('skipping mass-radius posteriors plot save')

        # ---------------------------------------------------------------------#
        if 'fig6' in options:

            # fig6data = np.column_stack((xip, xib, distance, base, Z, X))
            # fig6data = np.column_stack((X, Z, base, distance, xib, xip))

            # plot with chainconsumer:

            # cc = ChainConsumer()

            # cc.add_chain(fig6data, parameters=["X", "$Z$", "$Q_b$ (MeV)",
            #     "$d$ (kpc)", "$\\xi_b$", "$\\xi_p$"])\
            fig = self.cc.plotter.plot(
                parameters=[self.cc_parameters[x] for x in ['X','Z','Qb','d','xi_b','xi_p']],
                truth=truths, figsize="page", legend=title, display=False)

            fig.show()

            if savefig:
                # save the figure
                savefile = '{}_fig6.pdf'.format(self.run_id)
                logger.info ('saving restricted posteriors plot to {}'.format(savefile))
                plt.savefig(savefile)
            else:
                logger.info ('skipping restricted posteriors plot save')

        # ---------------------------------------------------------------------#
        if ('converge' in options):

            # do the summary plot, comparing the two halves of the burnin

            _cc = ChainConsumer()
            _n = int(np.shape(self.samples)[0]/2)
            _cc.add_chain(self.samples[:_n,:self.ndim],
                parameters=list(self.cc_parameters.keys())[:self.ndim],
                name='{}-{}'.format(self.samples_burnin,
                self.samples_burnin+int(_n/self.nwalkers)))
            _cc.add_chain(self.samples[_n:,:self.ndim],
                parameters=list(self.cc_parameters.keys())[:self.ndim],
                name='{}-{}'.format(self.samples_burnin+int(_n/self.nwalkers),
                self.nsteps_completed))

            fig.show()

            if savefig:
                savefile = '{}_converge.pdf'.format(self.run_id)
                logger.info ('saving convergence check plot to {}'.format(savefile))
                _cc.plotter.plot_summary(
                    filename=savefile)#,figsize="page")
            else:
                fig = _cc.plotter.plot_summary()#figsize="page")
                logger.info ('skipping convergence check plot save')

        # ---------------------------------------------------------------------#
        if ('comparison' in options) & ((self.models_burnin != burnin) | (part is not None)):

            # and finally read in the model realisations
            # This loop can take a LOOOOOONG time for long runs

            logger.info ("reading in and processing blobs, ignoring first {}...".format(burnin))
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
            logger.info ("... done")

            # here we split up the solutions based on the number of bursts
            # in the model train, and populate the model_pred attribute,
            # which we use also in the plotting part
            # This approach is not perfect, as some solutions with different
            # numbers of bursts have the same imatch, so are effectively
            # equivalent;
	    # TODO should combine the matching section in the plot segment
	    # with this analysis, along with checks to combine solutions
	    # with different burst numbers but the same imatch

            numbursts_pred = [len(x) for x in time]

            # special here for the base20 run, to 3500 steps at least;
            # more generally want to have a way of doing this on the fly
            # part = ['loX-36' if ((_n == 36) & (self.samples[i,0] < 0.3))
            #     else 'loX-37' if ((_n == 37) & (self.samples[i,0] < 0.3))
            #     else 'hiX' #if (self.samples[i,0] > 0.3)
            #     for i, _n in enumerate(numbursts_pred)]
            if (part is None): # & (len(set(numbursts_pred)) > 1):
                part = numbursts_pred

            # to get the parameter middle values and uncertainty use the
            # functions get_param_uncert_obs and get_param_uncert_part
            # e.g.

            times = get_param_uncert_part(time, partition=part)

            # Here we calculate the parameter uncertainties on the
            # predcted fluences, for comparison with the observations

            # this fails with inhomogeneous arrays; need to do something
            # a bit more complicated
            # fpred = (np.array(e_b).T*self.fluen_fac/np.array(xib)/np.array(distance)**2).T
            fpred = [list(np.array(_e_b)*self.fluen_fac/self.samples[_i,4]/self.samples[_i,3]**2) for _i, _e_b in enumerate(e_b)]
            ebs = get_param_uncert_part(fpred, partition=part)

            alphas = get_param_uncert_part(alpha, partition=part)

	    # store these parameters and flag it so we don't need to
	    # calculate them again (unless burnin changes)

            # part_stats = None
            # if part is not None:
            #     part_stats = {x: len(np.where(np.array(part) == x)[0]) for x in set(part)}
            part_stats = {x: len(np.where(np.array(part) == x)[0]) for x in times.keys()}
            self.model_pred = { 'mdot': mdot,
                'times': time, 'time_stats': times,
                'numbursts': numbursts_pred, 'partition': part,
                'part_stats': part_stats,
                'e_b': e_b, 'e_b_stats': ebs,
                'alpha': alpha, 'alpha_stats': alphas }
            self.models_burnin = burnin

            # Here also we modify the ChainConsumer object if we have
            # multiple models

            if (len(times) != self.cc_nchain) & (part is not None):
                self.cc = ChainConsumer()
                self.cc_nchain = 0
                for _n in set(part):
                    # need to check that there are sufficient samples here
                    _sel = np.array(part) == _n
                    _check = np.shape(self.samples[_sel])[0]
                    if _check > 1000:
                        self.cc.add_chain(self.samples[_sel],
                            # name='{} bursts'.format(_n),
                            name = _n if type(_n) == str else str(_n),
                            parameters=list(self.cc_parameters.values()))
                        self.cc_nchain += 1
                    else:
                      logger.info ('skipping walkers for n={}, too few samples ({})'.format(_n, _check))
                # make sure we ad the same configuration as for the single
                # chain object
                self.cc.configure(usetex=True, serif=True,
                    flip=False, summary=False,
                    bins=0.7, # has the effect of light smoothing of the histograms
                    diagonal_tick_labels=False, max_ticks=3, shade=True, \
                    shade_alpha=1.0 ,bar_shade=True, tick_font_size='xx-large', \
                    label_font_size='xx-large',smooth=True, \
                    sigma2d=False, sigmas=[1,2]) #np.linspace(0, 3, 4))
                logger.info ('updated chain object with {} model classes'.format(self.cc_nchain))

        # ---------------------------------------------------------------------#
        if 'fig8' in options:

            # here we read in data from the anisotropy models. There's
            # probably a better way to do this, via concord (if it's
            # available)

            counts, ybins, xbins, image = plt.hist2d(self.samples[:,5],
                self.samples[:,4], bins=50, norm=LogNorm(), cmap='OrRd')

            xi_p_model2 = np.arange(0, 2.5, 0.01)
            xi_b_model2 = np.empty(len(xi_p_model2))

            for i in range(0,250):

                xi_b_model2[i] = 1./((1./(2*xi_p_model2[i])) + 0.5)

            # overplot the various models

            plt.plot(xi_p_model2, xi_b_model2, color = 'black',ls='-', label = 'Fujimoto (1988)')

            if Beans.HAS_CONCORD:
                # setup dict with list of models, legend labels and linestyles
                he16_models = {'he16_a': ('He \& Keek (2016) model A', '--'),
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

                    if ilab is not None:
                        for _inc in ilab:
                            _i = np.argmin(abs(model_theta-_inc))
                            plt.annotate('{}'.format(_inc), (model_xip[_i], model_xib[_i]))

            else:
                logger.warning ('''install concord if you want to overplot model curves
              See https://github.com/outs1der/concord''')

            plt.xlabel(r'$\xi_{\mathrm{p}}$',fontsize='xx-large')
            plt.ylabel(r'$\xi_{\mathrm{b}}$',fontsize='xx-large')

            plt.legend(loc='best',fontsize='large')

            plt.axis([0.,2.1,0.,2.1])

            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)

            plt.show()

            if savefig:
                savefile = '{}_xipvsxib_models_contourlines.pdf'.format(self.run_id)
                logger.info ('saving xi_p vs xi_b plot to {}'.format(savefile))
                plt.savefig(savefile)
            else:
                logger.info ('skipping xi_p vs xi_b plot save')

        # ---------------------------------------------------------------------#
        if 'comparison' in options:

            # make plot of observed burst comparison with predicted bursts:

            plt.figure(figsize=(10,7))

            # copy colors from plot method for consistency

            bursts_colour = 'tab:blue'
            obs_colour = 'tab:grey'

            # plt.scatter(self.bstart, self.fluen, color = 'black', marker = '.', label='Observed', s =200)
            times = self.model_pred['time_stats']
            ebs = self.model_pred['e_b_stats']
            if self.train:
                # 2-panel plot like in plot

                if self.continuous:
                    itoff = 1  # time offset for models produced with generate_burst_train
                else:
                    itoff = 0  # time offset for models produced with punkt_train

                fig, axs = plt.subplot_mosaic([['main'],['main'],['resid']], sharex=True)
                ax1 = axs['main']

                ax1.errorbar(self.bstart[self.ifluen], self.fluen[self.ifluen],
                    yerr=self.fluene[self.ifluen],
                    color=obs_colour, linestyle='', marker='.', ms=13,
                    label='observed')
                # non-redundant option for missing fluence values
                for i in range(self.numburstsobs):
                    if (i not in self.ifluen) & (i != self.ref_ind):
                        ax1.axvline(self.bstart[i], color=obs_colour, ls='--')

                if self.continuous:
                    # only show reference burst for continuous runs
                    ax1.axvline(self.bstart[self.ref_ind], c='k', ls='--')

                if self.gti_checking:
                    #plot satellite gtis
                    for i in range(1,len(self.st)):
                        ax1.axvspan(self.st[i],self.et[i],facecolor='0.5', alpha=0.2)
                    ax1.axvspan(self.st[0],self.et[0], facecolor='0.5',alpha=0.2,label='Satellite GTIs')

                # loop over all the different solutions
                _label = 'matched'
                for i, numburstssim in enumerate(times.keys()):
                    # print (i, numburstssim, times, ebs)
                    timepred = [x[0] for x in times[numburstssim]]
                    timepred_errup = [x[1] for x in times[numburstssim]]
                    timepred_errlow = [x[2] for x in times[numburstssim]]
                    ebpred = [x[0] for x in ebs[numburstssim]]
                    ebpred_errup = [x[1] for x in ebs[numburstssim]]
                    ebpred_errlow = [x[2] for x in ebs[numburstssim]]
                    ax1.errorbar(timepred[itoff:], ebpred,
                        yerr=[ebpred_errlow, ebpred_errup],
                        # xerr=[timepred_errup[1:], timepred_errlow[1:]],
                        marker='*', ms=11, linestyle='', color='C{}'.format(i),
                        label='predicted ({})'.format(numburstssim))

                    if self.continuous:
                        # only relevant for train mode with generate_burst_train
                        ref_tpred = np.argmin(np.abs(self.bstart[self.ref_ind]-timepred))
                        imatch = burst_time_match(self.ref_ind, self.bstart,
                            ref_tpred, np.array(timepred))
                    else:
                        imatch = burst_time_match(0, self.bstart,
                            0, np.array(timepred))

                    print (numburstssim, imatch)
                    # Not sure this will work so well if there are
                    # multiple sets of solutions
                    # if len(times.keys()) == 1:
                    imatchm1 = [x-itoff for x in imatch if x-itoff >= 0]
                    ax1.plot(np.array(timepred[itoff:])[imatchm1],
                        np.array(ebpred)[imatchm1],
                        marker='*', ms=5, linestyle='', color='tab:red',
                        label=_label,zorder=99)
                    _label = None # only give the label the first time

                    resid = -(self.bstart-np.array(timepred)[imatch])*24.
                    axs['resid'].errorbar(self.bstart, resid,
                        yerr=[np.array(timepred_errlow)[imatch]*24.,
                        np.array(timepred_errup)[imatch]*24.],
                        marker='*', ms=11, linestyle='', color='C{}'.format(i))
                    logger.info ('RMS obs-model offset ({}, {:.2f}%) = {:.4f} hr'.format(
                        numburstssim, 100.*self.model_pred['part_stats'][numburstssim]/len(self.samples),
                        np.sqrt(np.mean(resid**2))))

                ax1.set_ylabel("Fluence ($10^{-6}\\,{\\rm erg\\,cm^{-2}}$)")

                if title is not False:
                    ax1.set_title('beansp v{} run {} last {}/{}'.format(
                        self.version, self.run_id, self.nsteps_completed-self.samples_burnin, self.nsteps_completed) if title is None else title,loc='right')

                axs['resid'].axhline(0.0, color=obs_colour, ls='--')
                axs['resid'].set_ylabel('Time offset (hr)')
                axs['resid'].set_xlabel("Time (days after MJD {})".format(self.tref))
                ax1.legend(loc=2)

            else:
                # different style plot for the ensemble mode; the bstart
                # still records the burst time (epoch) but now we prefer
                # to plot vs. recurrence time
                timepred = [x[0] for x in times[self.numburstssim]]
                timepred_errup = [x[1] for x in times[self.numburstssim]]
                timepred_errlow = [x[2] for x in times[self.numburstssim]]
                ebpred = [x[0] for x in ebs[self.numburstssim]]
                ebpred_errup = [x[1] for x in ebs[self.numburstssim]]
                ebpred_errlow = [x[2] for x in ebs[self.numburstssim]]

                fig = plt.figure()
                plt.errorbar(self.tdel, self.fluen, xerr=self.tdele, yerr=self.fluene,
                    color='black', linestyle='', marker='.', ms=13, label='Observed')
                plt.scatter(timepred, ebpred, marker='*', color=bursts_colour, s=100, label='Predicted')
                plt.errorbar(timepred, ebpred,
                    yerr=[ebpred_errlow, ebpred_errup],
                    xerr=[timepred_errup, timepred_errlow], fmt='.',
                    color=bursts_colour)
                plt.xlabel("Recurrence time (hr)")

                plt.ylabel("Fluence ($10^{-6}\\,{\\rm erg\\,cm^{-2}}$)")
                plt.legend(loc=2)

            fig.show()

            if savefig:
                savefile = '{}_predictedburstscomparison.pdf'.format(self.run_id)
                logger.info ('saving burst comparison plot to {}'.format(savefile))
                fig.savefig(savefile)
            else:
                logger.info ('skipping burst comparison plot save')


    def compare(self, alt, burnin=None, parameters=None, label='result 2'):
        '''
        This method will update the cc attribute to include data from a
        different Beans object, or a different set of analyses, to be able
        to plot both posteriors simultaneously

        :param alt: object to compare with
        :param label: label to give the second set of results

        :returns: nothing, but updates the cc attribute
        '''

        # first re-define the main cc object
        self.cc = ChainConsumer()
        # augment with some additional parameters, e.g. as used in the
        # 1826 work
        # TODO probably should check if these have already been added
        self.cc_parameters['dxi_b'] = '$d\\xi_b$ (kpc)'
        self.cc_parameters['xi_p/xi_b'] = '$\\xi_p/\\xi_b$'
        _samples = np.column_stack((self.samples,
            self.samples[:,3]*self.samples[:,4],
            self.samples[:,5]/self.samples[:,4]))

        self.cc.add_chain(_samples,
            parameters=list(self.cc_parameters.values()),
            name=self.run_id)
        self.samples = _samples

        if type(alt) == Beans:
            # don't really need to worry about common parameters here,
            # just use whatever is in the other bean

            if not hasattr(alt, 'samples'):
                logger.error ('comparison object has no samples, run do_analysis first')
                return

            self.cc.add_chain(alt.samples,
                parameters=list(alt.cc_parameters.values()),
                name=label)

        elif type(alt) == str:
            # read in the parameters from a file
            if alt[-6:] == 'npy.gz':
                # this option will read in parameters from the files
                # distributed with Johnston et al. 2020 (MNRAS 494, 4576);
                # see https://data.mendeley.com/datasets/nmb24z6jrp/2
                f = gzip.GzipFile(alt, 'r')
                chain = np.load(f)
                logger.info ('read in array with {} walkers, {} steps and {} parameters from\n  {}'.format(*np.shape(chain), alt))
                # shape is (nwalkers, nsteps, ndim)
                assert np.shape(chain)[2] == 12 # won't work for the other files
                if burnin is not None:
                    chain_flat = chain[:, burnin:, :].reshape((-1, 12))
                else:
                    chain_flat = chain[:, :, :].reshape((-1, 12))
                self.cc.add_chain(chain_flat,
                    parameters=['$\dot{m}_1','$\dot{m}_2','$\dot{m}_3',
                        '$Q_{b,1}$ (MeV)', '$Q_{b,2}$ (MeV)', '$Q_{b,3}$ (MeV)',
                        '$X$', '$Z$', '$g$ (cm s$^{-2}$)', '$M$ ($M_\odot$)',
                        '$d\\xi_b$ (kpc)', '$\\xi_p/\\xi_b$'],
                    name=label)
            else:
                breakpoint()
        else:
            print ('** ERROR ** can only compare with another Beans object or chains read in from a file')

        # and finally set all the plot options (have to do this last)
        self.cc.configure(usetex=True, serif=True,
            flip=False, summary=False,
            bins=0.7, # has the effect of light smoothing of the histograms
            diagonal_tick_labels=False, max_ticks=3, shade=True, \
            shade_alpha=1.0 ,bar_shade=True, tick_font_size='xx-large', \
            label_font_size='xx-large',smooth=True, \
            sigma2d=False, sigmas=[1,2]) #np.linspace(0, 3, 4))
        self.cc_nchain = 2

        print ('Added posteriors for run_id {} to the ChainConsumer object; re-run do_analysis\n  to show the posterior comparison'.format(alt.run_id))


    def prune(self, key=None, nwalkers=None, scale=0.0, savefile=None):
        '''
        This method will "prune" the walkers to keep only one set of
        solutions. It is necessary that :meth:`Beans.do_analysis` has
        already been called, with the `comparison` option, to identify the
        model classes.

        The method will return the new set of positions, and optionally
        save them to a pickle file

        :param key: model class to keep, usually labeled by the number of predicted bursts
        :param nwalkers: number of samples to generate, if not the current number of walkers
        :param scale: in case it's necessary to distribute the walker positions around the seed values, this parameter sets the Gaussian scale
        :param savefile: name of pickle file to save to

        :returns: array of walker positions, dimensions (nwalkers, ndim)
        '''

        if not hasattr(self, 'model_pred'):
            print ('** ERROR ** no model predictions available, run the comparison first')
            return None

        # TODO: if no key is supplied, could just save all the positions
        if key is None:
            print ("** ERROR ** no key supplied, don't know which set to keep")
            return None

        if not (key in set(self.model_pred['partition'])):
            print ("** ERROR ** key not present in model prediction set: {}".format(set(self.model_pred['partition'])))
            return None

        new_pos, d1, d2, d3 = self.reader.get_last_sample()

        # now redistribute the "bad" walkers

        bad = np.array(self.model_pred['partition'])[-self.nwalkers:] != key

        print ('Redistributing {} of {} walkers...'.format(len(np.where(bad)[0]), self.nwalkers))

        for i in (np.where(bad)[0]):
            new_pos[i] = new_pos[np.random.choice(np.where(~bad)[0])] + scale*np.random.randn(self.ndim)

        if savefile is not None:
            print ('Saving positions to {}'.format(savefile))
            pickle.dump(new_pos, open(savefile, 'wb'))

        return new_pos


# end of beans.py
