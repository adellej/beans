import numpy as np
import bilby
from beansp import beans,run_emcee

# This function is designed to implement sampling via bilby (see
#   https://github.com/bilby-dev/bilby
#
# Implemented initially with bilby v. 2.2.2.1

class BeansLikelihood(bilby.Likelihood):
    def __init__(self, data):
        """
        Adapted from the very simple Gaussian likelihood

        :param data: a Beans object
        """
        self.ndim = data.ndim
        if self.ndim == 6:
            super().__init__(parameters={"X": None, "Z": None, "Q_b": None, 
                                         # "f_a": None, "f_E": None, 
                                         # "r1": None, "r2": None, "r3": None, 
                                         "d": None, "xi_b": None, "xi_p": None})
        elif self.ndim == 8:
            super().__init__(parameters={"X": None, "Z": None, "Q_b": None, 
                                         "d": None, "xi_b": None, "xi_p": None,
                                         "M": None, "R": None})
        else:
            print ('** ERROR ** likelihood not implemented for ndim!=(6,8)')

        self.data = data # this is a beans object
        self.N = len(data.y) # don't know if this is necessary

        
    def log_likelihood(self):
        # here we use Beans.lnlike; this returns the likelihood AND the model, so have
        # to select only the first element
        
        X = self.parameters["X"]
        Z = self.parameters["Z"]
        Q_b = self.parameters["Q_b"]
        # f_a = self.parameters["f_a"]
        # f_E = self.parameters["f_E"]
        # r1 = self.parameters["r1"]
        # r2 = self.parameters["r2"]
        # r3 = self.parameters["r3"]
        d = self.parameters["d"]
        xi_b = self.parameters["xi_b"]
        xi_p = self.parameters["xi_p"]

        if self.ndim == 7:
            M = self.parameters["M"]
            theta = (X, Z, Q_b, d, xi_b, xi_p, M)
        elif self.ndim == 8:
            R = self.parameters["R"]
            theta = (X, Z, Q_b, d, xi_b, xi_p, M, R)
        else:
            theta = (X, Z, Q_b, d, xi_b, xi_p)

        # lnlike is the model likelihood only
        return self.data.lnlike(theta, None, self.data.y, self.data.yerr)[0]


def runbilby(bean, outdir='bilby_out', sampler='emcee'):
    """
    Function to implement sampling via bilby

    :param bean: Beans object, on which to do the sampling
    :param outdir: output directory

    :return:
    """

    bilby.utils.check_directory_exists_and_if_not_mkdir(outdir)

    theta = bean.theta

    priors = dict(
        # adapted from prior_func, for X, Z, Q_b, f_a, f_E, r1, r2, r3, mass, radius = theta_in
        X=bilby.core.prior.Uniform(1e-5, 0.76, "X"),
        Z=bilby.core.prior.Uniform(1e-5, 0.056, "Z"),
        Q_b=bilby.core.prior.Uniform(1e-6, 5.0, "Q_b"), # MeV/nucleon
        # f_a=bilby.core.prior.Uniform(1., 100, "f_a"),
        # f_E=bilby.core.prior.Uniform(1., 100, "f_E"),
        # r1=bilby.core.prior.Uniform(5e-3, 1.0, "r1"),
        # r2=bilby.core.prior.Uniform(5e-3, 3.0, "r2"),
        # r3=bilby.core.prior.Uniform(0, 1., "r3"),
        d=bilby.core.prior.Uniform(1,20,"d"),
        xi_b=bilby.core.prior.Uniform(0.1,2,"xi_b"),
        xi_p=bilby.core.prior.Uniform(0.1,2,"xi_p"),
    )
    if len(theta) > 6:
        M=bilby.core.prior.Uniform(1.15, 2.5, "M"), # solar masses

    if len(theta) > 7:
        R=bilby.core.prior.Uniform(9, 17, "R"), # km

    # set initial positions
    pos_i = run_emcee.set_initial_positions(theta, bean.nwalkers, beans.prior_func)

    # define likelihood
    likelihood = BeansLikelihood(bean)

    # run sampler
    result = bilby.run_sampler(
        likelihood=likelihood,
        priors=priors,
        pos0 = pos_i,
        outdir=outdir,
        label=bean.run_id,
        # bilby_mcmc:
        # sampler="bilby_mcmc",
        # nsamples=100,
        # emcee settings:
        sampler=sampler,
        nwalkers=bean.nwalkers, nsteps=bean.nsteps,
        # dynesty:
        # nlive=1000,
        npool=bean.threads, # 4 threads?
    )

    return result
