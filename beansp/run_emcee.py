import numpy as np
import emcee
import sys
import pickle
from multiprocessing import Pool
import time
import h5py

def set_initial_positions(theta, nwalkers, prior, scale=1e-4):
    '''
    Distribute the initial positions for the walkers around the supplied theta value, 
    making sure they're consistent with the prior

    :param theta: desired centroid for walker parameters
    :param nwalkers: number of initial positions to generate
    :param prior: prior function adopted
    :param scale: Gaussian spread of the values in each dimension

    :return: list of nwalkers positions
    '''
    
    ndim = len(theta)
    pos = np.array([theta + scale*np.random.randn(ndim) for i in range(nwalkers)])
    
    valid = np.array([prior(pos[i]) != -np.inf for i in range(nwalkers)])
    
    print ('Initial walker positions within {} of supplied parameter vector, checking for consistency with prior...'.format(scale))

    # print (len(np.where(~valid)[0]))
    while (len(np.where(~valid)[0])) > 0:
        pos[~valid] = [theta + scale*np.random.randn(ndim) for i in range(len(np.where(~valid)[0]))]
        valid = np.array([prior(pos[i]) != -np.inf for i in range(nwalkers)])
        # print (len(np.where(~valid)[0]))                                                                     

    return list(pos)


def runemcee(nwalkers, nsteps, theta, lnprob, prior, x, y, yerr, run_id,
    restart, threads, stretch_a, pos=None, **kwargs):
    """
    Function to initilise and run emcee.
    Removed the redundant parameter ndim, which can be determined from the
    theta parameter array

    :param nwalkers: number of walkers for the emcee run
    :param nsteps: number of MCMC steps to run
    :param theta_in: parameter tuple, with *X*, *Z*, *Q_b*, *d*, *xi_b*,
      *xi_p*, and (optionally) *mass*, *radius*, *f_E* & *f_a*
    :param lnprob: log-probability function to use with emcee
    :param prior: prior function to use with emcee, used to check the initial walker positions
    :param x: the "independent" variable, passed to lnlike
    :param y: the "dependent" variable (i.e. measurements), passed to lnlike
    :param yerr: erorr estimates on y
    :param run_id: string identifier for the run, used to label all the
      result files, and where you want output to be saved
    :param restart: set to True to continue a previously interrupted run
    :param threads: number of threads for emcee to use (e.g. number of
      cores your computer has). Set to None to use all available
    :param stretch_a: the Goodman & Weare (2010) stretch move scale parameter
    """

    ndim = len(theta)

    # This section now defines the initial walker positions and next defines the chain and runs emcee.

    print("\n# ---------------------------------------------------------------------------#")

    # define the dtype of the blobs
    # this refers to the 2nd and subsequent parameters returned by lnprob;
    # in our case the prior probability and the full model result, encoded
    # as ASCII

    dtype = [("lnprior", float), ("model", h5py.string_dtype(encoding='ascii'))]

    # use emcee backend to save as a HD5 file
    # see https://emcee.readthedocs.io/en/stable/user/backends
    reader = emcee.backends.HDFBackend(filename=run_id + ".h5")

    if restart == True:
        steps_so_far = np.shape(reader.get_chain())[0]
	# logic here is, that if nsteps > steps_so_far, then you're
	# restarting to try to complete the original nsteps target
        if nsteps-steps_so_far > 0:
            print('    Continuing run={} with {} walkers from {}/{} steps...'.format(run_id,nwalkers,steps_so_far,nsteps))
            nsteps -= steps_so_far
        else:
            print('    Extending run={} with {} walkers after {} steps done...'.format(run_id,nwalkers,steps_so_far))
    else:
        reader.reset(nwalkers, ndim)

        print('    Running run={} with {} walkers, target {} steps...'.format(run_id,nwalkers,nsteps))

    # implement multiprocessing following the simple example here:
    # https://emcee.readthedocs.io/en/stable/tutorials/parallel
    with Pool(processes=threads) as pool:

        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob,
            args=(x, y, yerr), backend=reader, blobs_dtype=dtype, pool=pool,
            moves=emcee.moves.StretchMove(a=stretch_a),
            **kwargs)
    
        # We'll track how the average autocorrelation time estimate changes
        index = 0
        autocorr = np.empty(nsteps)

        # This will be useful to testing convergence
        old_tau = np.inf

        # Earlier we had 2 identical sections here for the 2 restart
	# options, but the only difference was the choice of the starting
	# positions, so have now merged them

        if pos is not None:
            # warning is probably not necessary since we've already vetted
            # the positions in beans.do_run
            print ('\n    ** WARNING ** setting walkers at provided position vector')
        elif restart == True:
            pos = sampler.get_last_sample()
        else:
            # set the intial position of the walkers
            # pos = [theta + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
            pos = set_initial_positions(theta,  nwalkers, prior)

        print("# ---------------------------------------------------------------------------#")
            
        for sample in sampler.sample(pos, iterations=nsteps, progress=True):

            # Only check convergence every 100 steps
            if sampler.iteration % 100:
                continue
    
            # accumulate the acceptance fraction in the file
            accept_frac = np.mean(sampler.acceptance_fraction)
            with open(f'{run_id}_acceptancefraction.txt', 'a') as f:
                np.savetxt(f, np.c_[float(sampler.iteration),accept_frac],
                    fmt='%1.4e, %1.4f')
    
            # Compute the autocorrelation time so far
            # Using tol=0 means that we'll always get an estimate even
            # if it isn't trustworthy
            tau = sampler.get_autocorr_time(tol=0)
            autocorr[index] = np.mean(tau)
            index += 1
    
            # Check convergence
            converged = np.all(tau * 100 < sampler.iteration)
            converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
            if converged:
                print('...sampling complete! Chains are converged')
                break
            old_tau = tau

        print('...sampling complete! ** WARNING ** chains may not be converged.')

    return sampler

