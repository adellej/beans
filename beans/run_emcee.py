""" function to initiliase and run emcee. Please note this software uses the latest stable version of emcee (v2.2.1). """
import numpy as np
import emcee
import h5py

def runemcee(nwalkers, nsteps, ndim, theta, lnprob, x, y, yerr, run_id, restart):

    # This section now defines the initial walker positions and next defines the chain and runs emcee.

    print("# -------------------------------------------------------------------------#")

    # define the dtype of the blobs
    #dtype = [("lnprob", float), ("model", "S1000")]
    dtype = [("lnprob", float), ("model", h5py.string_dtype(encoding='ascii'))]



    if restart == True:
        #pos = [theta + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
        print('Restarting',run_id,'with',nwalkers,'walkers')
    else:
        # set the intial position of the walkers
        pos = [theta + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
        print('Ready to run',run_id,'with',nwalkers,'walkers')

    # with Pool() as pool: # for some reason the multi-threaded isn't working with the new emcee. Revert back to not using MPI.

    print("Beginning sampling..")
    # use emcee backend to save as a HD5 file
    reader = emcee.backends.HDFBackend(filename=run_id+".h5")
    if restart == False:
        reader.reset(nwalkers, ndim)
        # sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr), pool=pool, backend=reader, blobs_dtype=dtype, moves=emcee.moves.StretchMove(a=1.5))
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr), backend=reader, blobs_dtype=dtype, moves=emcee.moves.StretchMove(a=1.5))



    # We'll track how the average autocorrelation time estimate changes
    index = 0
    autocorr = np.empty(nsteps)

    # This will be useful to testing convergence
    old_tau = np.inf
    if restart == False:
        for sample in sampler.sample(pos, iterations=nsteps, progress=True, store=True):
            # Only check convergence every 100 steps
            if sampler.iteration % 100:
                continue

            # Compute the autocorrelation time so far
            # Using tol=0 means that we'll always get an estimate even
            # if it isn't trustworthy
            accept_frac = np.mean(sampler.acceptance_fraction)
            with open(f'{run_id}_acceptancefraction.txt', 'a') as f:
                np.savetxt(f, [[float(sampler.iteration)],[accept_frac]], delimiter=',', newline='\n')

            tau = sampler.get_autocorr_time(tol=0)
            autocorr[index] = np.mean(tau)
            index += 1

            # Check convergence
            converged = np.all(tau * 100 < sampler.iteration)
            converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
            if converged:
                print('Complete! Chains are converged')
                break
            old_tau = tau
    else:
        pos = sampler.get_last_sample()
        for sample in sampler.sample(pos, iterations=nsteps, progress=True, store=True):
            print(sampler.acceptance_fraction)
            if sampler.iteration % 100:
                continue
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
                print('Complete! Chains are converged')
            break
            old_tau = tau


        #tau = sampler.get_autocorr_time()
        print('Complete! WARNING max number of steps reached but chains may or may not be converged.')


    #     print('Samples complete. Took {0:.1f} seconds'.format(multi_time))

    # with Pool() as pool:
    #     print("Beginning sampling..")
    #     # use emcee backend to save as a HD5 file
    #     reader = emcee.backends.HDFBackend(f"chains_{run_id}.h5")
    #     reader.reset(nwalkers, ndim)
    #     sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr), pool=pool, backend=reader, blobs_dtype=dtype)
    #     start = time.time()
    #     sampler.run_mcmc(pos, nsteps, progress=True, store=True)
    #     end = time.time()
    #     multi_time = end-start

    #     print('Samples complete. Took {0:.1f} seconds'.format(multi_time))


    return sampler
