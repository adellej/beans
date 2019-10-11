""" function to initiliase and run emcee. Please note this software uses the latest stable version of emcee (v2.2.1). """
import numpy as np
import emcee
import sys
import pickle

def runemcee(nwalkers, nsteps, ndim, theta, lnprob, x, y, yerr, threads, run_id):

    # This section now defines the initial walker positions and next defines the chain and runs emcee.

    pos = [theta + 1e-4 * np.random.randn(ndim) for i in range(nwalkers)]

    print("# -------------------------------------------------------------------------#")
    print("Ready to run", run_id, "with", len(pos), "walkers")

    # Define the sampler for emcee. Threads argument should be set to how many cores your computer has.
    # sampler.reset()
    sampler = emcee.EnsembleSampler(
        nwalkers, ndim, lnprob, args=(x, y, yerr), threads=threads
    )

    print("Running sampler...")
    # Option to run without progress bar:
    # run = sampler.run_mcmc(pos, 1000)

    # Option to run with progress bar:
    width = 60

    # open file to incrementally save progress in case run crashes 
    # f = open("chain.dat", "w")
    # f.close()
    
    for i, result in enumerate(sampler.sample(pos, iterations=nsteps, storechain=False)):
        n = int((width + 1) * float(i) / nsteps)
        sys.stdout.write("\r[{0}{1}]".format("#" * n, " " * (width - n)))

        # position = result[0]
        # f = open("chain.dat", "a")
        # for k in range(position.shape[0]):
        #     np.savetxt(f, np.column_stack((k, position[k])))
        # f.close()
    sys.stdout.write("\n")

    print("...sampler run " + run_id + " complete")

    # -------------------------------------------------------------------------#
    # Save multi-threaded sampler:
    print(f"Now saving the chains to file chains/{run_id}_chain.p")


    def __getstate__(sampler):
        sampler_dict = sampler.__dict__.copy()
        del sampler_dict["pool"]
        return sampler_dict


    sampler_nopool = __getstate__(sampler)
    pickle.dump(sampler_nopool, open("chains/" + run_id + "_chain.p", "wb"))


    def __setstate__(self, state):
        self.__dict__.update(state)

    print("All done!")

    return sampler

