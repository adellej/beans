""" Reads data in from ascii files and puts in right format for emcee """

from astropy.io import ascii
import numpy as np

def get_obs(ref_ind, bc, gti_checking, obsname, burstname, gtiname):
    """
    Reads data in from ascii files and puts in right format for emcee
    Requires persistent flux observations, burst observations and satellite
    telescope gtis

    The following data is required, in ascii format, with columns in the given order:
    observations (start time (MJD), stop time (MJD) persistent flux measurements (in 3-25keV 1e-9 erg/cm^2/s), pflux error)
    burst observations (time (MJD), fluence (in 1e-9 erg/cm^2/s) fluence error, alpha, alpha error)
    The following data is required in a tab separated file with just 2 columns:
    satellite gti observing data (start time of obs, stop time of obs (both in MJD))

    bc is the bolometric correction to convert the persistent fluxes from 3-25keV into bolometric fluxes.

    Example screenshots of how to generate this ascii format using the MINBAR web interface for SAX J1808.4-3658 is included in this directory.

    """

    global st, et

    print("Reading in observation data files ...")

    # Read in the observation data that contains initial conditions and burst observation parameters:


    # -------------------------------------------------------------------------#
    # Need len(tobs) to intialise emcee:
    # Get the observing times and peak flux arrays:
    if obsname is not 'none':
        obsdata = ascii.read(obsname)
        ta_1 = obsdata['col1']
        ta_2 = obsdata['col2']

        ssa_1 = ta_1
        ssa_2 = ta_2

        tobs = 0.5*(ta_1 + ta_2) # get the average of the start and stop times
        tobs = np.array(tobs)
        tobs_err = 0.5*(ta_2-ta_1) #calculate an error
        tobs_err = np.array(tobs_err)
        pflux = obsdata['col3']
        pflux = np.array(pflux)
        pfluxe = obsdata['col4']
        pfluxe = np.array(pfluxe)
    else:
        burstdata = ascii.read(burstname)
        tobs = burstdata['col1']
        tobs_err = np.ones(len(tobs)) * 0.5
        pflux = np.array(burstdata['col6'])
        pfluxe = np.array(burstdata['col7'])
    iobs = np.arange(len(tobs))
    good = iobs
    # Check the arrays are sorted here
    _i = np.argsort(tobs)
    if not np.all(_i == iobs):
        print ('** WARNING ** input observation data is not sorted, fixing')
        tobs = tobs[_i]
        pflux = pflux[_i]
        pfluxe = pfluxe[_i]
        good = good[_i]

    # -------------------------------------------------------------------------#

    if burstname is not None:
        burstdata = ascii.read(burstname)

        # Get the burst time, fluence and alpha arrays:

        bstart = burstdata['col1']
        bstart = np.array(bstart)
        bstart_err = np.full(len(bstart), 0.0069)
        fluen = burstdata['col2']
        fluen = np.array(fluen)
        fluene = burstdata['col3']
        fluene = np.array(fluene)
        alpha = np.array(burstdata['col4'])
        alphae = np.array([burstdata['col5']])

        # If alphas are pre-calculated, read them in here. If not the initialisation routine will do it below
        #alpha = bursts['col9']
        #alpha = np.array(alpha)
        #alpha_err = bursts['col10']
        #alpha_err = np.array(alpha_err)

        # Define reference time as start of first burst
        bstart0 = bstart[0]
        bstart = bstart - bstart0

        # Define reference burst:
        tref = bstart[ref_ind]

        # Assemble obs and obs_err arrays

        obs = []
        obs.extend(bstart)
        obs.extend(fluen)
        obs.extend(alpha[1:])  # always exclude the first alpha because we do not know the recurrence time of the first burst to calculate this
        obs_err = []
        obs_err.extend(bstart_err)
        obs_err.extend(fluene)
        obs_err.extend(alphae[0][1:])

        del obs[ref_ind]  # delete the time of the reference burst because we do not model for this
        del obs_err[ref_ind]

    else:
        print ('** WARNING ** skipping read of burst data, assuming no bursts observed')
        # if there's no burst data, define the start time as the first observation instead

        bstart, fluen, obs, obs_err = None, None, None, None
        # bstart0 = min(tobs)
        bstart0 = tobs[np.argmax(pflux)]

    # Shift the observations and gtis so that they start at bstart0:

    tobs = tobs - bstart0

    if gti_checking == 1:

        # Read in the gtis (required arrays are st (start time) and et (end time) of times telescope IS observing (indexes need to match)
        # gtis should be in MJD
        gtidata = np.loadtxt(gtiname)

        # Now extract and scale gti data:
        st = gtidata[:,0]
        st = np.array(st)
        et = gtidata[:,1]
        et = np.array(et)
        st = st-bstart0
        et = et-bstart0

    pflux = pflux[good]
    pfluxe = pfluxe[good]

    # Bolometric correction:

    pflux=pflux*bc
    pfluxe=pfluxe*bc

    if burstname is not None:
        if gti_checking == 1:
            print(f"Observations found: \n burst times = {bstart}\n fluences = {fluen}\n alphas = {alpha}\n persistent fluxes = {pflux}\n You have chosen gti_checking = {gti_checking}, so the gti data is st = {st}")

        else:
            print(f"Observations found: \n burst times = {bstart}\n fluences = {fluen}\n alphas = {alpha}\n persistent fluxes = {pflux}\n You have chosen gti_checking = {gti_checking}, so no gti data will be used.")
    else:
        if gti_checking == 0:
            print(f"Observations found: \n persistent fluxes = {pflux}\n You have chosen gti_checking = {gti_checking}, so no gti data will be used.")
            # Need the GTIs to run the simulations in the absence of bursts, so should let the user know
            # here
            pass
        else:
            print(f"Observations found: \n persistent fluxes = {pflux}\n You have chosen gti_checking = {gti_checking}, so the gti data is st = {st}")

    if gti_checking == 1:
        return bstart0, bstart, fluen, obs, obs_err, pflux, pfluxe, tobs, st, et
    else:
        return bstart0, bstart, fluen, obs, obs_err, pflux, pfluxe, tobs
