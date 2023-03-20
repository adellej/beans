""" Reads data in from ascii files and puts in right format for emcee """

from astropy.io import ascii
import numpy as np
from os.path import exists
import sys

def get_obs(ref_ind, bc, gti_checking, obsname, burstname, gtiname):
    """
    Reads data in from ascii files and returns the data structures required for emcee
    Requires persistent burst observations and, optionally, flux observations & satellite
    telescope GTIs; or just persistent flux measurements and GTIs.

    The following data is required, in tab-separated ascii format, with columns in the given order:
    burst observations:
        time (MJD), fluence (in 1e-9 erg/cm^2/s) fluence error, alpha, alpha error; and optionally for ensemble mode,
            persistent flux (in 1e-9 erg/cm^2/s) and error, and recurrence time (hr) and error

    observations:
        start time (MJD), stop time (MJD) persistent flux measurements (in 3-25keV 1e-9 erg/cm^2/s), pflux error

    satellite gti observing data:
        start time of obs, stop time of obs (both in MJD)

    bc is the bolometric correction to convert the persistent fluxes from 3-25keV into bolometric fluxes.

    Example screenshots of how to generate this ascii format using the MINBAR web interface for SAX J1808.4-3658
    is included in this directory.
    """

    # Is this statement redundant? st and et are returned when the GTI file is read in
    global st, et

    print("Reading input data files ...")

    # Read in the burst and observation data that contains the measurements to match

    if burstname is not None:
        # if not exists(burstname):
        #     print("** ERROR ** burst data file {} not found".format(burstname))
        #     sys.exit()

        # The "ensemble" mode tables, with additional columns, seem to require these additional
        # parameters to guarantee they are read in correctly
        burstdata = ascii.read(burstname, format='tab', header_start=None, data_start=0)
        assert (len(burstdata.columns) >=5) or ((obsname is None) & (len(burstdata.columns) >=9))

        # Get the burst time, fluence and alpha arrays:

        bstart = np.array(burstdata['col1'])
        fluen = np.array(burstdata['col2'])
        fluene = np.array(burstdata['col3'])
        alpha = np.array(burstdata['col4'])
        alphae = np.array(burstdata['col5'])

        # Define reference time as start of first burst/epoch
        bstart0 = bstart[0]

    else:
        print ('** WARNING ** skipping read of burst data, assuming no bursts observed')
        # if there's no burst data, define the start time as the first observation instead

    # -------------------------------------------------------------------------#
    # Need len(tobs) to intialise emcee:
    # Get the observing times and peak flux arrays:
    if obsname is not None:
        # if not exists(obsname):
        #     sys.exit("** ERROR ** observation file {} not found".format(obsname))
        obsdata = ascii.read(obsname)
        ta_1 = obsdata['col1']
        ta_2 = obsdata['col2']

        ssa_1 = ta_1
        ssa_2 = ta_2

        tobs = 0.5*(ta_1 + ta_2) # get the average of the start and stop times
        tobs = np.array(tobs)
        tobs_err = 0.5*(ta_2-ta_1) #calculate an error
        tobs_err = np.array(tobs_err)
        pflux = np.array(obsdata['col3'])
        pfluxe = np.array(obsdata['col4'])

        # Check the arrays are sorted here
        _i = np.argsort(tobs)
        if not np.all(_i == np.arange(len(tobs))):
            print('** WARNING ** input observation data is not sorted, fixing')
            tobs = tobs[_i]
            pflux = pflux[_i]
            pfluxe = pfluxe[_i]

        # Bolometric correction:

        pflux = pflux * bc
        pfluxe = pfluxe * bc

        if burstname is not None:
            # We're using the start time in the MCMC so have to assign an error
            bstart_err = np.full(len(bstart), 0.0069)
            bstart = bstart - bstart0

            # Define reference burst:
            tref = bstart[ref_ind]

            # Set the obs and obs_err arrays, which are what the MCMC code uses for comparison

            obs = np.concatenate((bstart, fluen, alpha[1:]), axis=0)
                                                # always exclude the first alpha because we do not know
                                                # the recurrence time of the first burst to calculate this
            obs_err = np.concatenate((bstart_err, fluene, alphae[1:]), axis=0)

            obs = np.delete(obs, ref_ind)  # delete the time of the reference burst because we do not model for this
            obs_err = np.delete(obs_err, ref_ind)
        else:
            bstart0 = tobs[np.argmax(pflux)]

            bstart, fluen, obs, obs_err = None, None, None, None
            # bstart0 = min(tobs)

    else:
        # If there's no observation data that indicates an "ensemble" mode run
        # Here we define some additional/alternate parameters for the "ensemble" mode run

        tobs = burstdata['col1']
        tobs_err = np.ones(len(tobs)) * 0.5
        pflux = np.array(burstdata['col6'])
        pfluxe = np.array(burstdata['col7'])
        tdel = np.array(burstdata['col8'])
        tdele = np.array(burstdata['col9'])

        # Set the obs and obs_err arrays, which are what the MCMC code uses for comparison

        obs = np.concatenate((tdel, fluen, alpha), axis=0)
        obs_err = np.concatenate((tdele, fluene, alphae), axis=0)

    # -------------------------------------------------------------------------#

    print ('...done')

    # Shift the observations and gtis so that they start at bstart0:

    tobs = tobs - bstart0

    if gti_checking:

        # Read in the gtis (required arrays are st (start time) and et (end time) of times telescope IS observing (indexes need to match)
        # gtis should be in MJD
        # if not exists(gtidata):
        #     sys.exit("** ERROR ** GTI data file {} not found".format(gtidata))
        gtidata = np.loadtxt(gtiname)

        # Now extract and scale gti data:
        st = np.array(gtidata[:,0])
        et = np.array(gtidata[:,1])
        st = st-bstart0
        et = et-bstart0

        if burstname is None:
            print(f"Observations found: \n persistent fluxes = {pflux}\n You have chosen gti_checking = {gti_checking}, so the gti data is st = {st}")
        else:
            print(f"Observations found: \n burst times = {bstart}\n fluences = {fluen}\n alphas = {alpha}\n persistent fluxes = {pflux}\n You have chosen gti_checking = {gti_checking}, so the gti data is st = {st}")

        # not sure if you need to return st, et if they're also global variables (defined at start)
        return bstart0, bstart, fluen, obs, obs_err, pflux, pfluxe, tobs, st, et

    else:
        if burstname is None:
            print(f"Observations found: \n persistent fluxes = {pflux}\n You have chosen gti_checking = {gti_checking}, so no gti data will be used.")
        else:
            print(f"Observations found: \n burst times = {bstart}\n fluences = {fluen}\n alphas = {alpha}\n persistent fluxes = {pflux}\n You have chosen gti_checking = {gti_checking}, so no gti data will be used.")

    return bstart0, bstart, fluen, obs, obs_err, pflux, pfluxe, tobs
