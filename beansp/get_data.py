""" Reads data in from ascii files and puts in right format for emcee """

from astropy.io import ascii
import numpy as np
from os.path import exists
import sys

def get_obs(beans_obj):
    """
    Reads data in from ascii files and returns the data structures
    required for emcee Requires persistent burst observations and,
    optionally, flux observations & satellite telescope GTIs; or just
    persistent flux measurements and GTIs.

    The following data is required, in tab-separated ascii format, with
    columns in the given order: burst observations:
        time (MJD), fluence (in 1e-9 erg/cm^2/s) fluence error, alpha, alpha error
    and optionally for ensemble mode,
        persistent flux (in 1e-9 erg/cm^2/s) and error, and recurrence time (hr) and error

    observations:
        start time (MJD), stop time (MJD) persistent flux measurements (in 3-25keV 1e-9 erg/cm^2/s), pflux error

    satellite gti observing data:
        start time of obs, stop time of obs (both in MJD)

    bc is the bolometric correction to convert the persistent fluxes (e.g.
    from the 3-25keV range used for MINBAR) into bolometric estimates.

    Example screenshots of how to generate this ascii format using the
    MINBAR web interface for SAX J1808.4-3658 is included in this directory.

    In the earlier version of this routine we passed the parameters
    ref_ind, bc, obsname, burstname, gtiname, and returned x, y, yerr,
    tref, bstart, pflux, pfluxe, tobs, fluen, fluen_err, and optionally
    st, et; but now we just get the required parameters from the Beans
    object (beans_obj), and set the attributes directly

    :param beans_obj: Beans object from which ref_ind, bc, obsname,
      burstname & gtiname will be used to find the input data

    :return: nothing
    """

    print("Reading input data files ...")

    # Read in the burst and observation data that contains the measurements to match

    if beans_obj.burstname is not None:
        # if not exists(burstname):
        #     print("** ERROR ** burst data file {} not found".format(burstname))
        #     sys.exit()

        # The "ensemble" mode tables, with additional columns, seem to require these additional
        # parameters to guarantee they are read in correctly
        burstdata = ascii.read(beans_obj.burstname, format='tab', header_start=None, data_start=0)
        assert (len(burstdata.columns) >=5) or ((beans_obj.obsname is None) & (len(burstdata.columns) >=9))

        # Get the burst time, fluence and alpha arrays:

        beans_obj.bstart = np.array(burstdata['col1'])
        beans_obj.fluen = np.array(burstdata['col2'])
        beans_obj.fluene = np.array(burstdata['col3'])
        alpha = np.array(burstdata['col4'])
        alphae = np.array(burstdata['col5'])

        # Define reference time as start of first burst/epoch
        # bstart0 = bstart[0]
        beans_obj.tref = beans_obj.bstart[0]

    else:
        print ('** WARNING ** skipping read of burst data, assuming no bursts observed')
        # if there's no burst data, define the start time as the first observation instead

    # -------------------------------------------------------------------------#
    # Need len(tobs) to intialise emcee:
    # Get the observing times and peak flux arrays:
    if beans_obj.obsname is not None:
        # if not exists(obsname):
        #     sys.exit("** ERROR ** observation file {} not found".format(obsname))
        obsdata = ascii.read(beans_obj.obsname)
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

        beans_obj.pflux = pflux * beans_obj.bc
        beans_obj.pfluxe = pfluxe * beans_obj.bc

        if beans_obj.burstname is not None:
            # We're using the start time in the MCMC so have to assign an error
            # This parameter is not returned to the init func, except
            # through the obs_err -> yerr parameters
            # bstart_err = np.full(len(beans_obj.bstart), 0.0069)
            bstart_err = np.full(len(beans_obj.bstart), beans_obj.bstart_err)
            # bstart = bstart - bstart0
            beans_obj.bstart = beans_obj.bstart - beans_obj.tref

            # Define reference burst:
            tref = beans_obj.bstart[beans_obj.ref_ind]

	    # Set the y and yerr arrays (previously obs, obs_err), which
	    # are what the MCMC code uses for comparison

            beans_obj.y = np.concatenate((beans_obj.bstart, beans_obj.fluen, alpha[1:]), axis=0)
            # always exclude the first alpha because we do not know
            # the recurrence time of the first burst to calculate this
            beans_obj.yerr = np.concatenate((bstart_err, beans_obj.fluene, alphae[1:]), axis=0)

            beans_obj.y = np.delete(beans_obj.y, beans_obj.ref_ind)  # delete the time of the reference burst because we do not model for this
            beans_obj.yerr = np.delete(beans_obj.yerr, beans_obj.ref_ind)
        else:
            beans_obj.tref = tobs[np.argmax(beans_obj.pflux)]
            # bstart0 = tobs[np.argmax(pflux)]

            beans_obj.bstart, beans_obj.fluen, beans_obj.y, beans_obj.yerr = None, None, None, None
            # bstart0 = min(tobs)

    else:
        # If there's no observation data that indicates an "ensemble" mode run
	# Here we define some additional/alternate parameters for the
	# "ensemble" mode run
        # Note that no bolometric correction is applied to the persistent
        # fluxes

        tobs = burstdata['col1']
        tobs_err = np.ones(len(tobs)) * 0.5
        beans_obj.pflux = np.array(burstdata['col6'])
        beans_obj.pfluxe = np.array(burstdata['col7'])
        tdel = np.array(burstdata['col8'])
        tdele = np.array(burstdata['col9'])

	# Set the y and yerr arrays (previously obs and obs_err), which
	# are what the MCMC code uses for comparison

        beans_obj.y = np.concatenate((tdel, beans_obj.fluen, alpha), axis=0)
        beans_obj.yerr = np.concatenate((tdele, beans_obj.fluene, alphae), axis=0)

    # -------------------------------------------------------------------------#

    print ('...done')

    # Shift the observations and gtis so that they start at bstart0:

    # tobs = tobs - bstart0
    beans_obj.tobs = tobs - beans_obj.tref

    gti_checking = beans_obj.gtiname is not None
    if gti_checking:

        # Read in the gtis (required arrays are st (start time) and et (end time) of times telescope IS observing (indexes need to match)
        # gtis should be in MJD
        # if not exists(gtidata):
        #     sys.exit("** ERROR ** GTI data file {} not found".format(gtidata))
        gtidata = np.loadtxt(beans_obj.gtiname)

        # Now extract and scale gti data:
        beans_obj.st = np.array(gtidata[:,0])-beans_obj.tref
        beans_obj.et = np.array(gtidata[:,1])-beans_obj.tref
        # st = st-bstart0
        # et = et-bstart0

        if beans_obj.burstname is None:
            print(f"""
Observations found:
  persistent fluxes = {beans_obj.pflux}
  gti data file = {beans_obj.gtiname}
GTI checking will be performed""")
        else:
            print(f"""
Observations found:
  burst times = {beans_obj.bstart}
  fluences = {beans_obj.fluen}
  alphas = {alpha}
  persistent fluxes = {beans_obj.pflux}
  gti data file = {beans_obj.gtiname}
GTI checking will be performed""")

        # not sure if you need to return st, et if they're also global variables (defined at start)
        return # bstart0, bstart, fluen, fluene, obs, obs_err, pflux, pfluxe, tobs, st, et

    else:
        if beans_obj.burstname is None:
            print(f"""
Observations found:
  persistent fluxes = {beans_obj.pflux}
You have not supplied a GTI file, so no GTI checking will be performed.""")
        else:
            print(f"""
Observations found:
  burst times = {beans_obj.bstart}
  fluences = {beans_obj.fluen}
  alphas = {alpha}
  persistent fluxes = {beans_obj.pflux}
You have not supplied a GTI file, so no GTI checking will be performed.""")

        beans_obj.st, beans_obj.et = None, None

    return # bstart0, bstart, fluen, fluene, obs, obs_err, pflux, pfluxe, tobs
