""" Reads data in from ascii files and puts in right format for emcee """

from astropy.io import ascii
import numpy as np
from os.path import exists
import sys

def get_obs(bean, alpha=True, fluen=True):
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
    object (bean), and set the attributes directly

    :param bean: Beans object from which ref_ind, bc, obsname,
      burstname & gtiname will be used to find the input data
    :param alpha: set to True (default) to include the alphas in the
      data for comparison, or False to omit
    :param fluen: set to True (default) to include the fluences in the
      data for comparison, or False to omit

    :return: nothing (everything comes back via the Beanobject)
    """

    print("Reading input data files ...")

    # Read in the burst and observation data that contains the measurements to match

    if bean.burstname is not None:
        # if not exists(burstname):
        #     print("** ERROR ** burst data file {} not found".format(burstname))
        #     sys.exit()

        # The "ensemble" mode tables, with additional columns, seem to require these additional
        # parameters to guarantee they are read in correctly
        burstdata = ascii.read(bean.burstname, format='tab', header_start=None, data_start=0)
        assert (len(burstdata.columns) >=5) or ((bean.obsname is None) & (len(burstdata.columns) >=9))

        # Get the burst time, fluence and alpha arrays:

        bean.bstart = np.array(burstdata['col1'])
        bean.tdel = (bean.bstart[1:]-bean.bstart[:-1])*24.
        bean.numburstsobs = len(bean.bstart)
        bean.cmpr_fluen = fluen
        bean.fluen = np.array(burstdata['col2'])
        bean.fluene = np.array(burstdata['col3'])
        bean.ifluen = bean.fluen > 0.
        bean.cmpr_alpha = alpha
        bean.alpha = np.array(burstdata['col4'])
        bean.alphae = np.array(burstdata['col5'])
        if np.any(~bean.ifluen):
            if np.any(bean.alpha[~bean.ifluen] > 0.):
                print('** WARNING ** nonzero alphas despite missing fluences will be ignored in fit')
        bean.ifluen = np.where(bean.ifluen)[0]

    else:
        print ('** WARNING ** skipping read of burst data, assuming no bursts observed')
        # if there's no burst data, define the start time as the first observation instead

    # -------------------------------------------------------------------------#
    # Need len(tobs) to intialise emcee:
    # Get the observing times and peak flux arrays:
    if bean.obsname is not None:
        # if not exists(obsname):
        #     sys.exit("** ERROR ** observation file {} not found".format(obsname))
        obsdata = ascii.read(bean.obsname)# format='tab', header_start=None, data_start=0)
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

        bean.pflux = pflux * bean.bc
        bean.pfluxe = pfluxe * bean.bc

        if bean.burstname is not None:

            # Define reference time as day of first burst/epoch, or day of
            # first observation, unless it's already been defined (e.g. read
            # in from the config file)
            # Changed this to be the day of the first burst/obs instead of the
            # time of the first burst from v2.11.0 onwards, so older runs need
            # to specify the value in the config file
            # bstart0 = bstart[0]

            if not hasattr(bean, 'tref'):
                bean.tref = np.min([np.floor(bean.bstart[0]),
                    np.floor(np.min(ta_1))])

            # We're using the start time in the MCMC so have to assign an error
            # This parameter is not returned to the init func, except
            # through the obs_err -> yerr parameters
            # bstart_err = np.full(len(beans_obj.bstart), 0.0069)
            bstart_err = np.full(len(bean.bstart), bean.bstart_err)
            # bstart = bstart - bstart0
            bean.bstart = bean.bstart - bean.tref

	    # Set the y and yerr arrays (previously obs, obs_err), which
	    # are what the MCMC code uses for comparison
            # always exclude the first alpha because we do not know
            # the recurrence time of the first burst to calculate this

            bean.y = bean.bstart
            bean.yerr = bstart_err
            if bean.cmpr_fluen:
                bean.y = np.concatenate((bean.y, bean.fluen[bean.ifluen]), axis=0)
                bean.yerr = np.concatenate((bean.yerr, bean.fluene[bean.ifluen]), axis=0)
            if bean.cmpr_alpha:
                bean.y = np.concatenate((bean.y, bean.alpha[1:][bean.ifluen[1:]-1]), axis=0)
                bean.yerr = np.concatenate((bean.yerr, bean.alphae[1:][bean.ifluen[1:]-1]), axis=0)

            bean.y = np.delete(bean.y, bean.ref_ind)  # delete the time of the reference burst because we do not model for this
            bean.yerr = np.delete(bean.yerr, bean.ref_ind)
        else:
            # if there's no burst data, define the tref in a different way
            # (time of peak of outburst)
            if not hasattr(bean, 'tref'):
                bean.tref = tobs[np.argmax(bean.pflux)]
            # bstart0 = tobs[np.argmax(pflux)]

            bean.bstart, bean.fluen, bean.y, bean.yerr = None, None, None, None
            # bstart0 = min(tobs)

    else:
        # If there's no observation data that indicates an "ensemble" mode run
	# Here we define some additional/alternate parameters for the
	# "ensemble" mode run
        # Note that no bolometric correction is applied to the persistent
        # fluxes

        tobs = burstdata['col1']
        tobs_err = np.ones(len(tobs)) * 0.5
        bean.pflux = np.array(burstdata['col6'])
        bean.pfluxe = np.array(burstdata['col7'])
        bean.tdel = np.array(burstdata['col8'])
        bean.tdele = np.array(burstdata['col9'])

        if not hasattr(bean, 'tref'):
            bean.tref = 0.

	# Set the y and yerr arrays (previously obs and obs_err), which
	# are what the MCMC code uses for comparison

        if bean.cmpr_fluen:
            bean.y = np.concatenate((bean.tdel, bean.fluen), axis=0)
            bean.yerr = np.concatenate((bean.tdele, bean.fluene), axis=0)
        if bean.cmpr_alpha:
            bean.y = np.concatenate((bean.y, bean.alpha), axis=0)
            bean.yerr = np.concatenate((bean.yerr, bean.alphae), axis=0)

    bean.ly = len(bean.y)
    assert bean.ly == len(bean.yerr)

    # -------------------------------------------------------------------------#

    print ('...done')

    # Shift the observations and gtis so that they start at bstart0:

    # tobs = tobs - bstart0
    bean.tobs = tobs - bean.tref

    gti_checking = bean.gtiname is not None
    if gti_checking:

        # Read in the gtis (required arrays are st (start time) and et (end time) of times telescope IS observing (indexes need to match)
        # gtis should be in MJD
        # if not exists(gtidata):
        #     sys.exit("** ERROR ** GTI data file {} not found".format(gtidata))
        gtidata = np.loadtxt(bean.gtiname)

        # Now extract and scale gti data:
        bean.st = np.array(gtidata[:,0])-bean.tref
        bean.et = np.array(gtidata[:,1])-bean.tref
        # st = st-bstart0
        # et = et-bstart0

        if bean.burstname is None:
            print(f"""
Observations found:
  persistent fluxes = {bean.pflux}
  gti data file = {bean.gtiname}
GTI checking will be performed""")
        else:
            print(f"""
Observations found:
  burst times = {bean.bstart}
  fluences = {bean.fluen}
  alphas = {bean.alpha}
  persistent fluxes = {bean.pflux}
  gti data file = {bean.gtiname}
GTI checking will be performed""")

        # not sure if you need to return st, et if they're also global variables (defined at start)
        return # bstart0, bstart, fluen, fluene, obs, obs_err, pflux, pfluxe, tobs, st, et

    else:
        if bean.burstname is None:
            print(f"""
Observations found:
  persistent fluxes = {bean.pflux}
You have not supplied a GTI file, so no GTI checking will be performed.""")
        else:
            print(f"""
Observations found:
  burst times = {bean.bstart}
  fluences = {bean.fluen}
  alphas = {bean.alpha}
  persistent fluxes = {bean.pflux}
You have not supplied a GTI file, so no GTI checking will be performed.""")

        bean.st, bean.et = None, None

    return # bstart0, bstart, fluen, fluene, obs, obs_err, pflux, pfluxe, tobs
