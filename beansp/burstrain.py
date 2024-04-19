""" This module defines the functions that generate the burst train """

import numpy as np
import random
import matplotlib.pyplot as plt

# load local  modules
from .settle import settle

# Next we define the different functions that will generate the burst train: generate_burst_train and next_burst

# ------------------------------------------------------------------------- #

run=1
debug=0

def next_burst( bean, base, x_0, z, t1, dist, xi_p, cfac, mass, radius,
    direction=1, debug=False ):
    """
    Routine to find the next burst in the series and return its properties
    Adapted from sim_burst.pro

    :param bean: Beans object, from which the other inputs are read in:
      tobs, a, b,
    :param base: base flux [MeV/nucleon]
    :param x_0: accreted H-fraction
    :param z: accreted CNO metallicity
    :param t1: time to start burst search at
    :param dist: source distance (kpc)
    :param xi_p: anisotropy of persistent emission
    :param cfac: scale factor for recurrence time, fluence
    :param mass: NS mass (M_sun)
    :param radius: NS radius (km)
    :param direction: forward (+1) or backward (-1) in time
    :param debug: set to True to show additional debugging information
    :return:
    """

    tobs = bean.tobs

    debug_plot=False # for now don't do the debug plots 

    mdot_res = 1e-6
    fn = "next_burst"
    assert direction in (1,-1)

    minmdot = 0.0
    maxmdot = 1.0

    # Determine the initial guess for mean mdot (linear)
    #  i0=min([n_elements(a)-1,max(where(t1 gt tobs))])
    itobs = np.where(t1 > tobs)[0]
    if (len(itobs) == 0) & (direction == -1):
        # the start time is before *any* of the observations; don't bother!
        return None
    if len(itobs) == 0:
        # this makes no sense to me; if the t1 value is < all the tobs values, then the
        # nearest element would be the zeroth
        # itobs = [-1]
        itobs = [0]
    # i0=max([0,min([len(a)-1,max([i for i, value in enumerate(tobs) if value < t1])])])
    # Now that we have a couple of options for interpolation, need to
    # remove the reliance on the linear interpolation parameters
    # i0 = max([0, min([len(a) - 1, max(itobs)])])
    # mdot0 = (0.67 / 8.8) * bean.pflux[itobs[0]] * r1 
    mdot0 = bean.flux_to_mdot(x_0, dist, xi_p, mass, radius, bean.pflux[itobs[-1]])
    if debug:
        print("{}: z={}, X_0={}, mdot_0={}, itobs={}".format(
            fn, z, x_0, mdot0, itobs[-1] ))

    # Calculate the burst properties for the trial mdot value
    trial = settle(base, z, x_0, mdot0, mass, radius, corr=bean.corr)

    if debug:
        print ('{}: initial guess mdot0={} @ t1={}, tdel={}, direction={}'.format(fn,mdot0,t1,trial.tdel,direction))

    # Now update the mdot with the value averaged over the trial interval
    # not sure why trial.tdel has suddenly become a vector
    if direction == 1:
        mdot = bean.flux_to_mdot(x_0, dist, xi_p, mass, radius,
            bean.mean_flux(t1, t1 + trial.tdel[0] / 24.0, bean) )
    else:
        mdot = bean.flux_to_mdot(x_0, dist, xi_p, mass, radius,
            bean.mean_flux(t1 - trial.tdel[0] / 24.0, t1, bean) )

    # Now retain the entire history of this iteration, so we can check for loops
    mdot_hist = [mdot0]
    tdel_hist = [trial.tdel[0]/24.]

    nreturn = 0
    nreturn_total = 0
    while (abs(mdot - mdot_hist[-1]) > mdot_res / 2.0) \
        and (((t1 + trial.tdel / 24.0 < 2.*max(tobs)) & (direction == 1)) \
            or ((t1 - trial.tdel / 24.0 > min(tobs)-(max(tobs)-min(tobs))) & (direction == -1))) \
        and (mdot > minmdot and mdot < maxmdot):

        trial = settle(base, z, x_0, mdot, mass, radius, corr=bean.corr)
        nreturn = nreturn + 1
        nreturn_total = nreturn_total + 1

        mdot_hist.append(mdot)
        tdel_hist.append(trial.tdel[0]/24.)

        if direction == 1:
            mdot = bean.flux_to_mdot(x_0, dist, xi_p, mass, radius, 
                bean.mean_flux(t1, t1 + (trial.tdel[0] / 24.0), bean) )

        else:
            mdot = bean.flux_to_mdot(x_0, dist, xi_p, mass, radius, 
                bean.mean_flux(t1 - (trial.tdel[0] / 24.0), t1, bean) )

        # Break out of the loop here, if necessary
        if nreturn > 10:
            e = random.random()
            mdot = mdot_hist[-1] * (1.0 - e) + mdot * e
            # Perhaps you should try to reset this randomly every 10 steps? - dkg
            # Yes, otherwise every trial above 10 steps will be random
            nreturn = 0

        if nreturn_total > 30:
            e = random.random()
            mdot = mdot_hist[-1] * (1.0 - e)+ mdot * e


    # save the final versions to the history arrays
    mdot_hist.append(mdot)
    tdel_hist.append(trial.tdel[0]/24.)
    if debug:
        print('{}: nreturn={}, nreturn_total={}, mdot_hist={}'.format(
            fn, nreturn, nreturn_total, mdot_hist))

        if debug_plot:
            # now produce a diagnostic plot with the debug flag
            # plt.plot(t1+np.array(tdel_hist), mdot_hist, '.', label='tdel history')
            for tdel in tdel_hist:
                plt.axvline(t1+tdel, color='k', ls='--')
            # also calculate a bunch of values to compare with
            t_arr = np.arange(t1, max(tobs), step=0.1)
            m_arr = [0]
            t_arr2 = [t1]
            for t in t_arr[1:]:
                _mdot = bean.flux_to_mdot(x_0, dist, xi_p, mass, radius,
                    bean.mean_flux(t1, t, bean) )
                _tmp = settle(base, z, x_0, _mdot, mass, radius, corr=bean.corr)
                t_arr2.append(t1+_tmp.tdel[0]/24.)
                m_arr.append(_mdot)
            plt.plot(t_arr, np.array(t_arr2), '-', label='tdel')
            plt.plot(t_arr, t_arr, '-', label='1:1')
            plt.xlim((0,1.1*max(t1+np.array(tdel_hist))))
            plt.ylim((0,1.1*max(t1+np.array(tdel_hist))))
            # plt.plot(np.array(t_arr2), np.array(m_arr), '.')
            plt.legend()
            plt.show()

    # if mdot < minmdot or mdot > maxmdot:
    if abs(mdot - mdot_hist[-2]) > mdot_res / 2.0:
        return None

        # create array
    #print(f'{fn}: mdot={mdot}, tdel={trial.tdel}')
    result = np.recarray(
        (1,), dtype=[("t2", np.float64), ("e_b", np.float64), ("alpha", np.float64)]
    )
    # assign elements
    result.t2 = t1 + direction * trial.tdel / 24.0
    # at one point we were multiplying eb by 0.8 to account for incomlpete
    # burning of fuel, as in Goodwin et al (2018).
    result.e_b = trial.E_b
    result.alpha = trial.alpha
    # result.qnuc = tmp.Q_nuc
    # result.xbar = tmp.xbar
    result.mdot = mdot

    return result


# -------------------------------------------------------------------------#

# -------------------------------------------------------------------------#


def generate_burst_train( bean, base, x_0, z, dist, xi_p, mass, radius,
    full_model=False, debug=False):
    """
    This routine generates a simulated burst train based on the model
    input parameters, and the mdot history inferred from the persistent
    flux measurements (tobs, pflux, pfluxe)

    :param bean: Beans object, from which the remaining parameters are drawn:
      bstart, pflux, pfluxe, tobs, numburstssim, ref_ind
    :param base: base flux [MeV/nucleon]
    :param x_0: accreted H-fraction
    :param z: accreted CNO metallicity
    :param dist: source distance (kpc)
    :param xi_p: anisotropy of persistent emission
    :param mass: NS mass (M_sun)
    :param radius: NS radius (km)
    :param full_model: if set to True, include all the parameters in the
      dict that is returned
    :param debug: set to True to show additional debugging information

    :return: a dictionary with the following keys: ``['time', 'mdot', 'alpha', 'e_b', 'iref']``, and optionally (if full_model is set to True) also: ``['base', 'z', 'x_0', 'dist', 'xi_p', 'mdot_max', 'mass', 'radius', 'forward', 'backward']`` (all but the last 2 are copies of the input parameters).

    The main outputs are the ``time``, ``mdot``, ``alpha``, and ``e_b`` elements
    (numpy arrays) which are the predicted properties of the bursts and
    the accretion rate inferred for each interval.

    We have one more element of the time array than the other arrays, because
    we can't determine the properties for that burst, as we may not have
    enough data to back project. So the ith element of e_b, alpha etc.
    belongs with the (i+1)th element of time:

    .. code-block:: text

        time  [ 0  1  2  3  4  5  6  7  ...  n ]
        mdot     [ 0  1  2  3  4  5  6  ... n-1 ]
        alpha    [ 0  1  2  3  4  5  6  ... n-1 ]
        e_b      [ 0  1  2  3  4  5  6  ... n-1 ]

    """

    fn = "generate_burst_train"
    forward, backward = True, True  # go in both directions at the start

    mdot_max = -1

    cfac = 1.0  # default to take the output of settle directly
    # cfac = 0.65 # temporary to try to match the objective of Goodwin et
                  #   al. (2019)

    # Now to go ahead and try to simulate the bursts train with the resulting
    # best set of parameters
    # Start from the *second* (observed) burst in the train
    # Originally this would have simulated four bursts following the reference,
    # and three preceding. However, the last burst in the train (the 8th) for
    # runs test17 were wildly variable, so now restrict the extent by one

    if bean.bstart is not None:
        sbt = bean.bstart[bean.ref_ind]
    else:
	# In the absence of any bursts, set the reference time to ref_ind
	# (can be any time within the outburst)
	# TODO check if this is consistent with the usage for initial
	# value of earliest/latest below; not sure of the use case
        # sbt = 0.0
        sbt = bean.ref_ind

    flag = 1  # Initially OK

    stime, iref = [sbt], 0  # initialise list to store simulated times
    salpha, se_b, smdot = [], [], [] # ditto for other lists
    earliest = sbt  # (running) time of the earliest burst in the train
    latest = sbt    # (running) time of the latest burst in the train

    for i in range(0, bean.numburstssim):

        assert earliest == stime[0] # I think these variables are redundant
        assert latest == stime[-1]
        if debug:
            print ("{}: simulating burst {} of {}".format(fn, i, bean.numburstssim))

        # Here we adopted recurrence time corrections for SAX
	# J1808.4--3658 ,since the accretion rate is not constant over the
	# extrapolated time, resulting in the recurrence time being
	# underestimated by settle. Correction factors are from Zac
	# Johnston, calculated using KEPLER

	# if i == 0:  # This is observed burst at 1.89 cfac1 = 1.02041
        #     cfac2 = 1.02041
        # if (
        #     i == 1
        # ):  # to the right this is 3rd observed burst, to left it is predicted burst
        #     cfac1 = 1.00
        #     cfac2 = 1.1905
        # if (
        #     i == 2
        # ):  # to the right this is 4th observed burst, to left is predicted burst
        #     cfac1 = 1.00
        #     cfac2 = 1.2346
        # if (
        #     i == 3
        # ):  # to the right this is final predicted burst, to the left is first observed burst (note that cfac = 1.25 is estimated interpolation)
        #     cfac1 = 1.00
        #     cfac2 = 1.25
        # if i == 4: # to the right this is final predicted burst, to the left is first observed burst (note that cfac = 1.25 is estimated interpolation)
        #    cfac1 = 0.98
        #    cfac2 = 1.27

        if backward:
            # Find the time for the *previous* burst in the train
            result_b = next_burst( bean, base, x_0, z, earliest, 
                dist, xi_p, cfac, mass, radius, direction=-1, debug=debug)

            if result_b is not None:
                # we have a result from the next_burst call going backward, so add its properties to the arrays
                stime.insert(0, result_b.t2[0])
                iref += 1
                salpha.insert(0, result_b.alpha[0])
                se_b.insert(0, result_b.e_b[0])
                smdot.insert(0, result_b.mdot)
                earliest = result_b.t2[0]
            else:
                # if the earlier burst has failed, we don't need to pursue any further
                if debug:
                    print("{}: abandoning backward search, step {}".format(fn, i))
                backward = False

        if forward:
            # Also find the time for the *next* burst in the train
            result_f = next_burst( bean, base, x_0, z, latest, 
                dist, xi_p, cfac, mass, radius, direction=1, debug=debug)

            if result_f is not None:
                # we have a result from the next_burst call going forward, so add its properties to the arrays
                stime.append(result_f.t2[0])
                salpha.append(result_f.alpha[0])
                se_b.append(result_f.e_b[0])
                smdot.append(result_f.mdot)
                latest = result_f.t2[0]
            else:
    
                if debug:
                    print("{}: abandoning backward search, step {}".format(fn, i))
                forward = False

        # Check the results here

        # I don't think t2 or t3 are ever set to these "dummy" values anymore
        # if abs(t2) == 99.99 or abs(t3) == 99.99:
        if not (forward or backward):
            break

    if (mdot_max == -1) & (len(stime) > 1):

        # only run this if len(stime) > 1, because it is initialised as a
        # list with a single value (sbt)
        mdot_max = max(smdot)

    result = dict()

    if full_model:
        # model parameters are redundant for the model returned
        result["base"] = [base]
        result["z"] = [z]
        result["x_0"] = [x_0]
        result["dist"] = [dist]
        result["xi_p"] = [xi_p]

        result["mdot_max"] = [mdot_max]

        result["mass"] = [mass]
        result["radius"] =  [radius]

        result["forward"] = forward     # to keep track of the outcome of each direction
        result["backward"] = backward

    # now the actual predictions

    result["time"] = np.array(stime)
    if len(stime) > 0:
        # The simulation might fail to generate any bursts, so only add the arrays if they exist
        result["mdot"] = np.array(smdot)
        result["iref"] = iref
        result["alpha"] = np.array(salpha)
        result["e_b"] = np.array(se_b)
        #print(f"In burstrain fluence is {se_b}")

    if debug:
        print ("{}: train complete, result={}".format(fn, result))

    return result


def burstensemble( bean, base, x_0, z, dist, xi_p, mass, radius, full_model=False ):
    """
    This routine generates as many burst predictions as there are burst
    measurements.
    Written initially by Luke Waterson, 2021

    :param bean: Beans object, from which the remaining parameters are drawn:
      bstart, pflux, numburstsobs
    :param base: base flux [MeV/nucleon]
    :param z: accreted CNO metallicity
    :param x_0: accreted H-fraction
    :param dist: source distance (kpc)
    :param xi_p: anisotropy of persistent emission
    :param mass: NS mass (M_sun)
    :param radius: NS radius (km)

    :return: a dictionary with the following keys: ``['base', 'z', 'x_0', 'dist', 'xi_p', 'time', 'mdot_max', 'mdot', 'iref', 'alpha', 'e_b', 'mass', 'radius']``

    The main outputs (as for :meth:`burstrain.generate_burst_train`) are
    the ``time``, ``e_b``, ``alpha``, and ``mdot``  elements (numpy
    arrays) which are the predicted properties of the bursts and the
    accretion rate inferred for each interval.
    """

    salpha = []
    stime = []
    smdot = []
    se_b = []

    mdot = bean.flux_to_mdot(x_0, dist, xi_p, mass, radius, bean.pflux)

    for i in range(0, bean.numburstsobs):

        tmp = settle(base, z, x_0, mdot[i], mass, radius, corr=bean.corr)

        # accumulate the predictions into the arrays here

        se_b.append(tmp.E_b[0])
        salpha.append(tmp.alpha[0])
        smdot.append(mdot[i])
        # stime.append(bstart[i])
        stime.append(tmp.tdel[0])

    mdot_max = max(smdot)

    result = dict()

    if full_model:
        # model parameters are redundant for the model returned
        result["base"] = [base]
        result["z"] = [z]
        result["x_0"] = [x_0]
        result["dist"] = [dist]
        result["xi_p"] = [xi_p]

        result["mdot_max"] = [mdot_max]

        result["mass"] = [mass]
        result["radius"] = [radius]

    # now the actual predictions

    result["time"] = np.array(stime)
    result["mdot"] = np.array(smdot)
    result["alpha"] = np.array(salpha)
    result["e_b"] = np.array(se_b)

    # omit the printing for now, as it prevents assessing the progress
    # print('ensemble')
    # print(f"In burstrain fluence is {se_b}")

    return result
