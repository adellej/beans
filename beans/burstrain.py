""" This module defines the functions that generate the burst train """

import random

import numpy as np
import matplotlib.pyplot as plt

from .settle import settle

def mean_flux(t1, t2, tobs, a, b):

    # Calculates the mean flux between t1 and t2 from the piecewise linear
    # interpolation of tobs,a,b

    na = len(a)

    if len(tobs) != na + 1 or len(b) != na:
        print("** ERROR ** some problem with tobs, a, b arrays")
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


# ------------------------------------------------------------------------- #

# ------------------------------------------------------------------------- #
# To account for the uncertainty in the persistent flux observations, we re-calculate the fit to the flux evolution each time,
# and vary the persistent flux by a random amount within the uncertainty of the measurement.


def get_a_b(pflux, pfluxe, tobs):

    # Do piecewise continuous fits to the flux evolution, here
    # determine the appropriate parameters for each interval:

    # Now actually calculate the coefficients for the flux fit

    # Linear fit
    ng = len(tobs)

    b0 = np.zeros(ng - 1)
    a0 = np.zeros(ng - 1)

    fpflux = np.zeros(ng)
    for i in range(0, ng):
        fpflux[i] = pflux[i]
    for i in range(1, ng):
        b0[i - 1] = (fpflux[i] - fpflux[i - 1]) / (tobs[i] - tobs[i - 1])
        a0[i - 1] = fpflux[i - 1] - b0[i - 1] * tobs[i - 1]

    a = a0
    b = b0

    result = dict()

    result["a"] = [a]
    result["b"] = [b]

    # No longer used
    # mflux = []
    # for i in range(1, n_burst):
        # mflux.append(mean_flux(bstart[i - 1], bstart[i], tobs, a0, b0))

    # plt.scatter(tobs, fpflux,color='black')
    # plt.scatter(tobs,pflux,color='magenta')
    # plt.scatter(bstart[1:len(bstart)],mflux,color='cyan')
    # plt.show()

    return result


# ------------------------------------------------------------------------- #

# Next we define the different functions that will generate the burst train: generate_burst_train and next_burst

# ------------------------------------------------------------------------- #

run=1
debug=0

def next_burst(
    base,
    z,
    x_0,
    t1,
    tobs,
    a,
    b,
    r1,
    cfac,
    mass,
    radius,
    direction=1,
    debug=False ):
    """
    Routine to find the next burst in the series and return its properties
    Adapted from sim_burst.pro
    """

    mdot_res = 1e-6
    fn = "next_burst"
    assert direction in (1,-1)

    minmdot = 0.0
    maxmdot = 1.0

    # a, b passed as an array of an array
    a = a[0]
    b = b[0]

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
    i0 = max([0, min([len(a) - 1, max(itobs)])])
    mdot0 = ( (0.67 / 8.8) * (a[i0] + b[i0] * t1) * r1 )
    if debug:
        print("{}: z={}, X_0={}, r1={}".format(fn, z, x_0, r1 ))

    # Calculate the burst properties for the trial mdot value
    trial = settle(base, z, x_0, mdot0, cfac, mass, radius)
    if debug:
        print ('{}: initial guess mdot0={} @ t1={}, tdel={}, direction={}'.format(fn,mdot0,t1,trial.tdel,direction))

    # Now update the mdot with the value averaged over the trial interval
    if direction == 1:
        mdot = (0.67 / 8.8) * mean_flux(t1, t1 + trial.tdel / 24.0, tobs, a, b) * r1
    else:
        mdot = (0.67 / 8.8) * mean_flux(t1 - trial.tdel / 24.0, t1, tobs, a, b) * r1

    # Now retain the entire history of this iteration, so we can check for loops

    mdot_hist = [mdot0]
    tdel_hist = [trial.tdel[0]/24.]

    nreturn = 0
    while (abs(mdot - mdot_hist[-1]) > mdot_res / 2.0) \
        and (((t1 + trial.tdel / 24.0 < 2.*max(tobs)) & (direction == 1)) \
            or ((t1 - trial.tdel / 24.0 > min(tobs)-(max(tobs)-min(tobs))) & (direction == -1))) \
        and (mdot > minmdot and mdot < maxmdot):

        trial = settle(base, z, x_0, mdot[0], cfac, mass, radius)
        nreturn = nreturn + 1

        mdot_hist.append(mdot[0])
        tdel_hist.append(trial.tdel[0]/24.)

        if direction == 1:
            mdot = (0.67 / 8.8) * mean_flux(t1, t1 + (trial.tdel / 24.0), tobs, a, b) * r1

        else:
            mdot = (0.67 / 8.8) * mean_flux(t1 - (trial.tdel / 24.0), t1, tobs, a, b) * r1

        # Break out of the loop here, if necessary
        if nreturn > 10:
            e = random.random()
            mdot = mdot_hist[-1] * (1.0 - e) + mdot * e
            # Perhaps you should try to reset this randomly every 10 steps? - dkg
            # Yes, otherwise every trial above 10 steps will be random
            nreturn = 0

    # save the final versions to the history arrays
    mdot_hist.append(mdot[0])
    tdel_hist.append(trial.tdel[0]/24.)
    if debug:
        print('{}: mdot_hist={}'.format(fn, mdot_hist))

        # now produce a diagnostic plot with the debug flag
        # plt.plot(t1+np.array(tdel_hist), mdot_hist, '.', label='tdel history')
        for tdel in tdel_hist:
            plt.axvline(t1+tdel, color='k', ls='--')
        # also calculate a bunch of values to compare with
        t_arr = np.arange(t1, max(tobs), step=0.1)
        m_arr = [0]
        t_arr2 = [t1]
        for t in t_arr[1:]:
            _mdot = (0.67 / 8.8) * mean_flux(t1, t, tobs, a, b) * r1
            _tmp = settle(base, z, x_0, _mdot, cfac, mass, radius)
            t_arr2.append(t1+_tmp.tdel[0]/24.)
            m_arr.append(_mdot)
        plt.plot(t_arr, np.array(t_arr2), '-', label='tdel')
        plt.plot(t_arr, t_arr, '-', label='1:1')
        plt.xlim((0,1.1*max(t1+np.array(tdel_hist))))
        plt.ylim((0,1.1*max(t1+np.array(tdel_hist))))
        # plt.plot(np.array(t_arr2), np.array(m_arr), '.')
        plt.legend()
        plt.show()
        breakpoint()

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
    result.e_b = trial.E_b # multiply eb by 0.8 to account for incomlpete burning of fuel, as in Goodwin et al (2018).
    result.alpha = trial.alpha
    # result.qnuc = tmp.Q_nuc
    # result.xbar = tmp.xbar
    result.mdot = mdot

    return result


# -------------------------------------------------------------------------#

# -------------------------------------------------------------------------#


def generate_burst_train(
    base,
    z,
    x_0,
    r1,
    r2,
    r3,
    mass,
    radius,
    # Now all the parameters below can be passed from a Beans object
    bstart,
    pflux,
    pfluxe,
    tobs,
    numburstssim,
    ref_ind,
    debug=False
):

    # This routine generates a simulated burst train. The output is a
    # structure with the following elements:
    #   BASE            FLOAT          0.175000
    #   Z               FLOAT         0.0100000
    #   X_0             FLOAT          0.440000
    #   R1              FLOAT          0.108533
    #   R2              FLOAT           1.00000
    #   R3              FLOAT           1.00000
    #   RUN             INT              1
    #   DOUBLE          INT              0
    #   FLAG            INT              1
    #   MDOT            DOUBLE    Array[7]
    #   MDOT_MAX        DOUBLE         0.043402292
    #   TIME            DOUBLE    Array[8]
    #   ALPHA           FLOAT     Array[7]
    #   E_B             FLOAT     Array[7]
    #   QNUC            FLOAT     Array[7]
    #   XBAR            FLOAT     Array[7]
    #
    # We have one more element of the time array than the other arrays, because
    # we can't determine the properties for that burst, as we don't have
    # enough data to back project. So the ith element of e_b, alpha etc.
    # belongs with the (i+1)th element of time

    # obs = obs
    # bstart = bstart
    # tobs = tobs
    # a = a
    # b = b

    # if len(double) == 0:
    #    double = 0
    # if len(debug) == 0:
    #    debug = 0
    # if len(run) == 0:
    #    run=1

    mdot_max = -1

    # Now to go ahead and try to simulate the bursts train with the resulting
    # best set of parameters
    # Start from the *second* (observed) burst in the train
    # Originally this would have simulated four bursts following the reference,
    # and three preceding. However, the last burst in the train (the 8th) for
    # runs test17 were wildly variable, so now restrict the extent by one

    if bstart is not None:
        sbt = bstart[ref_ind]
    else:
        # In the absence of any bursts, set the reference time to ref_ind (can be
        # any time within the outburst)
        # sbt = 0.0
        sbt = ref_ind

    salpha = -1
    flag = 1  # Initially OK

    #          for i=0,2*(1+double) do begin
    #  for i=0,3*(1+double) do begin # Do the 5th burst also, forward only

    # salpha = np.zeros(2*(1+double)+1)
    # se_b = np.zeros(2*(1+double)+1)
    # smdot = np.zeros(2*(1+double)+1)
    # sqnuc = np.zeros(2*(1+double)+1)
    # sxbar = np.zeros(2*(1+double)+1)
    # stime=np.zeros(2*(1+double)+1)

    # Get a and b for varying persistent flux:

    # n_burst = len(bstart)
    result0 = get_a_b(pflux, pfluxe, tobs)# , n_burst, bstart)
    a = result0["a"]
    b = result0["b"]

    stime = []  # initialise array to store simulated times
    forward, backward = True, True  # go in both directions at the start
    earliest = sbt  # this is the earliest burst in the train
    latest = sbt    # this is the time of the latest burst in the train
    # for i in range (0,2*(1+double)+1): # Do the 5th burst also, forward only
    for i in range(0, numburstssim):  # Do the 5th burst also, forward only
        # print "i = ",i
        # print 'stime = {}'.format(stime)
        # Now that we reduce the overall number of forward bursts, this if statement
        # is redundant

        #    if i lt 3*(1+double) then $

        # Here we introduce recurrence time corrections since the accretion rate is not flat over the extrapolated time, resulting in the recurrence time being underestimated by settle. Correction factors are from Zac, calculated using KEPLER for SAX J1808.4--3658

        # if i == 0:  # This is observed burst at 1.89
        #     cfac1 = 1.02041
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
            result2 = next_burst(
                base,
                z,
                x_0,
                earliest,
                tobs,
                a,
                b,
                r1,
                1.0,
                mass,
                radius,
                direction=-1,
                debug=debug
            )

        if forward:
            # Also find the time for the *next* burst in the train

            result3 = next_burst(
                base,
                z,
                x_0,
                latest,
                tobs,
                a,
                b,
                r1,
                1.0,
                mass,
                radius,
                direction=1,
                debug=debug
            )

        if result2 is not None:
            # we have a result from the next_burst call going backward, so add its properties to the arrays
            t2 = result2.t2[0]
            _alpha = result2.alpha[0]
            _e_b = result2.e_b[0]
            _mdot = result2.mdot[0]
            if salpha == -1:
                # create the arrays with which to accumulate the results
                stime = [t2, sbt]
                iref = 1    # index for reference burst
                salpha = [_alpha]
                se_b = [_e_b]
                smdot = [_mdot]
            else:
                stime.insert(0, t2)
                iref += 1
                salpha.insert(0, _alpha)
                se_b.insert(0, _e_b)
                smdot.insert(0, _mdot)
            earliest = t2
        else:
            # if the earlier burst has failed, we don't need to pursue any further
            backward = False

        if result3 is not None:
            # we have a result from the next_burst call going forward, so add its properties to the arrays
            t3 = result3.t2[0]
            _alpha2 = result3.alpha[0]
            _e_b2 = result3.e_b[0]
            _mdot2 = result3.mdot[0]
            if salpha == -1:
                # This shouldn't happen, as we should be able to get at least one earlier burst
                stime = [sbt, t3]
                iref = 0
                salpha = [_alpha2]
                se_b = [_e_b2]
                smdot = [_mdot2]
            else:
                salpha.append(_alpha2)
                se_b.append(_e_b2)
                smdot.append(_mdot2)
                stime.append(t3)
            latest = t3

        # Check the results here

        # I don't think t2 or t3 are ever set to these "dummy" values anymore
        # if abs(t2) == 99.99 or abs(t3) == 99.99:
        if not (forward or backward):
            break

    if (mdot_max == -1) & (len(stime) > 0):

        mdot_max = max(smdot)

    result = dict()

    result["base"] = [base]
    result["z"] = [z]
    result["x_0"] = [x_0]
    result["r1"] = [r1]
    result["r2"] = [r2]
    result["r3"] = [r3]
    result["time"] = stime
    result["mdot_max"] = [mdot_max]
    if len(stime) > 0:
        # The simulation might fail to generate any bursts, so only add the arrays if they exist
        result["mdot"] = smdot
        result["iref"] = iref
        result["alpha"] = salpha
        result["e_b"] = se_b
        #print(f"In burstrain fluence is {se_b}")

    # result["qnuc"] = sqnuc
    # result["xbar"] = sxbar
    result["mass"] = [mass]
    result["radius"] =  [radius]
    result["forward"] = forward     # to keep track of the outcome of each direction
    result["backward"] = backward

    return result


# ------------------------------------------------------------------------- #
