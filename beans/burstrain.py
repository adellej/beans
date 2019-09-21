""" This module defines the functions that generate the burst train """

import numpy as np
from settle import settle
import random

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


def get_a_b(pflux, pfluxe, tobs, n_burst, bstart):

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

    mflux = []
    for i in range(1, n_burst):
        mflux.append(mean_flux(bstart[i - 1], bstart[i], tobs, a0, b0))

    # plt.scatter(tobs, fpflux,color='black')
    # plt.scatter(tobs,pflux,color='magenta')
    # plt.scatter(bstart[1:len(bstart)],mflux,color='cyan')
    # plt.show()

    return result


# ------------------------------------------------------------------------- #

# Next we define the different functions that will generate the burst train: generate_burst_train and next_burst

# ------------------------------------------------------------------------- #
debug = 0
double = 0
run = 1


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
    direction,
    run=run,
    debug=debug,
):
    # t2,e_b,alpha,mdot,qnuc,xbar,

    # Routine to find the next burst in the series and return it's properties
    # Adapted from sim_burst.pro

    # common cube
    # if len(direction) == 0 then direction=1.
    # if len(dbg) == 0 then dbg=0
    # if len(run) == 0 then run=0
    mdot_res = 0.001
    fn = "next_burst"
    tol = 1.0e-6
    seed = 49221941
    # if fine == 0:
    #    mdot_res==0.001
    # if fine == 1:
    #    mdot_res==0.0004
    #  if fine eq 2 then mdot_res=0.00025
    # if fine == 2 or run == 1:
    # mdot_res==0.0002 # Latest newcube4/5

    minmdot = 0.0
    maxmdot = 1.0

    # change precision of mdot_res:
    # mdot_res = mdot_res/10.
    a = a[0]
    b = b[0]

    #  i0=min([n_elements(a)-1,max(where(t1 gt tobs))])
    itobs = max([np.where(t1 > tobs)])
    itobs = max(itobs)
    if len(itobs) == 0:
        itobs = [-1]
    # i0=max([0,min([len(a)-1,max([i for i, value in enumerate(tobs) if value < t1])])])
    i0 = max([0, min([len(a) - 1, max(itobs)])])
    mdot0 = (
        (0.67 / 8.8) * (a[i0] + b[i0] * t1) * r1
    )  # Initial guess at mean mdot (linear)

    # print,base,x_0,z,mdot0
    tmp = settle(base, z, x_0, mdot0, cfac, mass, radius)
    if direction == 1:
        mdot = (0.67 / 8.8) * mean_flux(t1, t1 + tmp.tdel / 24.0, tobs, a, b) * r1
    else:
        mdot = (0.67 / 8.8) * mean_flux(t1 - tmp.tdel / 24.0, t1, tobs, a, b) * r1

    #  while abs(mdot-mdot0)/mdot gt 0.01 do begin	; To 1%

    # Now retain the entire history of this iteration, so we can check for loops

    mdot_hist = [mdot0]

    # test = np.where(abs(mdot_hist-mdot) < tol)
    #  while abs(mdot-mdot0) gt mdot_res/2. $
    # while abs(mdot-mdot_hist[len(mdot_hist)-1]) > mdot_res/2. and (mdot > minmdot and mdot < maxmdot):

    #mdot = np.float64(mdot)

    nreturn = 0

    # while ((abs(mdot-mdot_hist) > mdot_res/2.) - (mdot > minmdot and mdot < maxmdot)).any():
    while abs(mdot - mdot_hist[len(mdot_hist) - 1]) > mdot_res / 2.0 and (
        mdot > minmdot and mdot < maxmdot
    ):
        #    mdotm1=mdot0
        #    mdot0=mdot

        mdot_hist.append(mdot)
        tmp = settle(base, z, x_0, mdot[0], cfac, mass, radius)

        if direction == 1:
            mdot = (0.67 / 8.8) * mean_flux(t1, t1 + (tmp.tdel / 24.0), tobs, a, b) * r1
            nreturn = nreturn + 1

        else:
            mdot = (0.67 / 8.8) * mean_flux(t1 - (tmp.tdel / 24.0), t1, tobs, a, b) * r1
            nreturn = nreturn + 1

        # Break out of the loop here, if necessary
        # print 'nreturn = {}'.format(nreturn)
        # test = np.where(abs(mdot_hist-mdot) < tol)
        #    if abs(mdotm1-mdot) lt tol then begin
        if nreturn > 10:
            #      mdot=mdot0+(mdot-mdot0)*randomu(seed)
            e = random.random()
            #      mdot=mdot0+(mdot-mdot0)*randomu(seed)
            mdot = mdot_hist[len(mdot_hist) - 1] * (1.0 - e) + mdot * e
        # if dbg:
        # print 'next_burst: randomizing mdot...'

    # print,abs(mdot-mdot0)/mdot
    # t2 = 0.0
    # if dbg:
    #    print mdot0,tmp

    if mdot < minmdot or mdot > maxmdot and run == 0:
        if mdot < minmdot:
            t2 = -99.99
        else:
            t2 = +99.99

        # create array
    result = np.recarray(
        (1,), dtype=[("t2", np.float64), ("e_b", np.float64), ("alpha", np.float64)]
    )
    # assign elements
    result.t2 = t1 + direction * tmp.tdel / 24.0
    result.e_b = (
        tmp.E_b
    )  # multiply eb by 0.8 to account for incomlpete burning of fuel, as in Goodwin et al (2018).
    result.alpha = tmp.alpha
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
    bstart,
    pflux,
    pfluxe,
    tobs,
    numburstssim,
    run=run,
    double=double,
    debug=debug,
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

    sbt = bstart[1]
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

    n_burst = len(bstart)
    result0 = get_a_b(pflux, pfluxe, tobs, n_burst, bstart)
    a = result0["a"]
    b = result0["b"]

    stime = []
    # for i in range (0,2*(1+double)+1): # Do the 5th burst also, forward only
    for i in range(0, numburstssim):  # Do the 5th burst also, forward only
        # print "i = ",i
        # print 'stime = {}'.format(stime)
        # Now that we reduce the overall number of forward bursts, this if statement
        # is redundant

        #    if i lt 3*(1+double) then $

        # Here we introduce recurrence time corrections since the accretion rate is not flat over the extrapolated time, resulting in the recurrence time being underestimated by settle. Correction factors are from Zac, calculated using KEPLER

        if i == 0:  # This is observed burst at 1.89
            cfac1 = 1.02041
            cfac2 = 1.02041
        if (
            i == 1
        ):  # to the right this is 3rd observed burst, to left it is predicted burst
            cfac1 = 1.00
            cfac2 = 1.1905
        if (
            i == 2
        ):  # to the right this is 4th observed burst, to left is predicted burst
            cfac1 = 1.00
            cfac2 = 1.2346
        if (
            i == 3
        ):  # to the right this is final predicted burst, to the left is first observed burst (note that cfac = 1.25 is estimated interpolation)
            cfac1 = 1.00
            cfac2 = 1.25
        # if i == 4: # to the right this is final predicted burst, to the left is first observed burst (note that cfac = 1.25 is estimated interpolation)
        #    cfac1 = 0.98
        #    cfac2 = 1.27

        if len(stime) == 0:

            # Find the time for the *previous* burst in the train, t2

            result2 = next_burst(
                base,
                z,
                x_0,
                sbt,
                tobs,
                a,
                b,
                r1,
                cfac1,
                mass,
                radius,
                direction=-1,
                run=run,
                debug=debug,
            )

            # Also find the time for the *next* burst in the train, t3

            result3 = next_burst(
                base,
                z,
                x_0,
                sbt,
                tobs,
                a,
                b,
                r1,
                cfac2,
                mass,
                radius,
                direction=1,
                run=run,
                debug=debug,
            )

        else:

            # Find the time for the *previous* burst in the train, t2
            result2 = next_burst(
                base,
                z,
                x_0,
                stime[0],
                tobs,
                a,
                b,
                r1,
                cfac1,
                mass,
                radius,
                direction=-1,
                run=run,
                debug=debug,
            )

            # Also find the time for the *next* burst in the train, t3

            result3 = next_burst(
                base,
                z,
                x_0,
                stime[len(stime) - 1],
                tobs,
                a,
                b,
                r1,
                cfac2,
                mass,
                radius,
                direction=1,
                run=run,
                debug=debug,
            )

        t2 = result2.t2
        t2 = t2[0]
        t3 = result3.t2
        t3 = t3[0]
        _alpha = result2.alpha
        _alpha = _alpha[0]
        _alpha2 = result3.alpha
        _alpha2 = _alpha2[0]
        _e_b = result2.e_b
        _e_b = _e_b[0]
        _e_b2 = result3.e_b
        _e_b2 = _e_b2[0]
        # _qnuc = result2.qnuc
        # _qnuc = _qnuc[0]
        # _qnuc2 = result3.qnuc
        # _qnuc2 = _qnuc2[0]
        # _xbar = result2.xbar
        # _xbar = _xbar[0]
        # _xbar2 = result3.xbar
        # _xbar2 = _xbar2[0]
        _mdot = result2.mdot
        _mdot = _mdot[0]
        _mdot2 = result3.mdot
        _mdot2 = _mdot2[0]

        # Check the results here

        if abs(t2) == 99.99 or abs(t3) == 99.99:
            flag = 0
            i = 3 * (1 + double)  # replace this with a break statement?
        if t2 == 99.99 or t3 == 99.99:
            mdot_max = 99.99
        if t2 == -99.99 or t3 == -99.99:
            mdot_max = -99.99
        else:

            # We have some good results, so append to the sbt array (and others)

            if salpha == -1:
                stime = [t2, sbt, t3]
                salpha = [_alpha, _alpha2]
                se_b = [_e_b, _e_b2]
                smdot = [_mdot, _mdot2]
                # sqnuc = [_qnuc,_qnuc2]
                # sxbar = [_xbar,_xbar2]

            else:
                salpha.insert(0, _alpha)
                salpha.append(_alpha2)

                se_b.insert(0, _e_b)
                se_b.append(_e_b2)

                smdot.insert(0, _mdot)
                smdot.append(_mdot2)
                #
                # sqnuc.insert(0,_qnuc)
                # sqnuc.append(_qnuc2)
                #
                # sxbar.insert(0,_xbar)
                # sxbar.append(_xbar2)

                stime.insert(0, t2)
                stime.append(t3)

    if mdot_max == -1:

        mdot_max = max(smdot)

    result = dict()

    result["base"] = [base]
    result["z"] = [z]
    result["x_0"] = [x_0]
    result["r1"] = [r1]
    result["r2"] = [r2]
    result["r3"] = [r3]
    result["mdot"] = smdot
    result["mdot_max"] = [mdot_max]
    result["time"] = stime
    result["alpha"] = salpha
    result["e_b"] = se_b
    # result["qnuc"] = sqnuc
    # result["xbar"] = sxbar
    result["mass"] = [mass]
    result["radius"] = [radius]

    return result


# ------------------------------------------------------------------------- #
