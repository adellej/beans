# -------------------------------------------------------------------------#
import numpy as np
from burstrain import *

def runmodel(theta_in, y, tref, bstart, pflux, pfluxe, tobs, numburstssim, numburstsobs, ref_ind, gti_checking,train,
             gti_start=None, gti_end=None, debug=False):
    """
    This routine calls one of two functions that generate the burst model
    predictions, either generate_burst_train or burstensemble, depending on
    which mode the analysis is in

    :param theta_in: parameter tuple: X, Z, Q_b, f_a, f_E, r1, r2, r3,
      mass & radius
    :param y:
    :param tref:
    :param bstart:
    :param pflux:
    :param pfluxe:
    :param tobs:
    :param numburstssim:
    :param numburstsobs:
    :param ref_ind:
    :param gti_checking:
    :param rain:
    :param gti_start:
    :param gti_end:
    :param debug: set to True to display more diagnostic information
    :return: model array with predicted burst parameters, Boolean giving
      the validity of the solution (i.e. consistent with the GTI information),
      and the full result dict from generate_burst_train/burstensemble
    """

    if debug:
        print('Calling runmodel')

    X, Z, Q_b, f_a, f_E, r1, r2, r3, mass, radius = theta_in
    #    X, Z, Q_b, s_t, f_a, f_E, r1, r2, r3 = theta

    if train:
        # Now call the function. From the code:
        # This routine generates a simulated burst train. The output is a
        # dict with the following keys:
        # ['base', 'z', 'x_0', 'r1', 'r2', 'r3', 'time', 'mdot_max', 'mdot',
        #  'iref', 'alpha', 'e_b', 'mass', 'radius', 'forward', 'backward']

        result = generate_burst_train(
            Q_b, Z, X, r1, r2, r3, mass, radius, bstart, pflux, pfluxe, tobs, numburstssim, ref_ind
        )

        tpred = result["time"]

        # bug testing sample model values
        # model = np.array([-0.04106, 2.77858, 3.99188, 3.73303, 3.68542, 4.16907, 4.71480, 115.000, 126.903, 138.070]) #to test the code we define the model as an array, where the numbers in the array are values for the parameters of y, in the same order as y. Values have been taken from Duncan's example code output

        if (y is not None) & (len(tpred) > 0):

            # Assemble the array for comparison with the data
            # First the burst times; we dynamically determine which
            # predicted burst is closest in time to each observed one, and
            # assign that predicted burst to it for comparison of the
            # properties

            i1 = [] # index mapping predicted bursts to observed 
            for i in range(0, ref_ind):
                i1.append(np.argmin(np.abs(tpred - y[i])))

            i1.append(np.argmin(np.abs(tpred - tref)))

            for i in range(ref_ind, len(bstart) - 1):
                i1.append(np.argmin(np.abs(tpred - y[i])))

            # We compare the fluences for all the bursts
            # We subtract 1 here, and also to the expression for ialpha, because the indexing is different for
            # the e_b and alpha arrays, compared to the times, coming out of generate_burst_train

            ie_b = [np.max([x - 1, 0]) for x in i1]

            # We only compare the times of the bursts for observed events #0, 2 & 3; #10 is a "reference"
            # from which the train is calculated

            i2 = i1
            i2.pop(ref_ind)
            itime = i2

            # We compare the alpha values only for observed events #1, 2 & 3 as we don't have the recurrence
            # time for event #0

            i11 = []
            for i in range(0, ref_ind):
                i11.append(np.argmin(np.abs(tpred - y[i])))

            i11.append(np.argmin(np.abs(tpred - tref)))

            for i in range(ref_ind, len(bstart)):
                i11.append(np.argmin(np.abs(tpred - y[i])))
            li11 = list(i11)
            li1m11 = [np.max([x - 1, 0]) for x in li11]

            i3 = li1m11

            i3.pop(0)
            ialpha = i3

            model = []
            for i in range(0, len(bstart) - 1):
                model.append(result["time"][itime[i]])
            for i in range(0, len(bstart)):
                model.append(result["e_b"][ie_b[i]])
            for i in range(0, len(bstart) - 1):
                model.append(result["alpha"][ialpha[i]])

            model = np.array(model)
        else:
	    # If you're not comparing to observed bursts, just return the
	    # result of generate_burst_train
            # This loop will (also?) be triggered if the call to 
            # generate_burst_train results in no bursts. That can happen if
            # the model parameters are nonsensical - e.g. Z<0. An "unhashable
            # type" error will then be triggered - dkg
            # TODO: fix handling of empty tpred array in runmodel
            model = result
            i1 = [] # created as empty list for the gti_checking below

    else:
        # If we're not generating a burst train, just run the ensemble

        result = burstensemble(
            base, x, z, r1,r2,r3,mass,radius,bstart,pflux,numburstsobs)

        model = np.concatenate((result['time'], result['e_b'], result['alpha']))

    # Check here if the model instance is valid, i.e. the bursts that are NOT
    # matched with the observed ones must fall in gaps
    # Originally we used the (global) arrays st, et defined by 1808-match, to
    # avoid copying them over from IDL each time; but now these are passed
    # as parameters gti_start, gti_end
    valid = True
    if gti_checking == 1:
        # if "st" not in globals():
        if (gti_start is None) or (gti_end is None):
            print ('** WARNING ** can''t access GTI information')
            return model, valid, result
        else:
            st, et = gti_start, gti_end

        for index, rt in enumerate(tpred):
            if index not in i1:
                # ok, not one of the known bursts. Is it an excluded time?#
                for i in range(len(st)):

                    if rt >= st[i] and rt <= et[i] - 10.0 / 86400.0:

                        valid = False
                        return model, valid, result

    # Check here if anisoptropy estimates are consistent with Fujimoto model
    #   sqrt = (r1*r2*r3*1e3)/(63.23*0.74816)
    #   xi_p = np.power(sqrt, 0.5)
    #   xi_b = (0.74816*xi_p)/r2
    #   if not (((1./xi_b)-0.5)+0.2) > 1/(2*xi_p) > (((1./xi_b)-0.5)-0.2):
    #       valid = False
    #       print('anisotropies not consistent')
    #       return model, valid
    #   else:
    #       valid = True

    if debug:
        print(f'model = {model}')

    return model, valid, result


# -------------------------------------------------------------------------#
