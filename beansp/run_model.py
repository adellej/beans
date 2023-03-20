# -------------------------------------------------------------------------#
import numpy as np

# load local module
from .burstrain import *

def runmodel(theta_in, y, tref, bstart, pflux, pfluxe, tobs, numburstssim, numburstsobs, ref_ind, gti_checking,train,
             gti_start=None, gti_end=None, debug=False):
    """
    This routine calls one of two functions that generate the burst model
    predictions, either generate_burst_train or burstensemble, depending on
    which mode the analysis is in
    It then assembles the predictions from the model into a form that can
    be compared to the observations (with appropriate scaling)

    :param theta_in: parameter tuple: X, Z, Q_b, f_a, f_E, r1, r2, r3,
      mass & radius
    :param y: observation array, comprising burst times, fluences and alpha
      values
    :param tref: reference time for bursts and observations, generally the
      first observed burst
    :param bstart: start times of observed bursts (MJD), relative to tref
    :param pflux: persistent flux, used to infer the accretion rate
    :param pfluxe: uncertainty on persistent flux
    :param tobs: times of observations for flux measurements (MJD), relative
      to tref
    :param numburstssim: number of bursts to simulate, for the "train" mode,
      both earlier and later than the reference burst
    :param numburstsobs: number of bursts observed (probably redundant)
    :param ref_ind: index (pointer to bstart) of reference burst
    :param gti_checking: flag for GTI checking (probably redundant)
    :param train: run mode, set to True for "train" (sequential) mode, False
      otherwise
    :param gti_start: array of start times for GTIs (good time intervals)
    :param gti_end: array of end times for GTIs
    :param debug: set to True to display more diagnostic information
    :return: model array with predicted burst parameters, Boolean giving
      the validity of the solution (i.e. consistent with the GTI information),
      and the full result dict from generate_burst_train/burstensemble
    """

    if debug:
        print('Calling runmodel')

    def next_nearest(target, trials, exclude):
        """
        Function to find the element of an array that is closest (in absolute
        terms) to a particular value, excluding some "masked" elements
        Used to match the predicted bursts with those observed

        :param target: the target value we're trying to match
        :param trials: the set of candidates
        :param exclude: a list of indices that are excluded
        :return: index of matching item, or None if all are excluded
        """

        isort = np.argsort(np.abs(trials-target))
        i = 0
        while (isort[i] in exclude) and (i < len(isort)-1):
            i += 1

        if isort[i] in exclude:
            return None

        return isort[i]


    X, Z, Q_b, f_a, f_E, r1, r2, r3, mass, radius = theta_in
    #    X, Z, Q_b, s_t, f_a, f_E, r1, r2, r3 = theta

    # by default we assume the model is valid, i.e. has sufficient bursts
    # to match the observations, AND doesn't violate the GTI conditions 
    # (if we are checking those)

    valid = True

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
            # this approach is a bit newer and hopefully more robust, also
            # reflects the simulation approach in generate_burst_train

            imatch = [np.argmin(np.abs(tpred - bstart[ref_ind]))]
            for i in range(1,numburstssim):
                # now looking for a match for bursts ref_ind-i, ref_ind+i:
                # (but also want to exclude any that have already been 
                # matched!)
                if ref_ind-i >= 0:
                    # imatch.insert(0,np.argmin(np.abs(tpred-bstart[ref_ind-i])))
                    imatch.insert(0, next_nearest(bstart[ref_ind-i], tpred, imatch))
                if ref_ind+i < numburstsobs:
                    # imatch.append(np.argmin(np.abs(tpred-bstart[ref_ind+i])))
                    imatch.append(next_nearest(bstart[ref_ind+i], tpred, imatch))

            if np.any(np.array(imatch) == None):
                # Not enough bursts to match, so count model as invalid
                return None, False, result

            assert len(imatch) == numburstsobs
            # print (imatch)
            assert len(set(imatch)) == numburstsobs # don't double up on bursts

            # We only compare the times of the bursts for the events excluding
            # the reference burst, from which the train is calculated
            itime = imatch.copy()
            itime.pop(ref_ind)

            # We compare the fluences for all the bursts
            # We subtract 1 here, and also to the expression for ialpha, because the indexing is different for
            # the e_b and alpha arrays, compared to the times, coming out of generate_burst_train
            ie_b = [np.max([x - 1, 0]) for x in imatch]

            # We compare the alpha values only for observed events #1, 2 & 3 as we don't have the recurrence
            # time for event #0
            ialpha = ie_b.copy()
            ialpha.pop(0)

            # Now assemble the whole array for comparison with the measurements
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

            # I don't think it's ever correct to swap model for result
            # anymore, particularly since we're also returning the result;
            # just set the model as invalid and bail out

            # model = result
            return None, False, result

    else:
        # If we're not generating a burst train, just run the ensemble

        result = burstensemble(
            Q_b, X, Z, r1,r2,r3,mass,radius,bstart,pflux,numburstsobs)

        model = np.concatenate((result['time'], result['e_b'], result['alpha']))

    # Check here if the model instance is valid, i.e. the bursts that are NOT
    # matched with the observed ones must fall in gaps
    # Originally we used the (global) arrays st, et defined by 1808-match, to
    # avoid copying them over from IDL each time; but now these are passed
    # as parameters gti_start, gti_end

    if gti_checking:
        # if "st" not in globals():
        if (gti_start is None) or (gti_end is None):
            print ('** WARNING ** can''t access GTI information')
            return model, valid, result
        else:
            st, et = gti_start, gti_end

        for index, rt in enumerate(tpred):
            if index not in imatch:
                # ok, not one of the known bursts. Is it an excluded time?#
                for i in range(len(st)):

                    if rt >= st[i] and rt <= et[i] - 10.0 / 86400.0:

                        return model, False, result

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
