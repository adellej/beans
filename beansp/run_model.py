# -------------------------------------------------------------------------#
import numpy as np

# load local module
from .burstrain import *

def burst_time_match(iref1, time1, iref2, time2):
    """
    Function to generate an array of indices of elements in time2 that are
    closest (overall) to the elements of time1. One of the elements of
    time1 (index iref1) is chosen as a "reference" that should be
    replicated exactly in time2 (index iref2)

    :param iref1: index of "reference" time in array time1
    :param time1: times to match. In run_model, this is the observed burst times
    :param iref2: index of "reference" time in array time2
    :param time2: times to be matched; i.e. the predicted burst times

    :return: list of indices having same length as time1
    """

    def match_left(ix, iref, time1, time2, first=None):
        """
        Routine to determine nearest values in array time2 for the 0 to ix[0]'th
        elements of array time1, with the first guess optionally given
        """

        for _i in np.arange(iref,0,-1)-1:
            # print (_i, iref, time1[_i], time2[:ix[0]], np.argmin(np.abs(time1[_i]-time2[:ix[0]])))
            if (_i == iref-1) & (first is not None): # & (first < ix[0]):
                ix.insert(0, first)
            elif ix[0] <= 0:
                break
            else:
                ix.insert(0, np.argmin(np.abs(time1[_i]-time2[:ix[0]])))   
    
        return ix

    def match_right(ix, iref, time1, time2, first=None):
        """
	Routine to determine nearest values in array time2 for the
	ix[-1]'th to last elements of array time1, with the first guess
        optionally given
        """
    
        for _i in np.arange(iref,len(time1)-1)+1:
            # print (_i, time1[_i], time2[ix[-1]+1:], ix[-1]+1, len(time2))#np.argmin(np.abs(time1[_i]-time2[ix[-1]+1:])))
            if (_i == iref+1) & (first is not None):
                ix.append(first)
            elif ix[-1]+1 == len(time2):
                break
            else:
                ix.append( np.argmin(np.abs(time1[_i]-time2[ix[-1]+1:]))+ix[-1]+1 )

        return ix

    assert np.isclose(time1[iref1],time2[iref2],rtol=1e-4)

    # special here for IGR J17511-3057, forcing the match solution
    # this for continuing the base20 run
    # ioff = [-17, -16, -13, -11,  -6,  -5,  -3,  -2,  -1,   0,   3,   4,   5,
    #      7,   9,  11,  12,  14,  15,  19]
    # & for IGR J17498-2921, for the forcepm9 run
    # ioff = [-9, -5, -2, -1,  0,  1,  8,  9]
    # if (iref2+ioff[0] >= 0) & (len(time2) > iref2+ioff[-1]):
    #     return [iref2+x for x in ioff]
    # else:
    #     return None

    # make sure we have enough bursts to match, on either side of the
    # reference
    if ((iref2 < iref1) | (len(time2)-iref2-1 < len(time1)-iref1-1)):
        return None

    ix = [iref2]
    ix = match_left(ix, iref1, time1, time2)

    ix = match_right(ix, iref1, time1, time2)

    if len(ix) < len(time1):
        return None

    # now do the refinement stage. When you have skipped bursts (i.e.
    # predicted bursts that are not observed) it's possible that the
    # nearest burst for one in the sequence away from the ends, will force
    # a sub-optimal choice further out (from the ref) in the train. So
    # these sections are designed to test for that

    # print (ix, np.sum(np.abs(time1-time2[ix])))

    _i = 0
    while _i < iref1-1:
        # print (_i, ix[_i+1]-ix[_i], ix[_i+2]-ix[_i+1])
        if (ix[_i+1]-ix[_i] == 1) & (ix[_i+2]-ix[_i+1] > 1):
            ix_try = match_left(ix[_i+2:], _i+2, time1, time2, first=ix[_i+1]+1)
            # print (ix_try, np.sum(np.abs(time1-time2[ix_try])))
            if np.sum(np.abs(time1-time2[ix_try])) < np.sum(np.abs(time1-time2[ix])):
                ix = ix_try
                _i = -1 # so as to start from the start again, after _i is incremented below
        _i += 1

    # print ('after leftward refinement: ',ix)

    _i = len(ix)-1
    while _i > iref1+1:
        # print (_i, ix[_i]-ix[_i-1], ix[_i-1]-ix[_i-2], ix[:_i-1])
        if (ix[_i]-ix[_i-1] == 1) & (ix[_i-1]-ix[_i-2] > 1):
            ix_try = match_right(ix[:_i-1], _i-2, time1, time2, first=ix[_i-1]-1)
            # print (ix_try, np.sum(np.abs(time1-time2[ix_try])))
            if np.sum(np.abs(time1-time2[ix_try])) < np.sum(np.abs(time1-time2[ix])):
                ix = ix_try
                _i = len(ix) # so as to start from the end again, after _i is decremented below
        _i -= 1

    # print ('after rightward refinement: ',ix)

    return ix


def runmodel(theta_in, bean, match=True, debug=False):
    """
    This routine calls one of two functions that generate the burst model
    predictions, either generate_burst_train or burstensemble, depending on
    which mode the analysis is in
    It then assembles the predictions from the model into a form that can
    be compared to the observations (with appropriate scaling)

    :param theta_in: parameter vector, with *X*, *Z*, *Q_b*, *d*, *xi_b*,
      *xi_p*, and (optionally) *mass*, *radius*, *f_E* & *f_a*
    :param bean: Beans object, from which the required parameters are drawn:
      y, tref, bstart, pflux, pfluxe, tobs, numburstssim, numburstsobs,
      ref_ind, gti_checking,train, gti_start=None, gti_end=None,
    :param match: set to False to skip the burst matching stage, which is useful for exploratory/diagnostic purposes
    :param debug: set to True to display more diagnostic information, passed also to generate_burst_train

    :return: model array with predicted burst parameters, Boolean giving
      the validity of the solution (i.e. consistent with the GTI information),
      and the full result dict from generate_burst_train/burstensemble
    """

    if debug:
        print('Calling runmodel')

    X, Z, Q_b, dist, xi_b, xi_p, *extra = theta_in
    mass, radius, f_E, f_a = extra+[bean.M_NS, bean.R_NS, 1.0, 1.0][len(extra):]

    # by default we assume the model is valid, i.e. has sufficient bursts
    # to match the observations, AND doesn't violate the GTI conditions
    # (if we are checking those)

    valid = True

    if bean.train:

        # Now call the function. From the code:
        # This routine generates a simulated burst train. The output is a
        # dict with keys corresponding to the model parameters and
        # predicted values.
        # We now pass on the debug flag to help with exploratory work

        result = generate_burst_train( bean,  Q_b, X, Z, dist, xi_p, mass, radius, debug=debug)

        # we need to reject unsuitable models here. The simplest (although
        # insufficiently strict) criterion is to at least simulate as many
        # bursts as are observed.

        tpred = result["time"]
        npred = len(tpred)
        if npred < bean.numburstsobs:
            if debug:
                print ('runmodel: insufficient number of bursts simulated ({} of {}'.format(npred, bean.numburstsobs))
            return None, False, result

        # bug testing sample model values
        # model = np.array([-0.04106, 2.77858, 3.99188, 3.73303, 3.68542, 4.16907, 4.71480, 115.000, 126.903, 138.070]) #to test the code we define the model as an array, where the numbers in the array are values for the parameters of y, in the same order as y. Values have been taken from Duncan's example code output

    else:
        # If we're not generating a burst train, just run the ensemble

        result = burstensemble( bean, Q_b, X, Z, dist, xi_p, mass, radius)

    # Here we convert the model-predicted values to observational
    # quantities, for comparison with the observations. Previously
    # this was achieved via the "ratios" of the mdot/luminosity,
    # fluence and alphas
    # The times are already relative to the reference bursts, so
    # nothing needs to be done to those

    result['fluen'] = result['e_b'] * (bean.fluen_fac/xi_b/dist**2).value

    result['alpha_obs'] = result['alpha'] * xi_b/xi_p

    # And finally form the model array for comparison with the data

    model = None
    if not bean.train:
        # For the burstensemble version, this is simple
        model = np.concatenate((result['time'], result['fluen'],
            result['alpha_obs']))

    elif match & (bean.y is not None) & (npred > 0):

        # Assemble the array for comparison with the data
        # First need to match the burst times; this is surprisingly
        # difficult to do robustly.

        imatch = burst_time_match(bean.ref_ind, bean.bstart, result['iref'], tpred)

        if imatch is None:
            return None, False, result

        result['imatch'] = imatch

        # print (imatch)
        if len(imatch) != bean.numburstsobs:
            print ('** ERROR ** likely insufficient number of bursts simulated ({} < {})'.format(
                len(imatch), bean.numburstsobs))

            return model, False, result

        else:

            assert len(set(imatch)) == bean.numburstsobs # don't double up on bursts

            # We only compare the times of the bursts for the events excluding
            # the reference burst, from which the train is calculated
            itime = imatch.copy()
            itime.pop(bean.ref_ind)

            # We compare the fluences for all the bursts
            # We subtract 1 here, and also to the expression for ialpha, because the indexing is different for
            # the e_b and alpha arrays, compared to the times, coming out of generate_burst_train
            ie_b = [np.max([x - 1, 0]) for x in imatch]

            # We compare the alpha values only for observed events #1, 2 & 3 as we don't have the recurrence
            # time for event #0
            ialpha = ie_b.copy()
            ialpha.pop(0)

            # Now assemble the whole array for comparison with the measurements
            # Includes the selection of non-zero fluences (set in get_obs)

            model =  np.concatenate((result['time'][itime],
                result['fluen'][ie_b][bean.ifluen],
                result['alpha_obs'][ialpha][bean.ifluen[1:]-1]))

    # Check here if the model instance is valid, i.e. the bursts that are NOT
    # matched with the observed ones must fall in gaps
    # Originally we used the (global) arrays st, et defined by 1808-match, to
    # avoid copying them over from IDL each time; but now these are stored
    # with the Beans object

    if match & (model is not None) & bean.gti_checking:
        # if "st" not in globals():
        if (bean.st is None) or (bean.et is None):
            print ('** WARNING ** can''t access GTI information')
            return model, valid, result

        # this is *very* inefficient, but might be the simplest option,
        # particularly if we have overlapping GTIs
        for index, rt in enumerate(tpred):
            if index not in imatch:
                # ok, not one of the known bursts. Is it an excluded time?#
                for i in range(len(bean.st)):

                    if (rt >= bean.st[i]) and (rt <= bean.et[i] - 10.0 / 86400.0):

                        # print (rt, i, bean.st[i], bean.et[i])
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
