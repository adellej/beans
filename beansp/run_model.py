# -------------------------------------------------------------------------#
import numpy as np

# load local module
from .burstrain import *

def runmodel(theta_in, bean, debug=False):
    """
    This routine calls one of two functions that generate the burst model
    predictions, either generate_burst_train or burstensemble, depending on
    which mode the analysis is in
    It then assembles the predictions from the model into a form that can
    be compared to the observations (with appropriate scaling)

    :param theta_in: parameter tuple: X, Z, Q_b, dist, xi_b, xi_p, mass,
      radius, and (optionally) f_a & f_E
    :param bean: Beans object, from which the required parameters are drawn:
      y, tref, bstart, pflux, pfluxe, tobs, numburstssim, numburstsobs,
      ref_ind, gti_checking,train, gti_start=None, gti_end=None,
    :param debug: set to True to display more diagnostic information
    :return: model array with predicted burst parameters, Boolean giving
      the validity of the solution (i.e. consistent with the GTI information),
      and the full result dict from generate_burst_train/burstensemble
    """

    if debug:
        print('Calling runmodel')

    X, Z, Q_b, dist, xi_b, xi_p, mass, radius = theta_in[:8]

    # by default we assume the model is valid, i.e. has sufficient bursts
    # to match the observations, AND doesn't violate the GTI conditions
    # (if we are checking those)

    valid = True

    if bean.train:

        # Now call the function. From the code:
        # This routine generates a simulated burst train. The output is a
        # dict with keys corresponding to the model parameters and
        # predicted values.

        result = generate_burst_train( bean,  Q_b, X, Z, dist, xi_p, mass, radius)

        # we need to reject unsuitable models here. The simplest (although
        # insufficiently strict) criterion is to at least simulate as many
        # bursts as are observed.

        tpred = np.array(result["time"])
        npred = len(tpred)
        if npred < bean.numburstsobs:
            return None, False, result

        # bug testing sample model values
        # model = np.array([-0.04106, 2.77858, 3.99188, 3.73303, 3.68542, 4.16907, 4.71480, 115.000, 126.903, 138.070]) #to test the code we define the model as an array, where the numbers in the array are values for the parameters of y, in the same order as y. Values have been taken from Duncan's example code output

        if (bean.y is not None) & (npred > 0):

            # Assemble the array for comparison with the data
            # First the burst times; we dynamically determine which
            # match offers the smallest error in absolute terms, from the 
            # available possibilities

            ref_tpred = np.argmin(np.abs(tpred - bean.bstart[bean.ref_ind]))
            if ((ref_tpred < bean.ref_ind) |
                (npred-ref_tpred-1 < bean.numburstsobs-bean.ref_ind-1)):
                return None, False, result

            # i, j, ref_ind, ref_tpred = 5, 10, 2, 6

            # it likely is not necessary to generate the match array every
            # time. You could store the matches array and check if it needs
            # updating, e.g. with the following:

            calc_matches = True
            if bean.matches is not None:
                if (set(bean.matches[:,bean.ref_ind]) == set([ref_tpred]) and 
                    (np.max(bean.matches) == npred-1)):
                    calc_matches = False

            if calc_matches:
                # initial match array, with the earliest bursts before and
                # after the reference
                imatch = np.concatenate( (np.arange(bean.ref_ind),
                    np.arange(bean.numburstsobs-bean.ref_ind)+ref_tpred ) )
                # final match array, with the *latest* possible bursts
                jmatch = np.concatenate( (np.arange(ref_tpred-bean.ref_ind, ref_tpred+1),
                    np.arange(npred-(bean.numburstsobs-bean.ref_ind-1), npred) ) )

                # accumulate the arrays
                matches = np.expand_dims(imatch.copy(), axis=0)

                # now loop and generate the other possible match combinations
                while ((imatch[0] < ref_tpred-bean.ref_ind) |
                    (imatch[bean.ref_ind+1] < npred-(bean.numburstsobs-bean.ref_ind-1))):
                    imatch[-1] += 1
                    while np.any(imatch > jmatch):
                        ex_ind = min(np.where(imatch > jmatch)[0])
                        # print (ex_ind, i-ex_ind, ix[ex_ind-1])
                        # print ('in loop:',imatch[ex_ind-1:], np.arange(bean.numburstsobs-ex_ind+1), imatch[ex_ind-1]+1)
                        imatch[ex_ind-1:] = np.arange(bean.numburstsobs-ex_ind+1) \
                            + imatch[ex_ind-1]+1
                        if ex_ind < bean.ref_ind:
                            imatch[bean.ref_ind:] = np.arange(bean.numburstsobs-bean.ref_ind) \
                                + ref_tpred
                        # special to keep the index preserved here
                        if imatch[bean.ref_ind] > ref_tpred:
                            imatch[bean.ref_ind-1] += 1
                            imatch[bean.ref_ind:] = np.arange(bean.numburstsobs-bean.ref_ind) \
                                + ref_tpred
                    # print (ix, sum(abs(bstart-np.array(tpred)[ix])))
                    matches = np.append(matches, [imatch], axis=0)
         
                bean.matches = matches

            # print ("tpred etc:",tpred, matches, bean.bstart)
            _err = np.sum((tpred[bean.matches]-bean.bstart)**2, axis=1)

            imatch = bean.matches[np.argmin(_err),:]
            result['imatch'] = imatch
 
            assert len(imatch) == bean.numburstsobs
            # print (imatch)
            assert len(set(imatch)) == bean.numburstsobs # don't double up on bursts

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

        result = burstensemble( bean, Q_b, X, Z, dist, xi_p, mass, radius)

    # Here we convert the model-predicted values to observational
    # quantities, for comparison with the observations. Previously
    # this was achieved via the "ratios" of the mdot/luminosity,
    # fluence and alphas
    # The times are already relative to the reference bursts, so
    # nothing needs to be done to those

    result['fluen'] = list( np.array(result['e_b']) 
        * (bean.fluen_fac/xi_b/dist**2).value )

    result['alpha_obs'] = list(np.array(result['alpha']) * xi_b/xi_p )

    # TODO here also you could modify the predicted values, e.g. to make them
    # more consistent with Kepler (for example). This could be tricky for
    # the times, which are already relative to the start

    pass

    # And finally form the model array for comparison with the data

    if bean.train:

        # We only compare the times of the bursts for the events excluding
        # the reference burst, from which the train is calculated
        itime = list(imatch.copy())
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
        # this would be simpler if the result elements were numpy arrays
        # and not lists!

        model =  np.concatenate((np.array(result['time'])[itime],
            np.array(result['fluen'])[ie_b],
            np.array(result['alpha_obs'])[ialpha]))

    else:
        # model = np.concatenate((result['time'], result['e_b'], result['alpha']))
        model = np.concatenate((result['time'], result['fluen'],
            result['alpha_obs']))

    # Check here if the model instance is valid, i.e. the bursts that are NOT
    # matched with the observed ones must fall in gaps
    # Originally we used the (global) arrays st, et defined by 1808-match, to
    # avoid copying them over from IDL each time; but now these are stored
    # with the Beans object

    if bean.gti_checking:
        # if "st" not in globals():
        if (bean.st is None) or (bean.et is None):
            print ('** WARNING ** can''t access GTI information')
            return model, valid, result

        for index, rt in enumerate(tpred):
            if index not in imatch:
                # ok, not one of the known bursts. Is it an excluded time?#
                for i in range(len(bean.st)):

                    if rt >= bean.st[i] and rt <= bean.et[i] - 10.0 / 86400.0:

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
