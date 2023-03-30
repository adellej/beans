""" Optional routine to calculate the alpha values if they have not been precalculated """

 # # Calculate the alphas:

    # alpha = np.zeros(len(bstart)-1)
    # alpha_err = np.zeros(len(bstart)-1)
    # tdel = np.zeros(len(bstart)-1)
    # for i in range (1,len(bstart)):

    #     #alpha = Fp cbol delta t /fluence
    #     if i == 1:
    #         alpha[i-1]=alpha[i-1]/3.  # special for first interval
    #         tdel[i-1] = bstart[i]-bstart[i-1]
    #     else:
    #         alpha[i-1] = (mflux[i-1]*(bstart[i]-bstart[i-1])*86400.00)/(fluen[i-1]*1000.0)
    #         tdel[i-1] = bstart[i]-bstart[i-1]

    # # Correct for missed bursts:

    # for (j, i) in zip(N_missed, burst_index):
    #     alpha[i] = alpha[i]/j

    # for i in range (1,len(bstart)):
    #     alpha_err[i-1]=alpha[i-1]*math.sqrt(np.mean(bce/bc)**2 + np.mean(pfluxe/pflux)**2 + (fluene[i-1]/fluen[i-1])**2)


    # # Not sure what's happening with alpha errors, for now just set alpha_err = 5

    # alpha_err = np.empty(len(alpha))
    # alpha_err.fill(5.0)
