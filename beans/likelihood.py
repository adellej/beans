""" module that defines the likelihood and priors required by MCMC """

## Now we define the functions that emcee requires
# define likelihood as a function of theta, x, y and yerr as this is what emcee expects as the inputs

def lnlike(theta, x, y ,yerr):

# define y = "data" parameters

    for x,i in zip([ x for x in range(0, len(bstart)-1) if x != ref_ind], [i for i in range(0, len(bstart)-1) if i != ref_ind]):
        globals()['t%s' % i] = y[x]
    for x,i in zip(range(len(bstart)-1, len(fluen)+len(bstart)-1),range(0,len(bstart))):
        globals()['Eb%s' % i] = y[x]
    for x,i in zip(range(len(fluen)+len(bstart)-1, len(y)),range(0, len(bstart-1))):
        globals()['a%s' % i] = y[x]

# define yerr as variance terms (errors) for our data parameters (listed in same order as for y)
# *note that we have to enter three time errors for the code to work however in reality the error should be the same for all of them (st0, st2 and st3 are really dummy parameters)

    for x,i in zip([ x for x in range(0, len(bstart)-1) if x != ref_ind], [i for i in range(0, len(bstart)-1) if i != ref_ind]):
        globals()['st%s' % i] = yerr[x]
    for x,i in zip(range(len(bstart)-1, len(fluen)+len(bstart)-1),range(0,len(bstart))):
        globals()['sEb%s' % i] = yerr[x]
    for x,i in zip(range(len(fluen)+len(bstart)-1, len(y)),range(0, len(bstart-1))):
        globals()['sa%s' % i] = yerr[x]

# define theta = model parameters, which we define priors for

    X, Z, Q_b, f_a, f_E, r1, r2, r3, mass, radius = theta

# Instead of treating s_t as a parameter, we just hardwire it here

    s_t = 10./1440.

# call model from IDL code defined as modeldata(base, z, x, r1, r2 ,r3)
    model, valid = runmodel(theta,y,tref,bstart,pflux, pfluxe, tobs,numburstssim,ref_ind,gti_checking)

    if not valid:
        return -np.inf, model

# multiplying by scaling factors to match with the data
    model[len(bstart)-1:len(fluen)+len(bstart)-1] *= r3
    model[len(fluen)+len(bstart)-1:len(y)] *= r2

# To simplify final likelihood expression we define inv_sigma2 for each data parameter that describe the error.
# The variance (eg sEb0) is underestimated by some fractional amount, f, for each set of parameters.

    sEb = yerr[len(bstart)-1:len(fluen)+len(bstart)-1]
    sa = yerr[len(fluen)+len(bstart)-1:len(yerr)]

    inv_sigma2 = []
    for i in range (0,len(bstart)-1):
        inv_sigma2.append(1.0/(s_t**2))
    for i in range(0,len(bstart)):
        inv_sigma2.append(1.0/((sEb[i]*f_E)**2))
    for i in range(0,len(bstart)-1):
        inv_sigma2.append(1.0/((sa[i]*f_a)**2))

# Final likelihood expression
    cpts = (y - (model))**2*inv_sigma2 - (np.log(inv_sigma2))

# Test if the result string is defined here. It is, so we return the selected elements of result
# instead of the downselection in model


    base=Q_b
    z=Z
    x=X
    r1=r1
    r2=r2
    r3=r3
    mass=mass
    radius=radius

    model2 = generate_burst_train(base,z,x,r1,r2,r3,mass, radius,bstart,pflux, pfluxe, tobs,numburstssim,run=run,double=double,debug=debug)

# Now also return the model
    return -0.5*np.sum( cpts ), model2

# -------------------------------------------------------------------------#
# This is the expression for the prior on the model parameters. These are pretty uninformative

# Define priors for theta. mr prior function is located in mrprior.py

def lnZprior(z):
    # This beta function for the metallicity prior is from Andy Casey and is an approximation of the metallicity of a mock galaxy
    # at 2.5-4.5 kpc for the location of 1808. Assuming ZCNO = 0.01 is average value.
    from scipy import stats
    import numpy as np
    beta = stats.beta
    ZCNO = 0.01
    #if z == 0.:
    #    return -np.inf
    #else:
    return(np.log(beta(10,3).pdf((np.log10(z/ZCNO)+3)/3.75)/(3.75*np.log(10)*z)))

def lnprior(theta):
    import numpy as np
    X, Z, Q_b, f_a, f_E, r1, r2, r3, mass, radius = theta


    if 0.00001 < X < 0.76 and 0.00001 < Z < 0.056 and 0.000001 <= Q_b < 5.0 and 0 < f_a < 100 and 0 < f_E < 100 and   \
      0.005 < r1 < 1.0 and 0.005 < r2 < 3.0 and 0 < r3*1e3 < 1000 and 1.15 < mass < 2.5 \
      and 9 < radius < 17: #upper bound and lower bounds of each parameter defined here. Bounds were found by considering an estimated value for each parameter then giving reasonable limits.
        return 0.0 + lnZprior(Z) + mr_prior(mass, radius)
    return -np.inf

# -------------------------------------------------------------------------#
# Finally we combine the likelihood and prior into the overall lnprob function, called by emcee

# define lnprob, which is the full log-probability function
def lnprob(theta, x, y, yerr):
    import numpy as np
    lp = lnprior(theta)

# Now also returns the model, to accumulate along with the likelihoods

    like, model = lnlike(theta, x, y ,yerr)

    if (not np.isfinite(lp)) or (not np.isfinite(like)):
        return -np.inf, model

    return lp + like, model
# -------------------------------------------------------------- #