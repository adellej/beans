# -------------------------------------------------------------------------#
import numpy as np
from burstrain import *

def runmodel(theta_in, y, tref, bstart, pflux, pfluxe, tobs, numburstssim, ref_ind, gti_checking):

    X, Z, Q_b, f_a, f_E, r1, r2, r3, mass, radius = theta_in
    #    X, Z, Q_b, s_t, f_a, f_E, r1, r2, r3 = theta

    # Set the imput parameters to generate_burst_train
    # These variables need to be passed to idl, hence the "idl." in the definition
    # statement

    base = Q_b
    z = Z
    x = X
    r1 = r1
    r2 = r2
    r3 = r3
    mass = mass
    radius = radius

    # Now call the function. From the code:
    # This routine generates a simulated burst train. The output is a
    # structure with the following elements:
    #   BASE            FLOAT          0.175000
    #   Z               FLOAT         0.0100000
    #   X_0             FLOAT          0.440000
    #   R1              FLOAT          0.108533
    #   R2              FLOAT           1.00000
    #   R3              FLOAT           1.00000
    #   MDOT            DOUBLE    Array[7]
    #   MDOT_MAX        DOUBLE         0.043402292
    #   TIME            DOUBLE    Array[8]
    #   ALPHA           FLOAT     Array[7]
    #   E_B             FLOAT     Array[7]
    #   MASS            FLOAT     1.4
    #   RADIUS          FLOAT     11.2

    result = generate_burst_train(
        base, z, x, r1, r2, r3, mass, radius, bstart, pflux, pfluxe, tobs, numburstssim
    )

    tpred = result["time"]

    # bug testing sample model values
    # model = np.array([-0.04106, 2.77858, 3.99188, 3.73303, 3.68542, 4.16907, 4.71480, 115.000, 126.903, 138.070]) #to test the code we define the model as an array, where the numbers in the array are values for the parameters of y, in the same order as y. Values have been taken from Duncan's example code output

    # Assemble the array for comparison with the data

    i1 = []
    for i in range(0, ref_ind):
        i1.append(np.argmin(np.abs(tpred - y[i])))

    i1.append(np.argmin(np.abs(tpred - tref)))

    for i in range(ref_ind, len(bstart) - 1):
        i1.append(np.argmin(np.abs(tpred - y[i])))

    li1 = list(i1)
    li1m1 = [np.max([x - 1, 0]) for x in li1]

    # We compare the fluences for all the bursts
    # We add 1 here, and also to the expression for ialpha, because the indexing is different for
    # the e_b and alpha arrays, compared to the times, coming out of generate_burst_train

    ie_b = li1m1

    # We only compare the times of the bursts for observed events #0, 2 & 3; #10 is a "reference"
    # from which the train is calculated

    i2 = li1
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

    # Check here if the model instance is valid, i.e. the bursts that are NOT matched with the
    # observed ones must fall in gaps
    # We use the (global) arrays st, et defined by 1808-match, to avoid copying them over from IDL
    # each time
    valid = True
    if gti_checking is 1:
        if "st" not in globals():
            return model, valid

        for index, rt in enumerate(tpred):
            if index not in i1:
                # ok, not one of the known bursts. Is it an excluded time?#
                for i in range(len(st)):

                    if rt >= st[i] and rt <= et[i] - 10.0 / 86400.0:

                        valid = False
                        return model, valid

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

    return model, valid


# -------------------------------------------------------------------------#
