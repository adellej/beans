# -------------------------------------------------------------------------#

def settle(base, z, x_0, mdot, cfac, mass, radius):
    import settler as se
    import numpy as np

    # initialize settle interface
    settl = se.settle()

    # run settle:
    res = settl.full(Z=z, X=x_0, M=mdot, F=base, C=0, R=radius, Ma=mass)

    # extract results for comparison with obs
    result = np.recarray((1,),dtype=[('tdel', np.float64), ('E_b', np.float64),('alpha', np.float64)])
    # assign elements
    result.tdel = res[0]*cfac
    result.E_b = res[2]*cfac/1e39
    result.alpha = res[1]
    result.mdot = mdot

    print(result.tdel, result.E_b, result.alpha, result.mdot)
    return result




# def settle(base, z, x_0, mdot, cfac, mass, radius, binpath):
#
#     import numpy as np
#     import subprocess
#
#     compress = 0 # Compressional heating option (0 is off 1 is on)
#     #dist=dist #(kpc) if this is set, convert fluence to 1e-9 erg/cm^2
#
#     # I think these are unnecessary
#
#     #if n_elements(compress) is 0:
#     #    compress=0
#     #if n_elements(_base) is 0:
#     #    base = 0.1
#     #if n_elements(_z) is 0:
#     #    z = 0.016
#     #if n_elements(_x_0) is 0:
#     #    _x_0 = 0.7
#
#     #if n_elements(dist) is 0:
#     dfac=1.0
#     #else:
#      #   dfac==1e9*8.356e-6/dist^2
#
#     args = ("%f %f %f %f %f %f %f" % (base, z, x_0, mdot, radius, mass, compress))
#
#     # Run settle using subprocess, settle takes arguments as follows:
#     # settle <Q_base> <z> <x> <mdot> <compress>
#     # e.g. settle 0.1 0.012 0.6 0.11 0
#
#     # Call settle at location binpath with arguments args and get output:
#
#     output = subprocess.check_output(["%s %s"%(binpath, args)],shell=True)
#     output = str(output)
#     # Returns a string array like this:
#     #  1: Gravity g=2.45e+14
#     #  2: Radius = 10 km
#     #  3: Q=900 in the crust
#     #  4: I won't integrate the magnetic profile
#     #  5: Setting base flux = 0.1
#     #  6: Setting metallicity = 0.016
#     #  7: Setting accreted H fraction X=0.7
#     #  8: I get mdot/Edd=0.1
#     #  9: Depletion column=3.98603e+08
#     # 10: No compressional heating!
#     # 11: in cgs, Fb=8.4796e+20
#     # 12:
#     # 13: Searching for ignition depth...
#     # 14:
#     # 15: ------------- Ignition conditions ----------------------------
#     # 16:          Z       mdot          T          y          P          Y         X        rho       flux
#     # 17:      0.016       8800  2.099e+08  1.637e+08  4.009e+22     0.5717 0.4123   7.87e+05  1.605e+22
#     # 18: ---------------------------------------------------------------
#     # 19: z=-802.029 I=-1.37085e+33 I_0=1.37099e+33  delta I/I_0=-0.000101115
#     # 20: Recurrence time =6.76714 hours
#     # 21: Xbar=0.556138, Q=3.82455, Energy=5.79449e+39
#     # 22: eps (14C+alpha) is 5.9416e+10
#     # 23: eps (3a) is 1.14051e+13
#     # 24: kappa=0.0514806
#     # 25: t_alpha=8.13199 hours
#     #
#     # So then just have to get the right numbers out:
#
#
#     # Separate settle output into lines to extract values, start line counting from "ignition conditions"
#     output=output.split('Ignition conditions')
#     output = output[1]
#     lines = str(output).split('\\n')
#
#
#     # Extract xbar, qnuc and energy from settle output:
#
#     # These are all found in line 6 of the output
#     l = lines[6].split('=')
#
#     l1 = str(l[1]).split(',')
#     l2 = str(l[2]).split(',')
#     l3 = str(l[3]).split(',')
#
#     xbar = l1[0]
#     Q_nuc = l2[0]
#     Energy = l3[0]
#
#     # need to make energy and xbar floats as they are currently strings and
#     # can't be multiplied with other numbers to make alpha and E_b
#     xbar = np.float64(xbar)
#     Energy = np.float64(Energy)
#
#     alpha = 290/(1.6 + 4.0*xbar)
#     E_b = (dfac*Energy)*cfac/1e39
#
#     # Now get recurrence time, found in line 5
#     t = lines[5].split('=')
#     t1 = str(t[1]).split(' ')
#
#     tdel = t1[0]
#     tdel = np.float64(tdel)
#     tdel = tdel*cfac
#
#     # Now get y8, t8, x_b and y_b, found in line 2
#
#     m = lines[2]
#     m = ' '.join(m.split())
#     m = m.split(' ')
#
#     y8 = m[3]
#     y8 = np.float64(y8)
#     y8 = y8/1e8
#     t8 = m[2]
#     t8 = np.float64(t8)
#     t8 = t8/1e8
#     y_b = m[5]
#     x_b = m[6]
#
#     #print y8, t8, y_b, x_b
#
#     # Now assemble result array so elements can be accessed using e.g. result.tdel (use numpy record array)
#
#     # create array
#     result = np.recarray((1,),dtype=[('tdel', np.float64), ('y8', np.float64), ('t8', np.float64), ('x_b', np.float64),                                     ('y_b', np.float64),('Q_nuc', np.float64), ('E_b', np.float64),('alpha', np.float64),                                     ('xbar', np.float64)])
#     # assign elements
#     result.tdel = tdel
#     result.y8 = y8
#     result.t8 = t8
#     result.x_b = x_b
#     result.y_b = y_b
#     result.Q_nuc = Q_nuc
#     result.E_b = E_b
#     result.alpha = alpha
#     result.xbar = xbar
#
#     # done :D
#     return result #result array is in the order: tdel, y8, t8, x_b, y_b, Q_nuc, E_b, alpha, xbar
