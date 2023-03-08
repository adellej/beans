from __future__ import division
from __future__ import print_function

import ctypes as ct
import numpy
import pathlib


class Settle(object):
    """
    Super basic interface to settle code
    """

    def __init__(self, F=0.1, C=0):
        """
        Will just return a convenient object to call settle

        Has methods:

        full: calls settle with all arguments, in case you want to change paramenters,
              it overrides the default F and C, but does not overwrite them
              Takes f, m, x, z, c
              Returns alpha, trec, fluence
        run : more compact call, which assumes f and c from initialization.
              Takes m, x, and z
              Returns alpha, trec, fluence
        """
        path_to_data_file = (
            pathlib.Path(__file__).resolve().parent.parent / "settle" / "libsettle.so"
        )

        self.libsettle = ct.cdll.LoadLibrary("libsettle.so")

        print("o", end="")

        self.mainer = self.libsettle.mainer

        self.mainer.argtypes = [
            ct.POINTER(ct.c_double),
            ct.POINTER(ct.c_double),
            ct.POINTER(ct.c_double),
            ct.POINTER(ct.c_double),
            ct.POINTER(ct.c_int),
            ct.POINTER(ct.c_double),
        ]
        self.mainer.returntype = ct.c_int

        self.init_vars(F, C)

    def init_vars(self, F=0.1, C=0):
        """
        Takes:

        F: [default: 0.1] the flux from the bottom
        C: [default: 0  ] 0 or 1: to do compressional heating or not

        Has attributes:

        F, C     : the default parameters
                   (so that you do not specify these all the times, just for lazyness)
        libsettle: link to the settle code library
        mainer   : the real callable of the library (should not be used directly by user)
        """
        self.F = ct.c_double(F)
        self.C = ct.c_int(C)

    def run(self, M, X, Z, R, Ma):
        """
        Runs settle, you CANNOT specify F and C,
        because it uses the defaults you specified at creation.

        Can pass either scalars, or equally long arrays.
        """

        T = ct.c_double()
        A = ct.c_double()
        E = ct.c_double()

        if hasattr(M, "__iter__"):

            resA = []
            resR = []
            resE = []

            for i in range(len(M)):

                ret = self.mainer(
                    ct.byref(self.F),
                    ct.byref(ct.c_double(Z[i])),
                    ct.byref(ct.c_double(X[i])),
                    ct.byref(ct.c_double(M[i])),
                    ct.byref(self.C),
                    ct.byref(T),
                    ct.byref(A),
                    ct.byref(E),
                    ct.byref(ct.c_double(R[i])),
                    ct.byref(ct.c_double(Ma[i])),
                )

                resA.append(A.value)
                resR.append(T.value)
                resE.append(E.value)

            return numpy.array(resA), numpy.array(resR), numpy.array(resE)

        else:
            ret = self.mainer(
                ct.byref(self.F),
                ct.byref(ct.c_double(Z)),
                ct.byref(ct.c_double(X)),
                ct.byref(ct.c_double(M)),
                ct.byref(self.C),
                ct.byref(T),
                ct.byref(A),
                ct.byref(E),
                ct.byref(ct.c_double(R)),
                ct.byref(ct.c_double(Ma)),
            )

            return A.value, T.value, E.value

    def full(self, F, M, X, Z, C, R, Ma):
        """
        Runs settle, needs the full set of parameters.

        Can pass either scalars, or equally long arrays.
        """

        T = ct.c_double()
        A = ct.c_double()
        E = ct.c_double()

        if hasattr(M, "__iter__"):

            resA = []
            resR = []
            resE = []

            for i in range(len(M)):

                ret = self.mainer(
                    ct.byref(ct.c_double(F[i])),
                    ct.byref(ct.c_double(Z[i])),
                    ct.byref(ct.c_double(X[i])),
                    ct.byref(ct.c_double(M[i])),
                    ct.byref(ct.c_double(C[i])),
                    ct.byref(T),
                    ct.byref(A),
                    ct.byref(E),
                    ct.byref(ct.c_double(R[i])),
                    ct.byref(ct.c_double(Ma[i])),
                )

                resA.append(A.value)
                resR.append(T.value)
                resE.append(E.value)

            return numpy.array(resA), numpy.array(resR), numpy.array(resE)

        else:
            #print(self.mainer(ct.c_double(F),ct.c_double(Z),ct.c_double(X),ct.c_double(M),ct.c_double(C),ct.c_double(T),ct.c_double(A),ct.c_double(E),ct.c_double(R),ct.c_double(Ma)))
            ret = self.mainer(
                ct.byref(ct.c_double(F)),
                ct.byref(ct.c_double(Z)),
                ct.byref(ct.c_double(X)),
                ct.byref(ct.c_double(M)),
                ct.byref(ct.c_int(C)),
                ct.byref(T),
                ct.byref(A),
                ct.byref(E),
                ct.byref(ct.c_double(R)),
                ct.byref(ct.c_double(Ma)),
            )

            return A.value, T.value, E.value
