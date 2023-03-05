# This is Short Functional test of settle.
# settle calls are done via settler.py interface module.
#
# Test type: integration, end-to-end, positive.
# Runs a limited number of settle solutions and check result
# for the last one only agains expected values.

from __future__ import print_function
from pySettle import settler as se
import numpy as np

print("Settle Short Functional Test (SFT)...")

try:
    # initialize settle interface
    settl = se.Settle()

    for i in range(10):
        print(settl.full(Z=0.02, X=0.5, M=0.1, F=0.1, C=0, R=11.2, Ma=1.4))
        print(settl.run(Z=0.02, X=0.5, M=0.1, R=11.2, Ma=1.4))

    set2 = se.Settle(F=0.5)
    # having created the link with different default Flux will generate different output with the same compact call
    print(settl.run(Z=0.02, X=0.5, M=0.1, R=11.2, Ma=1.4), set2.run(Z=0.02, X=0.5, M=0.1, R=11.2, Ma=1.4))
    print()

    # you can override the default Flux you created using the full call, this will generate the same output,
    # even if you created different links with different defaults for Flux (same for C)
    print(settl.full(F=0.7, C=0, Z=0.02, X=0.5, M=0.1, R=11.2, Ma=1.4), set2.full(F=0.7, C=0, Z=0.02, X=0.5, M=0.1, R=11.2, Ma=1.4))
    print()

    print(settl.full(Z=0.026, X=0.505, M=0.1, F=0.12, C=0, R=11.16, Ma=1.426))

    print("settl.full() test - check expected result.")
    res = settl.full(F=0.1, M=0.1, X=0.7, Z=0.02, C=0, R=11.2, Ma=1.4)
    print(res)
    
    result = np.allclose(res, [66.32432920153866, 4.630885096736736, 7.516383459074593])
    # alpha,tdel, E_b
except:
    raise AssertionError("Settle SFT failed - some error occurred!")

if result:
    print("PASSED")
else:
    print("FAILED")
    raise AssertionError("Settle SFT failed - result of settler.full() not as expected!")

# this is just a simple check, without any useful infor printed.
# assert result
