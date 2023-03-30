""" test to check if settle has been compiled correctly """

from pySettle  import settler as se
import numpy as np

def test_settle_location():
    se.Settle()
    return

def test_settle_output():
    # initialize settle interface
    settle = se.Settle()

    # run settle:
    res = settle.full(F=0.1,M=0.1,X=0.7,Z=0.02, C=0, R=11.2, Ma=1.4)
    result = np.allclose(res, [66.32432920153866, 4.630885096736736, 7.516383459074593]) #alpha,tdel, E_b
    if result:
        print("PASSED")
    else:
        print("FAILED")
    assert result

    return


test_settle_location()
test_settle_output()


