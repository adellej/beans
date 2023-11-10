""" test to check if settle has been compiled correctly """

# from pySettle  import settler as se
import pySettle
import numpy as np

def test_settle_location():
    # se.Settle()
    pySettle.settler.Settle()
    return

def test_settle_output():
    settle_version = pySettle.__version__
    # initialize settle interface
    # settle = se.Settle()
    settle = pySettle.settler.Settle()

    # run settle:
    res = settle.full(F=0.1,M=0.1,X=0.7,Z=0.02, C=0, R=11.2, Ma=1.4)
    # updated here with the new (corrected) alpha calculation, for
    # pySettle v0.1.3 and later
    ## result = np.allclose(res, [66.32432920153866, 
    # result = np.allclose(res, [43.846559833302855, 
    #     4.630885096736736, 7.516383459074593]) #alpha,tdel, E_b
    # and now with the removal of the recurrence time scaling, and
    # scaling alpha by the redshift for the observer frame, for
    # v1.3.0 and later
    #alpha,tdel, E_b
    result = np.allclose(res, [55.209529360603085, 7.124438610364355, 7.516383459074704])

    if result:
        print("PASSED")
    else:
        print("FAILED; pySettle v{}, res={}".format(settle_version,res))
    assert result

    return


test_settle_location()
test_settle_output()


