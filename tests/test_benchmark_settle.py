#!/usr/bin/env python
"""
Settle benchmark test

Is part of beansp tests, because uses settle.py wrapper,
which is part of beansp package

It is also demonstartion how pySettle can be parallelised,
providing the settling is not strictly iterative, but can be
run by chunks.

The speedup is a lot better than just dividing run time
by number of CPU cores, beacause of reducing python
overhead manipulating with data and/or initialising
the extension (external settle libarry coded in C/C++)
"""

import time
import multiprocessing as mp
import os

import numpy as np
# from astropy.io using table andascii
import astropy
from astropy.io import ascii
# from astropy.table.table_helpers import simple_table

# local modules
from beansp import settle
# imported like this it will not be identifiable where grabbed from
# from settle import settle


def settle_multiprocessing_wrapper(ft_list_item):
    # print("@")
    return settle.settle(ft_list_item['Q_b'],
                         ft_list_item['Z'],
                         ft_list_item['X'],
                         ft_list_item['mdot']/(1+ft_list_item['X']),
                         1.0,
                         ft_list_item['M_NS'],
                         ft_list_item['R_NS'])


def test_settle_benchmark():
    # MCu: figuring out What is imported and from where
    # print("sys.path=")
    # print(sys.path)
    print("--------------------")
    # MCu note - this olways works (well so far...)
    print("settle.__name__=", settle.__name__)

    # MCu note - the following one does not work for a sub-module
    if hasattr(settle, '__file__'):
        print("settle.__file__=", settle.__file__)
    else:
        print("settle does not have attribute __path__")

    # MCu note - the following one does not work for a module
    if hasattr(settle, '__path__'):
        print("settle.__path__=", settle.__path__)
    else:
        print("settle does not have attribute __path__")

    # using local relative data directory, as input testdata is part of the repo
    local_data_dir = 'data'
    data_dir = os.path.join(os.path.dirname(
        os.path.abspath(__file__)), local_data_dir)

    # just print newline chatacter
    print(" ")

    # ----- Comparison of settle with Kepler

    # Here we want to use the concord table and run a settle model
    # with each set to compare the results.
    # First need to read in the table

    # ft = astropy.io.ascii.read('benchmark/table1.mrt')
    ft = ascii.read(os.path.join(data_dir, 'benchmark/table1.mrt'))
    print("ft table rows = ", len(ft))
    if hasattr(ft, '__iter__'):
        print("astropy.table.Table ft is iterable")
    else:
        print("astropy.table.Table ft is NOT iterable")

    # print(ft)

    # add your path to the mrt file below
    # fixed parameters
    M_NS, R_NS, Q_b = 1.4, 10., 0.1

    # print (M_NS, R_NS)
    # for run in ft:
    #     print (run['run'], run['mdot'], run['X'], run['Z'])

    # This section was run on xray, where settle works

    tdel, E_b, alpha = [], [], []

    num_runs = 5
    print("num_runs = ", num_runs)

    print("======== serial run ==========")

    ts_start = time.process_time()
    ts1_sum = 0.0
    ts2_sum = 0.0

    serial_count_settle_runs = 0

    # ----- serial settle run loop ------
    # the inner cycle is replaced by parallelisation
    # in the parallel version

    for j in range(num_runs):
        ts2_start = time.process_time()
        for i, run in enumerate(ft['run']):
            serial_count_settle_runs += 1
            # print ('Running settle for run #{}...'.format(run))
            # need to convert the mdot here, I think this is right
            # In the MRT file accretion rate is given as a fraction of the Eddington rate, i.e.
            # Mdot_Edd = 8.8e4/(1+X) g/cm^2/s; and since settle uses fraction of 8.8e4, we have
            # an extra factor of (1+X) in the MRT values that we need to divide by
            ts1_start = time.process_time()
            # res = settl.full(Q_b, ft[i]['Z'], ft[i]['X'], ft[i]['mdot']/(1+ft[i]['X']), 1, R_NS, M_NS)
            res = settle.settle(Q_b, ft[i]['Z'],
                                ft[i]['X'],
                                ft[i]['mdot']/(1+ft[i]['X']),
                                1.0, M_NS, R_NS)
            ts1_end = time.process_time()
            ts1_sum += (ts1_end-ts1_start)
            tdel.append(res['tdel'][0])
            E_b.append(res['E_b'][0])
            alpha.append(res['alpha'][0])
        ts2_end = time.process_time()
        ts2_sum += (ts2_end-ts2_start)
        print("Cycle #", j, " time = ", (ts2_end-ts2_start))

    ts_end = time.process_time()

    print("serial_count_settle_runs =", serial_count_settle_runs)
    print("total process time (", num_runs, "loops) = ", ts_end - ts_start)
    print("settle sum time (", num_runs, "loops) = ", ts1_sum)
    print("loop sum time (", num_runs, "loops) = ", ts2_sum)
    print("average", len(ft),
          "row data table one loop run settle sum time = ", (ts1_sum/num_runs))

    ts = astropy.table.Table([ft['run'],
                              tdel[: len(ft)],
                              E_b[: len(ft)],
                              alpha[: len(ft)]])
    # dump results produced by serial form of settle run
    ts.write(os.path.join(data_dir, 'serial.ecsv'),
             format='ascii.ecsv', overwrite=True)

    print("======== end of serial run ==========")

    # ------ clear lists ----------

    tdel.clear()
    E_b.clear()
    alpha.clear()

    print("======== multiprocessing parallel run ==========")

    cpu_count = mp.cpu_count()

    print("Detected CPUs: ", cpu_count)

    chunksize = round(len(ft)/cpu_count)

    print("Chunk size: ", chunksize)

    t_start = time.process_time()

    t2_sum = 0.0

    # ----- parallel settle run loop ------

    for j in range(num_runs):
        t2_start = time.process_time()
        # add input parameters into the array
        # yes, there is a ridiculous redundancy,
        # but multiprocessing.Pool.map() takes only one parameter
        ft['Q_b'] = Q_b
        ft['M_NS'] = M_NS
        ft['R_NS'] = R_NS
        # need to convert astropy.table.Table to list (of dictionaries)
        # because table is not iterable in Python >:-(
        ft_list = [dict(zip(ft.colnames, row)) for row in ft]
        # print(ft_list)

        # Simply configured like this, the speedup is excellent
        # the magic is done mostly by removing some internal Python
        # memory management inefficiency
        pool = mp.Pool(cpu_count)

        # The following is pure speed-up only by CPUs, removing the
        # beneficial effect of bypassing some internal Python
        # memory management inefficiency
        # pool = mp.Pool(1, maxtasksperchild=1)

        multi_res = pool.map(settle_multiprocessing_wrapper,
                             ft_list, chunksize=chunksize)

        # free up the resources used by process pool
        # WTF dunno why I have to call this before being
        # able to wait for processes to join!
        # (makes very little logical sense...)
        pool.close()
        # wait for all print()arallel processes to finish
        # they are identical, just with different
        # numerical (float) parameters, so shall take about same
        # time to execute, providing the selected chunk size make sense
        pool.join()
        t2_end = time.process_time()
        t2_sum += (t2_end-t2_start)
        print("Cycle #", j, " time = ", (t2_end-t2_start))

    # err print(multi_res['tdel'])

    # ----- unwrapp the typeless blob .......
    # this is more complicated than the parallelisation!

    # just debug prints
    print("len(ft)=", len(ft))
    print(type(multi_res))
    print(type(multi_res[0]))
    print(type(multi_res[0: len(ft)][0]))

    # it is impossible to subscript a list -> convert to np.array
    multi_res_array = np.array(multi_res)
    tdel_wrap1_vector = multi_res_array['tdel']
    E_b_wrap1_vector = multi_res_array['E_b']
    alpha_wrap1_vector = multi_res_array['alpha']

    # print("tdel_wrap1_vector =", tdel_wrap1_vector)
    print("type(tdel_wrap1_vector) =", type(tdel_wrap1_vector))
    # print("tdel_wrap1_vector[0] =", tdel_wrap1_vector[0])
    print("type(tdel_wrap1_vector[0]) =", type(tdel_wrap1_vector[0]))

    # this does not do nothing
    # tdel_float_vector = tdel_vector.astype(float)
    
    # tdel_wrap2_vector = tdel_wrap1_vector[0]
    # E_b_wrap2_vector = E_b_wrap1_vector[0]
    # alpha_wrap2_vector = alpha_wrap1_vector[0]

    # print("tdel_wrap2_vector =", tdel_wrap2_vector)
    # print("type(tdel_wrap2_vector) =", type(tdel_wrap2_vector))
    # print("tdel_wrap2_vector[0] =", tdel_wrap2_vector[0])
    # print("type(tdel_wrap2_vector[0]) =", type(tdel_wrap2_vector[0]))

    tdel_float64_vector = np.array([tdel_float64[0] for tdel_float64 in tdel_wrap1_vector])
    E_b_float64_vector = np.array([E_b_float64[0] for E_b_float64 in E_b_wrap1_vector])
    alpha_float64_vector = np.array([alpha_float64[0] for alpha_float64 in alpha_wrap1_vector])

    # does not work :-(
    # for i in range(len(ft)):
    #    tdel_float_vector[i]=float(tdel_string_vector[i])

    t_end = time.process_time()
    print("total process time (", num_runs, "loops) = ", t_end - t_start)
    print("loop sum time (", num_runs, "loops) = ", t2_sum)
    print("average", len(ft),
          " row data table one loop average time =", t2_sum/num_runs)
    t = astropy.table.Table([ft['run'],
                             tdel_float64_vector,
                             E_b_float64_vector,
                             alpha_float64_vector])

    print(type(t))
    # dump results produced by serial form of settle run
    t.write(os.path.join(data_dir, 'multiprocessing.ecsv'),
            format='ascii.ecsv', overwrite=True)

    print("======== end of multiprocessing parallel run ==========")

    print("Speed up is :", (ts_end - ts_start)/(t_end - t_start))

    print("======== compare files if result as expected ==========")

    ref = astropy.table.Table.read(
        os.path.join(data_dir, 'benchmark/settle_expected_result.ecsv'),
        format='ascii.ecsv')

    result = False
    for col in ref.colnames:
        result = np.allclose(ref[col], ts[col])
        if not result:
            print("serial run column " + col + " not identical!")
            break
        result = np.allclose(ref[col], t[col])
        if not result:
            print("multiprocessing run column " + col + " not identical!")
            break

    if result:
        print("PASSED")
    else:
        print("FAILED")
        raise AssertionError("Settle benchmark test failed - \
        result not as expected!")


if __name__ == "__main__":
    # introspect and run all the functions starting with 'test'
    for f in dir():
        if f.startswith('test'):
            print(f)
            globals()[f]()
