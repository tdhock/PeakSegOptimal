import numpy as np
import ctypes as ct 
import os

import sys

sys_path = sys.path
success = False

for path in sys_path:
    try:
        lib = np.ctypeslib.load_library('FastLZeroSpikeInference', path)
        success = True
        break
    except OSError as err:
        if '{0}'.format(err) == 'no file with expected extension':
            e = err
            continue
        else:
            print("OSError: {0}".format(err))
            raise
    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise

assert success, "Failed to import 'FastZeroSpikeInference': {}".format(e)



def arfpop(dat, gam, penalty, constraint):
	dat = np.ascontiguousarray(dat, dtype = float)
	dat_count = dat.shape[0]
	cost_mat = np.ascontiguousarray(np.zeros(dat_count, dtype=float))
	end_vec = np.ascontiguousarray(np.zeros(dat_count, dtype=np.int32))
	mean_vec = np.ascontiguousarray(np.zeros(dat_count, dtype=float))
	intervals_mat = np.ascontiguousarray(np.zeros(dat_count, dtype=np.int32))
	constraint = constraint
	success = False

	lib.ARFPOP_interface(dat.ctypes.data_as(ct.POINTER(ct.c_double)), # data ptr
	                                ct.pointer(ct.c_int(dat_count)),  # data count
	                                ct.pointer(ct.c_double(penalty)), # penalty
	                                ct.pointer(ct.c_double(gam)),  # gamma
	                                cost_mat.ctypes.data_as(ct.POINTER(ct.c_double)), # cost mat
	                                end_vec.ctypes.data_as(ct.POINTER(ct.c_int)), # end_vec
	                                mean_vec.ctypes.data_as(ct.POINTER(ct.c_double)), # mean_vec
	                                intervals_mat.ctypes.data_as(ct.POINTER(ct.c_int)), # int_vec
	                                ct.pointer(ct.c_bool(constraint)), # constraint
	                                ct.pointer(ct.c_bool(success)))

	out = {}
	out['mean_vec'] = np.flip(mean_vec, 0)
	out['intervals_mat'] = intervals_mat
	out['changePts'] = np.unique(end_vec) + 1
	out['spikes'] = out['changePts'][1:] + 1
	
	padded = np.array([0])

	out['spike_mag'] = np.concatenate((padded, out['mean_vec'][1:] - out['mean_vec'][0:-1]))
	out['pos_spike_mag'] = np.maximum(out['spike_mag'], np.zeros(out['spike_mag'].shape))

	## TODO: add error catching here! 

	return out 
