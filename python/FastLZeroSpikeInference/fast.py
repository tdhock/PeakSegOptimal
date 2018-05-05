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



def arfpop(dat, gam, penalty, constraint = False, compute_fitted_values = False, EPS = 1e-10):
	"""
	Estimate spike train, underlying calcium concentration, and changepoints based on fluorescence
	trace.
	
	Args:

	dat (np.array): fluorescence data
	gam (double): a scalar value for the AR(1) decay parameter; 0 < gam <= 1
	penalty (double): tuning parameter lambda > 0
	constraint (bool): constrained (true) or unconstrained (false) optimization
	compute_fitted_values (bool): specifying whether fitted values are calculated
	EPS (double): minimum calcium value
	
	Returns: 
	
	dict with elements 

	spikes: the set of spikes
	mean_vec: estimated calcium concentration
	changePts: the set of changepoints
	spike_mag: mag. of change to calcium concentration at a spike
	pos_spike_mag: positive subset of spike_mag

	as well as a few internal elements. 

	This algorithm solves the optimization problems
 	AR(1) model:
	  minimize_{c1,...,cT} 0.5 sum_{t=1}^T ( y_t - c_t )^2 + lambda sum_{t=2}^T 1_[c_t != max(gam c_{t-1}, EPS)]
	  
	  for the global optimum, where y_t is the observed fluorescence at the tth timepoint.
	
	Constrained AR(1) model:
		minimize_{c1,...,cT} 0.5 sum_{t=1}^T ( y_t - c_t )^2 + lambda sum_{t=2}^T 1_[c_t != max(gam c_{t-1}, EPS)]
		subject to			 c_t >= max(gam c_{t-1}, EPS), t = 2, ..., T

	For additional information see:  
	1. Jewell, Hocking, Fearnhead, and Witten (2018) <arXiv:1802.07380> and  
	2. Jewell and Witten (2017) <arXiv:1703.08644> 

	Examples: 

	import numpy as np

	gam = 0.98
	y = np.power(gam, np.concatenate([np.arange(100), np.arange(100)]))
	
	# Basic fit, no calcium concentration estimated
	fit = fast.arfpop(y, gam, 1, False)

	# Determine calcium concentration from fit
	fit = fast.fit_ar1_model(fit)

	# Estimate both spikes and calcium concentration
	fit = fast.arfpop(y, gam, 1, False, True)


	"""
	if (gam > 1 or gam <= 0):
		print "gam must be in interval (0, 1]"
		return 0
	if (dat.shape[0] < 3):
		print "number of data points must be >2"
		return 0

	dat = np.ascontiguousarray(dat, dtype = float)
	dat_count = dat.shape[0]
	cost_mat = np.ascontiguousarray(np.zeros(dat_count, dtype=float))
	end_vec = np.ascontiguousarray(np.zeros(dat_count, dtype=np.int32))
	mean_vec = np.ascontiguousarray(np.zeros(dat_count, dtype=float))
	intervals_mat = np.ascontiguousarray(np.zeros(dat_count, dtype=np.int32))
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
	                                ct.pointer(ct.c_int(success)),
	                                ct.pointer(ct.c_bool(compute_fitted_values)), # fitted values	
	                                ct.pointer(ct.c_double(EPS)))

	out = {}
	out['changePts'] = np.unique(end_vec) + 1
	out['spikes'] = out['changePts'][1:] + 1
	out['end_vec'] = end_vec + 1
	out['dat'] = dat
	out['gam'] = gam 
	out['EPS'] = EPS

	if (compute_fitted_values):
		out['mean_vec'] = mean_vec
		out['spike_mag'] = out['mean_vec'][out['spikes'] - 1] - gam * out['mean_vec'][out['spikes'] - 2]
		out['pos_spike_mag'] = np.maximum(out['spike_mag'], 0)
	else:
		out['mean_vec'] = None
		out['spike_mag'] = None
		out['pos_spike_mag'] = None

	## TODO: add error catching here! 

	return out 


def fit_ar1_model(fit):
	"""
	Estimate underlying calcium concentration based on fit from arfpop
	
	Args:

	fit: object returned from running arfpop
	
	Returns: 
	
	dict with elements 

	spikes: the set of spikes
	mean_vec: estimated calcium concentration
	changePts: the set of changepoints
	spike_mag: mag. of change to calcium concentration at a spike
	pos_spike_mag: positive subset of spike_mag

	as well as a few internal elements. 

	This algorithm solves the optimization problems
 	AR(1) model:
	  minimize_{c1,...,cT} 0.5 sum_{t=1}^T ( y_t - c_t )^2 + lambda sum_{t=2}^T 1_[c_t != max(gam c_{t-1}, EPS)]
	  
	  for the global optimum, where y_t is the observed fluorescence at the tth timepoint.
	
	Constrained AR(1) model:
		minimize_{c1,...,cT} 0.5 sum_{t=1}^T ( y_t - c_t )^2 + lambda sum_{t=2}^T 1_[c_t != max(gam c_{t-1}, EPS)]
		subject to			 c_t >= max(gam c_{t-1}, EPS), t = 2, ..., T

	For additional information see:  
	1. Jewell, Hocking, Fearnhead, and Witten (2018) <arXiv:1802.07380> and  
	2. Jewell and Witten (2017) <arXiv:1703.08644> 

	Examples: 

	import numpy as np

	gam = 0.98
	y = np.power(gam, np.concatenate([np.arange(100), np.arange(100)]))
	
	# Basic fit, no calcium concentration estimated
	fit = fast.arfpop(y, gam, 1, False)

	# Determine calcium concentration from fit
	fit = fast.fit_ar1_model(fit)

	# Estimate both spikes and calcium concentration
	fit = fast.arfpop(y, gam, 1, False, True)

	"""
	dat = np.ascontiguousarray(fit['dat'], dtype = float)
	dat_count = fit['dat'].shape[0]
	end_vec = np.ascontiguousarray(np.array(fit['end_vec'] - 1, dtype=np.int32))
	mean_vec = np.ascontiguousarray(np.zeros(dat_count, dtype=float))

	lib.FitSegmentModel_interface(dat.ctypes.data_as(ct.POINTER(ct.c_double)), # data ptr
	                                ct.pointer(ct.c_int(dat_count)),  # data count
	                                ct.pointer(ct.c_double(fit['gam'])),  # gamma
	                                end_vec.ctypes.data_as(ct.POINTER(ct.c_int)), # end_vec
	                                mean_vec.ctypes.data_as(ct.POINTER(ct.c_double)), # mean_vec
	                                ct.pointer(ct.c_double(fit['EPS'])))

	out = fit

	out['mean_vec'] = mean_vec
	out['spike_mag'] = out['mean_vec'][out['spikes'] - 1] - out['gam'] * out['mean_vec'][out['spikes'] - 2]
	out['pos_spike_mag'] = np.maximum(out['spike_mag'], 0)


	## TODO: add error catching here! 

	return out 



