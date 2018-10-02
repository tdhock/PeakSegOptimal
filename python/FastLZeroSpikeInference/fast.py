import numpy as np
import ctypes as ct
import os
from utils import update_path_stats, get_num_changepts, get_cost
import warnings
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


def estimate_spikes(dat, gam, penalty, constraint=False, estimate_calcium=False, EPS=1e-10):
    """
    Estimate spike train, underlying calcium concentration, and changepoints based on fluorescence
    trace.

    Args:

    dat (np.array): fluorescence data
    gam (double): a scalar value for the AR(1) decay parameter; 0 < gam <= 1
    penalty (double): tuning parameter lambda > 0
    constraint (bool): constrained (true) or unconstrained (false) optimization
    estimate_calcium (bool): boolean specifying whether to estimate the calcium
    EPS (double): minimum calcium value (>= 0)

    Returns:

    dict with elements

    spikes: the set of spikes
    estimated_calcium: estimated calcium concentration
    change_pts: the set of changepoints
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

    We introduce the constant EPS > 0, typically on the order of 10^-10, to avoid
    arbitrarily small calcium concentrations that would result in numerical
    instabilities. In practice, this means that the estimated calcium concentration
    decays according to the AR(1) model for values greater than EPS and is equal to EPS thereafter.

    When estimating the spikes, it is not necessary to explicitly compute the
    calcium concentration. Therefore, if only the spike times are required, the user
    can avoid this computation cost by setting the compute_fitted_values boolean to false.
    By default, the calcium concentration is not estimated.

    Given the set of estimated spikes produced from the estimate_spike, the calcium concentration
    can be estimated with the estimate_calcium function.

    For additional information see:
    1. Jewell, Hocking, Fearnhead, and Witten (2018) <arXiv:1802.07380> and
    2. Jewell and Witten (2017) <arXiv:1703.08644>

    Examples:

    import numpy as np

    gam = 0.98
    y = np.power(gam, np.concatenate([np.arange(100), np.arange(100)]))

    # Basic fit, no calcium concentration estimated
    fit = fast.estimate_spikes(y, gam, 1, False)

    # Determine calcium concentration from fit
    fit = fast.estimate_calcium(fit)

    # Estimate both spikes and calcium concentration
    fit = fast.estimate_spikes(y, gam, 1, False, True)


    """
    if (gam > 1 or gam <= 0):
        print "gam must be in interval (0, 1]"
        return 0
    if (dat.shape[0] < 3):
        print "number of data points must be >2"
        return 0
    if (EPS < 0):
        print "EPS must be >= 0"
        return 0

    dat = np.ascontiguousarray(dat, dtype=float)
    dat_count = dat.shape[0]
    cost_mat = np.ascontiguousarray(np.zeros(dat_count, dtype=float))
    end_vec = np.ascontiguousarray(np.zeros(dat_count, dtype=np.int32))
    estimated_calcium = np.ascontiguousarray(np.zeros(dat_count, dtype=float))
    intervals_mat = np.ascontiguousarray(np.zeros(dat_count, dtype=np.int32))
    success = False

    lib.ARFPOP_interface(dat.ctypes.data_as(ct.POINTER(ct.c_double)),  # data ptr
                         ct.pointer(ct.c_int(dat_count)),  # data count
                         ct.pointer(ct.c_double(penalty)),  # penalty
                         ct.pointer(ct.c_double(gam)),  # gamma
                         cost_mat.ctypes.data_as(ct.POINTER(ct.c_double)),  # cost mat
                         end_vec.ctypes.data_as(ct.POINTER(ct.c_int)),  # end_vec
                         estimated_calcium.ctypes.data_as(ct.POINTER(ct.c_double)),  # estimated_calcium
                         intervals_mat.ctypes.data_as(ct.POINTER(ct.c_int)),  # int_vec
                         ct.pointer(ct.c_bool(constraint)),  # constraint
                         ct.pointer(ct.c_int(success)),
                         ct.pointer(ct.c_bool(estimate_calcium)),  # fitted values
                         ct.pointer(ct.c_double(EPS)))

    out = {}
    out['change_pts'] = np.unique(end_vec) + 1
    out['spikes'] = out['change_pts'][1:] + 1
    out['end_vec'] = end_vec + 1
    out['dat'] = dat
    out['gam'] = gam
    out['EPS'] = EPS
    out['penalty'] = penalty
    out['cost'] = cost_mat

    if (estimate_calcium):
        out['estimated_calcium'] = estimated_calcium
        out['spike_mag'] = out['estimated_calcium'][out['spikes'] - 1] - gam * out['estimated_calcium'][
            out['spikes'] - 2]
        out['pos_spike_mag'] = np.maximum(out['spike_mag'], 0)
    else:
        out['estimated_calcium'] = None
        out['spike_mag'] = None
        out['pos_spike_mag'] = None

    ## TODO: add error catching here!

    return out


def estimate_calcium(fit):
    """
    Estimate underlying calcium concentration based on fit from estimate_spikes

    Args:

    fit: object returned from running estimate_spikes

    Returns:

    dict with elements

    spikes: the set of spikes
    estimated_calcium: estimated calcium concentration
    change_pts: the set of changepoints
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

    We introduce the constant EPS > 0, typically on the order of 10^-10, to avoid
    arbitrarily small calcium concentrations that would result in numerical
    instabilities. In practice, this means that the estimated calcium concentration
    decays according to the AR(1) model for values greater than EPS and is equal to EPS thereafter.

    When estimating the spikes, it is not necessary to explicitly compute the
    calcium concentration. Therefore, if only the spike times are required, the user
    can avoid this computation cost by setting the compute_fitted_values boolean to false.
    By default, the calcium concentration is not estimated.

    Given the set of estimated spikes produced from the estimate_spike, the calcium concentration
    can be estimated with the estimate_calcium function.

    For additional information see:
    1. Jewell, Hocking, Fearnhead, and Witten (2018) <arXiv:1802.07380> and
    2. Jewell and Witten (2017) <arXiv:1703.08644>

    Examples:

    import numpy as np

    gam = 0.98
    y = np.power(gam, np.concatenate([np.arange(100), np.arange(100)]))

    # Basic fit, no calcium concentration estimated
    fit = fast.estimate_spikes(y, gam, 1, False)

    # Determine calcium concentration from fit
    fit = fast.estimate_calcium(fit)

    # Estimate both spikes and calcium concentration
    fit = fast.estimate_spikes(y, gam, 1, False, True)



    """
    dat = np.ascontiguousarray(fit['dat'], dtype=float)
    dat_count = fit['dat'].shape[0]
    end_vec = np.ascontiguousarray(np.array(fit['end_vec'] - 1, dtype=np.int32))
    estimated_calcium = np.ascontiguousarray(np.zeros(dat_count, dtype=float))

    lib.FitSegmentModel_interface(dat.ctypes.data_as(ct.POINTER(ct.c_double)),  # data ptr
                                  ct.pointer(ct.c_int(dat_count)),  # data count
                                  ct.pointer(ct.c_double(fit['gam'])),  # gamma
                                  end_vec.ctypes.data_as(ct.POINTER(ct.c_int)),  # end_vec
                                  estimated_calcium.ctypes.data_as(ct.POINTER(ct.c_double)),  # estimated_calcium
                                  ct.pointer(ct.c_double(fit['EPS'])))

    out = fit

    out['estimated_calcium'] = estimated_calcium
    out['spike_mag'] = out['estimated_calcium'][out['spikes'] - 1] - out['gam'] * out['estimated_calcium'][
        out['spikes'] - 2]
    out['pos_spike_mag'] = np.maximum(out['spike_mag'], 0)

    ## TODO: add error catching here!

    return out


def estimate_spike_paths(dat, gam, lambda_min=1e-2, lambda_max=1e1, constraint=False, EPS=1e-10, max_iters=10):
    lambdas_used = (lambda_min, lambda_max)
    path_fits = []
    path_stats = []
    approximate_path = False

    # 1. Run CPD for penalty values lambda_min and lambda_max;
    path_fits.append(estimate_spikes(dat, gam, lambda_min, constraint, EPS))
    path_stats = update_path_stats(path_stats, path_fits[0])
    path_fits.append(estimate_spikes(dat, gam, lambda_max, constraint, EPS))
    path_stats = update_path_stats(path_stats, path_fits[1])

    n_fits = 1

    # 2. Set lambda_star = {[lambda_min, lambda_max]};
    lambda_star = [lambdas_used]

    while (len(lambda_star) > 0 and n_fits < max_iters):
        # 3. Choose an element of lambda_star; denote this element as [lambda_0,lambda_1];
        # here always take the biggest interval

        interval_sizes = [interval[1] - interval[0] for interval in lambda_star]
        max_interval = np.argsort(interval_sizes)[-1]

        current_interval = lambda_star[max_interval]
        if (get_num_changepts(current_interval[0], path_stats) > (get_num_changepts(current_interval[1], path_stats) + 1)):
            lambda_int = (get_cost(current_interval[1], path_stats) - get_cost(current_interval[0], path_stats)) / \
                         (get_num_changepts(current_interval[0], path_stats) - get_num_changepts(current_interval[1],
                                                                                                path_stats))

            if (lambda_int in path_stats['penalty']):
                print(lambda_int)

            n_fits = n_fits + 1
            path_fits.append(estimate_spikes(dat, gam, lambda_int, constraint, EPS))
            path_stats = update_path_stats(path_stats, path_fits[n_fits])

            if (get_num_changepts(lambda_int, path_stats) != get_num_changepts(current_interval[0], path_stats)):
                # Set lambda_star = {lambda_star,[lambda_0,lambda_int),[lambda_int,lambda_1]}.;
                lambda_star.append((current_interval[0], lambda_int))
                lambda_star.append((lambda_int, current_interval[1]))

        del lambda_star[max_interval]

        if n_fits == max_iters:
            approximate_path = True
            warnings.warn("Full search path terminated early since maximum number of iterations (%d) reached. "
                    "Rerun with larger 'max_iter' parameter for full path." % max_iters)


    out = {'path_stats': path_stats, 'path_fits': path_fits, 'approximate_path': approximate_path}

    return (out)
