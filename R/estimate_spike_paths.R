# Created by: jewellsean
# Created on: 2018-09-26

## implement crops algo
## https://www.tandfonline.com/doi/abs/10.1080/10618600.2015.1116445

# CROPS algorithm
# input: A dataset y1:n = (y1,y2,...,yn);
# Minimum and maximum values of the penalty,
# lambda_min and lambda_max;
# output: The details of optimal segmentations

## get arfpop stats
## input: arfpop object
## output: dataframe with cols (lambda, m(lambda), Q_m(lambda)(y1:n))
arfpop_stats <- function(fit) {
    df <- data.frame(lambda = fit$lambda,
    changepoints_n = length(fit$spikes),
    cost = fit$cost[length(fit$cost)] - length(fit$spikes) * fit$lambda)
    return(df)
}

update_path_stats <- function(df, new_fit) {
    return(rbind(df, arfpop_stats(new_fit)))
}

get_num_changepts <- function(lambda, df) {
    df$changepoints_n[df$lambda == lambda][1]
}

get_cost <- function(lambda, df) {
    df$cost[df$lambda == lambda][1]
}


#' Estimate spike train, underlying calcium concentration, and changepoints
#' based on a fluorescence trace. Automatic tuning parameter selection within
#' a range of values [lambda_min, lambda_max]
#'
#' @param dat fluorescence data
#' @param gam a scalar value for the AR(1) decay parameter
#' @param lambda_min minimum lamba value
#' @param lambda_max maximum lamba value
#' @param constraint boolean specifying constrained or unconstrained
#'   optimization problem (see below)
#' @param EPS double specifying the minimum calcium value
#' @param max_iters maximum number of tuning parameters attempted
#'
#' @return Returns a list with elements:
#' @return \code{path_stats} a dataframe with summary statistics (number of spikes, tuning parameters, cost)
#' @return \code{path_fits} a list with estimated_spikes object for each tuning parameter
#' @return \code{approximate_path} a boolean indicating whether an early stopping criterion condition occurred
#'
#' @details
#'
#' This algorithm solves the optimization problems
#'
#' \strong{AR(1) model:}
#'
#' minimize_{c1,...,cT} 0.5 sum_{t=1}^T ( y_t - c_t )^2 + lambda sum_{t=2}^T
#' 1_[c_t != max(gam c_{t-1}, EPS)]
#'
#' for the global optimum, where y_t is the observed fluorescence at the tth
#' timestep.
#'
#' \strong{Constrained AR(1) model:}
#'
#' minimize_{c1,...,cT} 0.5 sum_{t=1}^T ( y_t - c_t )^2 + lambda sum_{t=2}^T
#' 1_[c_t != max(gam c_{t-1}, EPS)]
#'
#' subject to c_t >= max(gam c_{t-1}, EPS), t = 2, ..., T
#'
#' We introduce the constant EPS > 0, to avoid
#' arbitrarily small calcium concentrations that would result in numerical
#' instabilities. In practice, this means that the estimated calcium
#' concentration decays according to the AR(1) model for values greater than EPS and
#' is equal to EPS thereafter.
#'
#' When estimating the spikes, it is not necessary to explicitly compute the
#' calcium concentration. Therefore, if only the spike times are required, the
#' user can avoid this computation cost by setting the estimate_calcium
#' boolean to false. Because estimating the calcium requires additional computation time,
#' we suggest estimating the calcium only if it is needed.
#'
#' Given the set of estimated spikes produced from the estimate_spike, the
#' calcium concentration can be estimated with the estimate_calcium function
#' (see examples below).
#'
#' For additional information see:
#'
#' 1. Jewell, Hocking, Fearnhead, and Witten (2018) <arXiv:1802.07380> and
#'
#' 2. Jewell, Sean; Witten, Daniela. Exact spike train inference via l0 optimization. 
#' Ann. Appl. Stat. 12 (2018), no. 4, 2457--2482. doi:10.1214/18-AOAS1162. 
#' https://projecteuclid.org/euclid.aoas/1542078052
#'
#' @examples
#'
#' sim <- simulate_ar1(n = 500, gam = 0.95, poisMean = 0.009, sd = 0.05, seed = 1)
#' plot(sim)
#'
#' ## Fits for tuning parameters between [0.1, 10]
#' fits <- estimate_spike_paths(dat = sim$fl, gam = 0.95, lambda_min = 0.1, lambda_max = 10)
#' print(fits)
#' plot(fits)
#' print(fits$path_fits[[1]])
#' plot(fits$path_fits[[1]])
#'
#' ## Fits for a single tuning parameter
#'
#' # AR(1) model
#' fit <- estimate_spikes(dat = sim$fl, gam = 0.95, lambda = 1)
#' print(fit)
#'
#' # compute fitted values from prev. fit
#' fit <- estimate_calcium(fit)
#' plot(fit)
#'
#' # or
#' fit <- estimate_spikes(dat = sim$fl, gam = 0.95, lambda = 1, estimate_calcium = TRUE)
#' plot(fit)
#'
#' # Constrained AR(1) model
#' fit <- estimate_spikes(dat = sim$fl, gam = 0.95, lambda = 1, constraint = TRUE,
#'                                                     estimate_calcium = TRUE)
#' print(fit)
#' plot(fit)
#'
#' @seealso
#' \strong{Estimate spikes:}
#' \code{\link{estimate_spikes}}
#' \code{\link{estimate_calcium}}
#'
#' \strong{Simulate:}
#' \code{\link{simulate_ar1}}
#'
#' @useDynLib FastLZeroSpikeInference
#' @export
estimate_spike_paths <- function(dat, gam, lambda_min = 1e-2, lambda_max = 1e1, constraint = FALSE, EPS = 1e-4, max_iters = 10) {
    lambdas_used <- c(lambda_min, lambda_max)
    path_fits <- list()
    path_stats <- NULL
    approximate_path <- FALSE

    ## input validity
    stopifnot(lambda_min > 0)
    stopifnot(lambda_max > 0)
    stopifnot(max_iters > 4)
    stopifnot(gam > 0 & gam < 1)
    stopifnot(is.numeric(dat))

    ## 1. Run CPD for penalty values lambda_min and lambda_max
    path_fits[[1]] <- estimate_spikes(dat, gam, lambda_min, constraint, EPS)
    path_stats <- update_path_stats(path_stats, path_fits[[1]])
    path_fits[[2]] <- estimate_spikes(dat, gam, lambda_max, constraint, EPS)
    path_stats <- update_path_stats(path_stats, path_fits[[2]])

    n_fits <- 2

    ## 2. Set lambda_star = {[lambda_min, lambda_max]}
    lambda_star <- list(lambdas_used)

    while (length(lambda_star) > 0 &&
        n_fits <= max_iters) {
        # 3. Choose an element of lambda_star; denote this element as [lambda_0,lambda_1];
        # here always take the first element of list
        max_interval_size <- 0
        for (interval_i in 1:length(lambda_star)) {
            interval = lambda_star[[interval_i]]
            if ((interval[2] - interval[1]) > max_interval_size) {
                ind <- interval_i
            }
        }

        current_interal <- lambda_star[[ind]]

        if (get_num_changepts(current_interal[1], path_stats) >
        get_num_changepts(current_interal[2], path_stats) + 1) {
            lambda_int <- (get_cost(current_interal[2], path_stats) -
            get_cost(current_interal[1], path_stats)) / (
            get_num_changepts(current_interal[1], path_stats) -
            get_num_changepts(current_interal[2], path_stats)
            )

            n_fits <- n_fits + 1
            path_fits[[n_fits]] <- estimate_spikes(dat, gam, lambda_int, constraint, EPS)
            path_stats <- update_path_stats(path_stats, path_fits[[n_fits]])

            if (get_num_changepts(lambda_int, path_stats) !=
            get_num_changepts(current_interal[2], path_stats)) {
                ## Set lambda_star = {lambda_star,[lambda_0,lambda_int),[lambda_int,lambda_1]}.;
                n_intervals <- length(lambda_star)
                lambda_star[[n_intervals + 1]] <- c(current_interal[1], lambda_int)
                lambda_star[[n_intervals + 2]] <- c(lambda_int, current_interal[2])
            }
        }
        lambda_star[[ind]] <- NULL

        if (n_fits == max_iters) {
            warning(paste0("Full search path terminated early since maximum number of iterations (", max_iters,
            ") reached. Rerun with larger 'max_iter' parameter for full path."))
            approximate_path = TRUE
        }
    }

    out <- list(path_stats = path_stats, path_fits = path_fits, approximate_path = approximate_path)
    class(out) <- "estimated_spike_paths"

    return(out)
}

