#' Estimate spike train, underlying calcium concentration, and changepoints
#' based on a fluorescence trace.
#'
#' @param dat fluorescence data
#' @param gam a scalar value for the AR(1) decay parameter
#' @param lambda tuning parameter lambda
#' @param constraint boolean specifying constrained or unconstrained
#'   optimization problem (see below)
#' @param estimate_calcium boolean specifying whether to estimate the calcium
#' @param EPS double specifying the minimum calcium value
#'
#' @return Returns a list with elements:
#' @return \code{spikes} the set of estimated spikes
#' @return \code{estimated_calcium} estimated calcium concentration
#' @return \code{change_pts} the set of changepoints
#' @return \code{cost} the cost at each time point
#' @return \code{n_intervals} the number of intervals at each point
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

estimate_spikes <- structure(function(dat, gam, lambda, constraint = FALSE, estimate_calcium = FALSE, EPS = 1e-4) {
  stopifnot(gam > 0 && gam <= 1)
  stopifnot(EPS > 0)
  stopifnot(!is.null(gam))
  weight.vec <- rep(1, length(dat))
  n.data <- length(dat)
  stopifnot(3 <= n.data)
  stopifnot(is.numeric(lambda))
  stopifnot(length(lambda) == 1)
  stopifnot(0 <= lambda)
  cost.mat <- double(n.data)
  ends.vec <- integer(n.data)
  mean.vec <- double(n.data)
  intervals.mat <- integer(n.data)
  success <- 1
  result.list <- .C(
    "ARFPOP_interface",
    dat.vec = as.numeric(dat),
    n.data = as.integer(n.data),
    penalty = as.numeric(lambda),
    gam = as.numeric(gam),
    cost.mat = as.double(cost.mat),
    ends.vec = as.integer(ends.vec),
    mean.vec = as.double(mean.vec),
    intervals.mat = as.integer(intervals.mat),
    constraint = constraint,
    success = as.integer(success),
    compute_fitted_values = estimate_calcium, 
    EPS = as.numeric(EPS), 
    PACKAGE = "FastLZeroSpikeInference"
  )
  
  ## 1-indexed segment ends!
  result.list$ends.vec <-
    result.list$ends.vec + 1L
  result.list$cost.mat <- matrix(result.list$cost.mat, 1, n.data, byrow = TRUE)
  result.list$intervals.mat <- 
    matrix(result.list$intervals.mat, 1, n.data, byrow = TRUE)
  
  changePts <- rev(unique(result.list$ends.vec))
  spikes <- changePts[-1] + 1
  
  constraint_str <- ""
  if (constraint) {
    constraint_str <- "-pos-constrained"  
  }
  
  if (estimate_calcium) {
    fitted_values = result.list$mean.vec  
  } else {
    fitted_values = NA 
  }
  
  
  if (result.list$success) {
    out <-
      list(
        spikes = spikes,
        estimated_calcium = fitted_values,
        dat = dat,
        type = paste0("ar1", constraint_str), 
        change_pts = changePts,
        call = match.call(),
        gam = gam,
        lambda = lambda,
        cost = as.numeric(result.list$cost.mat),
        n_intervals = as.numeric(result.list$intervals.mat),
        end_vec = result.list$ends.vec,
        EPS = EPS
      )
    class(out) <- "estimated_spikes"
    return(out)
  } else {
    stop("Numerical issues!")
  }

  
})