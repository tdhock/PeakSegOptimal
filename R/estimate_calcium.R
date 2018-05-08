#' Estimate underlying calcium concentration based on estimated spikes
#' 
#' @param fit object created by running estimate_spikes
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
#'  \strong{AR(1) model:}
#'  
#'  minimize_{c1,...,cT} 0.5 sum_{t=1}^T ( y_t - c_t )^2 + lambda sum_{t=2}^T 1_[c_t != max(gam c_{t-1}, EPS)]
#'  
#'  for the global optimum, where y_t is the observed fluorescence at the tth timepoint.
#'
#' \strong{Constrained AR(1) model:}
#' 
#' minimize_{c1,...,cT} 0.5 sum_{t=1}^T ( y_t - c_t )^2 + lambda sum_{t=2}^T 1_[c_t != max(gam c_{t-1}, EPS)]
#' 
#' subject to c_t >= max(gam c_{t-1}, EPS), t = 2, ..., T
#'
#' We introduce the constant EPS > 0, typically on the order of 10^-10, to avoid 
#' arbitrarily small calcium concentrations that would result in numerical  
#' instabilities. In practice, this means that the estimated calcium concentration 
#' decays according to the AR(1) model for values greater than EPS and is equal to EPS thereafter.
#'
#' When estimating the spikes, it is not necessary to explicitly compute the 
#' calcium concentration. Therefore, if only the spike times are required, the user
#' can avoid this computation cost by setting the estimate_calcium boolean to false. 
#' Because estimating the calcium requires additional computation time, 
#' we suggest estimating the calcium only if it is needed.
#'
#' Given the set of estimated spikes produced from the estimate_spike, the calcium concentration
#' can be estimated with the estimate_calcium function (see examples below).
#'
#' For additional information see: 
#' 
#' 1. Jewell, Hocking, Fearnhead, and Witten (2018) <arXiv:1802.07380> and 
#' 
#' 2. Jewell and Witten (2017) <arXiv:1703.08644> 
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
#' fit <- estimate_spikes(dat = sim$fl, gam = 0.95, lambda = 1, estimate_calcium = T)
#' plot(fit)
#' 
#' # Constrained AR(1) model
#' fit <- estimate_spikes(dat = sim$fl, gam = 0.95, lambda = 1, constraint = T, estimate_calcium = T)
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
estimate_calcium <- structure(function(fit) {
  
  n.data <- length(fit$dat)
  ends.vec <- fit$end_vec - 1
  mean.vec <- double(n.data)
  result.list <- .C(
    "FitSegmentModel_interface",
    dat.vec = as.numeric(fit$dat),
    n.data = as.integer(n.data),
    gam = as.numeric(fit$gam),
    ends.vec = as.integer(ends.vec),
    mean.vec = as.double(mean.vec),
    EPS = as.numeric(fit$EPS), 
    PACKAGE = "FastLZeroSpikeInference"
  )
  
  out <-
    list(
      spikes = fit$spikes,
      estimated_calcium = result.list$mean.vec,
      dat = fit$dat,
      type = fit$type, 
      change_pts = fit$change_pts,
      call = fit$call,
      gam = fit$gam,
      lambda = fit$lambda,
      cost = fit$cost,
      n_intervals = fit$n_intervals,
      end_vec = fit$end_vec,
      EPS = fit$EPS
    )
  class(out) <- "estimated_spikes"
  return(out)
})