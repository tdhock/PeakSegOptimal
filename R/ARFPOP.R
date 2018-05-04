#' Estimate spike train, underlying calcium concentration, and changepoints based on fluorescence
#' trace.
#'
#' @param dat fluorescence data
#' @param gam a scalar value for the AR(1) decay parameter.
#' @param lambda tuning parameter lambda
#' @param constraint boolean specifying constrained or unconstrained optimization
#' problem (see below)
#' @param compute_fitted_values boolean specifying whether fitted values are calculated
#' @param EPS double specfying the minimum calcium value
#'
#' @return Returns a list with elements:
#' @return \code{spikes} the set of spikes
#' @return \code{fittedValues} estimated calcium concentration
#' @return \code{changePts} the set of changepoints
#' @return \code{cost} the cost at each time point 
#' @return \code{nIntervals} the number of intervals at each point
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
#' c_t >= max(gam c_{t-1}, EPS), t = 2, ..., T
#'
#' For additional information see: 
#' 
#' 1. Jewell, Hocking, Fearnhead, and Witten (2018) <arXiv:1802.07380> and 
#' 
#' 2. Jewell and Witten (2017) <arXiv:1703.08644> 
#' 
#' @examples
#'
#' library(LZeroSpikeInference)
#' sim <- simulateAR1(n = 500, gam = 0.95, poisMean = 0.009, sd = 0.05, seed = 1)
#' plot(sim)
#'
#' ## Fits for a single tuning parameter
#'
#' # AR(1) model
#' fit <- ARFPOP(dat = sim$fl, gam = 0.95, lambda = 1)
#' print(fit)
#'
#' # compute fitted values from prev. fit
#' fit <- fit_ar1_model(fit)
#' plot(fit)
#'
#' # or 
#' fit <- ARFPOP(dat = sim$fl, gam = 0.95, lambda = 1, compute_fitted_values = T)
#' plot(fit)
#'
#' # Constrained AR(1) model
#' fit <- ARFPOP(dat = sim$fl, gam = 0.95, lambda = 1, constraint = T, compute_fitted_values = T)
#' print(fit)
#' plot(fit)
#' 
#' @seealso
#' \strong{Estimate spikes:}
#' \code{\link{ARFPOP}}
#' \code{\link{fit_ar1_model}}
#'
#' @export

ARFPOP <- structure(function(dat, gam, lambda, constraint = FALSE, compute_fitted_values = FALSE, EPS = 1e-10) {
  stopifnot(gam > 0 && gam <= 1)
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
    compute_fitted_values = compute_fitted_values, 
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
  
  if (compute_fitted_values) {
    fitted_values = result.list$mean.vec  
  } else {
    fitted_values = NA 
  }
  
  
  if (result.list$success) {
    out <-
      list(
        spikes = spikes,
        fittedValues = fitted_values,
        dat = dat,
        type = paste0("ar1-fpop", constraint_str), 
        changePts = changePts,
        call = match.call(),
        gam = gam,
        lambda = lambda,
        cost = as.numeric(result.list$cost.mat),
        nIntervals = as.numeric(result.list$intervals.mat),
        end_vec = result.list$ends.vec,
        EPS = EPS
      )
    class(out) <- "estimatedSpikes"
    return(out)
  } else {
    stop("Numerical issues!")
  }

  
})