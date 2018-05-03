#' Estimate spike train, underlying calcium concentration, and changepoints based on fluorescence
#' trace.
#' 
#' @param fit object created by running ARFPOP
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
#'  minimize_{c1,...,cT} 0.5 sum_{t=1}^T ( y_t - c_t )^2 + lambda sum_{t=2}^T 1_{c_t != max(gam c_{t-1}, EPS) }
#'  
#'  for the global optimum, where y_t is the observed fluorescence at the tth timepoint.
#'
#' \strong{Constrained AR(1) model}
#' minimize_{c1,...,cT} 0.5 sum_{t=1}^T ( y_t - c_t )^2 + lambda sum_{t=2}^T 1_{c_t != max(gam c_{t-1}, EPS) }
#' 
#' c_t >= max(gam c_{t-1}, EPS), t = 2, ..., T
#'
#' See Jewell, Hocking, Fearnhead, and Witten (2018) <arXiv:1802.07380> and 
#' Jewell and Witten (2017) <arXiv:1703.08644> 
#'   
#' @export
#'
fit_ar1_model <- structure(function(fit) {
  
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
      fittedValues = result.list$mean.vec,
      dat = fit$dat,
      type = fit$type, 
      changePts = fit$changePts,
      call = fit$call,
      gam = fit$gam,
      lambda = fit$lambda,
      cost = fit$cost,
      nIntervals = fit$nIntervals,
      end_vec = fit$end_vec,
      EPS = fit$EPS
    )
  class(out) <- "estimatedSpikes"
  return(out)
})