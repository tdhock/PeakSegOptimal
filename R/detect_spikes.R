#' Estimate spike train, underlying calcium concentration, and changepoints based on 
#' fluorescence trace for a range of tuning parameters. 
#' 
#' Given [lambda_min, lambda_max] this procedure selects lambdas such that each
#' lambda is associated with a different number of spike events. 
#'
#' @param dat fluorescence data
#' @param gam a scalar value for the AR(1) decay parameter.
#' @param constraint boolean specifying whether constrained or unconstrained optimization
#' problem
#' @param EPS double specfying the minimum calcium value
#' @param min_lambda smallest tuning parameter
#' @param max_lambda largest tuning parameter
#'
#' @return Returns a list with elements:
#' @return \code{spikes} the set of spikes
#' @return \code{fittedValues} estimated calcium concentration
#' @return \code{changePts} the set of changepoints
#' @return \code{cost} the cost at each time point (vector)
#' @return \code{nIntervals} the number of intervals at each point (vector)
#'
#' @details
#'
#' This algorithm solves the optimization problems
#' 
#'  \strong{AR(1) model:}
#'  
#'  minimize_{c1,...,cT} 0.5 sum_{t=1}^T ( y_t - c_t )^2 + lambda sum_{t=2}^T 1[c_t != max(gam c_{t-1}, EPS)]
#'  
#'  for the global optimum, where y_t is the observed fluorescence at the tth timepoint.
#'
#' \strong{Constrained AR(1) model}
#' 
#' minimize_{c1,...,cT} 0.5 sum_{t=1}^T ( y_t - c_t )^2 + lambda sum_{t=2}^T 1[c_t != max(gam c_{t-1}, EPS)]
#' 
#' c_t >= max(gam c_{t-1}, EPS), t = 2, ..., T
#'
#' See 
#' 
#' 1. Jewell, Hocking, Fearnhead, and Witten (2018) <arXiv:1802.07380>  
#' 
#' 2. Jewell and Witten (2017) <arXiv:1703.08644> 
#'
#'
#' @examples
#'
#'
#' @seealso
#'
#'
#' @export
detect_spikes <- function(dat, gam, constraint, EPS, lambda_min, lambda_max, sensitivity = 1) { 
  out <- crops(y, gam, lambda_min, lambda_max, constraint, EPS, sensitivity)
  return(out)
}
