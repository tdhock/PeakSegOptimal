#' FastLZeroSpikeInference: FastLZeroSpikeInference: A package for estimating spike
#' times from calcium imaging data using an L0 penalty
#'
#' This package implements an algorithm for deconvolving calcium imaging data
#' for a single neuron in order to estimate the times at which the neuron
#' spikes.
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
#' We introduce the constant EPS > 0, typically on the order of 10^-10, to avoid 
#' arbitrarily small calcium concentrations that would result in numerical  
#' instabilities. In practice, this means that the estimated calcium concentration 
#' decays according to the AR(1) for values greater than EPS and is equal to EPS thereafter.
#'
#' When estimating the spikes, it is not necessary to explicitly compute the 
#' calcium concentration. Therefore, if only the spike times are required, the user
#' can avoid this computation cost by setting the compute_fitted_values boolean to false. 
#' By default, the calcium concentration is not estimated. 
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
#' library(LZeroSpikeInference)
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
#' fit <- estimate_spikes(dat = sim$fl, gam = 0.95, lambda = 1, compute_fitted_values = T)
#' plot(fit)
#'
#' # Constrained AR(1) model
#' fit <- estimate_spikes(dat = sim$fl, gam = 0.95, lambda = 1, constraint = T, compute_fitted_values = T)
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
#' @export
#'
#' @docType package
#' @name FastLZeroSpikeInference
NULL
