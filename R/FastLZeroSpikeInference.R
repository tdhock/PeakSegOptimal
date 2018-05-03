#' FastLZeroSpikeInference: FastLZeroSpikeInference: A package for estimating spike
#' times from calcium imaging data using an L0 penalty
#'
#' This package implements an algorithm for deconvolving calcium imaging data
#' for a single neuron in order to estimate the times at which the neuron
#' spikes.
#'
#'
#' @seealso
#' \strong{Estimate spikes:}
#' \code{\link{ARFPOP.R}}

#'
#'
#' @details
#'
#' This package implements an algorithm for deconvolving calcium imaging data for a single 
#' neuron in order to estimate the times at which the neuron spikes (Jewell et al. 2018). 
#' This algorithim is an extension of the constrained functional pruning algorithm of 
#' Hocking et al. (2017). 
#'
#' This algorithm solves the optimization problems
#' ### AR(1) model
#' minimize_{c1,...,cT} 0.5 sum_{t=1}^T ( y_t - c_t )^2 + lambda sum_{t=2}^T 1_{c_t != max(gam c_{t-1}, EPS) }
#'
#' for the global optimum, where y_t is the observed fluorescence at the tth timepoint.
#'
#' ### Constrained AR(1) model 
#' minimize_{c1,...,cT} 0.5 sum_{t=1}^T ( y_t - c_t )^2 + lambda sum_{t=2}^T 1_{c_t != max(gam c_{t-1}, EPS) }
#'
#' subject to c_t >= max(gam c_{t-1}, EPS), t = 2, ..., T
#'
#' See Jewell, Hocking, Fearnhead, and Witten (2018) <arXiv:1802.07380> and 
#' Jewell and Witten (2017) <arXiv:1703.08644> and 
#'  
#'
#' @examples
#'
#' @docType package
#' @name FastLZeroSpikeInference
NULL
