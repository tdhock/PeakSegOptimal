#' Unconstrained optimal partitioning for Poisson loss
#'
#' This function implements the FPOP algorithm for finding optimal
#' segmentation of count data using Poisson loss with no constraints.
#'
#' @param count.vec integer vector of non-negative counts to segment
#' @param weight.vec numeric vector of positive weights (default: all 1)
#' @param penalty non-negative numeric: penalty parameter for segment creation
#'   (smaller = more segments, larger = fewer segments)
#'
#' @return A list containing:
#'   \item{count.vec}{input count vector}
#'   \item{weight.vec}{input weight vector}
#'   \item{n.data}{number of data points}
#'   \item{penalty}{penalty parameter used}
#'   \item{cost.mat}{optimal Poisson loss}
#'   \item{ends.vec}{1-indexed segment end positions}
#'   \item{mean.vec}{optimal segment means}
#'   \item{intervals.mat}{number of pieces in cost function}
#'
#' @examples
#' # example: detect peaks and valleys in stock price data
#' data <- c(10, 12, 15, 18, 14, 8, 5, 9, 7, 3)
#' result <- PeakSegFPOPLogUnconstrained(
#'   as.integer(data),
#'   penalty = 2.5
#' )
#'
#' @export
PeakSegFPOPLogUnconstrained <- function(
  count.vec,
  weight.vec = rep(1, length(count.vec)),
  penalty = NULL){
  
  n.data <- length(count.vec)
  
  stopifnot(3 <= n.data)
  stopifnot(is.integer(count.vec))
  stopifnot(0 <= count.vec)
  stopifnot(is.numeric(weight.vec))
  stopifnot(n.data == length(weight.vec))
  stopifnot(0 < weight.vec)
  stopifnot(is.numeric(penalty))
  stopifnot(length(penalty) == 1)
  stopifnot(0 <= penalty)
  
  # allocating output arrays (data_count × 1, not × 2)..
  cost.mat <- double(n.data)
  ends.vec <- integer(n.data)
  mean.vec <- double(n.data)
  intervals.mat <- integer(n.data)
  
  result.list <- .C(
    "PeakSegFPOPLogUnconstrained_interface",
    count.vec = as.integer(count.vec),
    weight.vec = as.numeric(weight.vec),
    n.data = as.integer(n.data),
    penalty = as.numeric(penalty),
    cost.mat = as.double(cost.mat),
    ends.vec = as.integer(ends.vec),
    mean.vec = as.double(mean.vec),
    intervals.mat = as.integer(intervals.mat),
    PACKAGE = "PeakSegOptimal"
  )
  
  result.list$ends.vec <- result.list$ends.vec + 1L  # convert to 1-indexed
  result.list$cost.mat <- result.list$cost.mat * cumsum(result.list$weight.vec)
  
  result.list
}