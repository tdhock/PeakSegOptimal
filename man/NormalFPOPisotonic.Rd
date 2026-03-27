\name{NormalFPOPisotonic}
\alias{NormalFPOPisotonic}
\title{Isotonic Regression via Normal FPOP}
\description{Regularized isotonic regression using the Normal
(Gaussian/quadratic) loss with a non-decreasing constraint on
segment means. Uses the functional pruning optimal partitioning (FPOP)
algorithm.}
\usage{NormalFPOPisotonic(data.vec,
    penalty = NULL)}
\arguments{
  \item{data.vec}{numeric vector of data to segment.}
  \item{penalty}{non-negative numeric scalar: penalty for each changepoint.
  Use \code{penalty=0} for unpenalized isotonic regression (equivalent to
  \code{\link{isoreg}}).}
}
\value{List with elements \code{seg.mean} (segment means in
chronological order, non-decreasing), \code{seg.end} (1-indexed segment
ends), \code{n.segments} (number of segments), \code{cost.vec},
\code{end.vec}, \code{mean.vec}, \code{intervals.vec}.}
\examples{
library(PeakSegOptimal)
y <- c(5, 3, 4, 1, 2, 8, 7, 9)
fit <- NormalFPOPisotonic(y, penalty=0)
fit$seg.mean
}
