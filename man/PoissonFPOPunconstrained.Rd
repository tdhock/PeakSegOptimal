\name{PoissonFPOPunconstrained}
\alias{PoissonFPOPunconstrained}
\title{Unconstrained Poisson FPOP}
\description{Optimal changepoint detection using the Poisson loss
without any constraint on segment mean changes (no up-down PeakSeg
constraint). Uses the functional pruning optimal partitioning (FPOP)
algorithm which is O(N log N) time and memory.}
\usage{PoissonFPOPunconstrained(count.vec,
    weight.vec = rep(1,
        length(count.vec)),
    penalty = NULL)}
\arguments{
  \item{count.vec}{integer vector of non-negative count data to segment.}
  \item{weight.vec}{numeric vector (same length as count.vec) of positive weights.}
  \item{penalty}{non-negative numeric scalar: penalty for each changepoint.}
}
\value{List with elements \code{seg.mean} (segment means in
chronological order), \code{seg.end} (1-indexed segment ends),
\code{n.segments} (number of segments), \code{cost.vec},
\code{end.vec}, \code{mean.vec}, \code{intervals.vec}.}
\examples{
library(PeakSegOptimal)
data.vec <- as.integer(c(3,5,7,4,6,2, 18,22,15,20,19,25, 4,6,3,5,7,2))
fit <- PoissonFPOPunconstrained(data.vec, penalty=10)
fit$n.segments
fit$seg.mean
}
