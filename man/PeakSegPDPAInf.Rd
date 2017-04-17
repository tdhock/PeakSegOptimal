\name{PeakSegPDPAInf}
\alias{PeakSegPDPAInf}
\title{PeakSegPDPAInf}
\description{Find the optimal change-points using the Poisson loss and the
PeakSeg constraint. This function is an interface to the C++ code
which always uses -Inf for the first interval's lower limit and
Inf for the last interval's upper limit -- it is for testing the
number of intervals between the two implementations.}
\usage{PeakSegPDPAInf(count.vec, weight.vec = rep(1, length(count.vec)), 
    max.segments = NULL)}
\arguments{
  \item{count.vec}{integer vector of count data.}
  \item{weight.vec}{numeric vector (same length as count.vec) of positive weights.}
  \item{max.segments}{integer of length 1: maximum number of segments (must be >= 2).}
}

\value{List of model parameters. count.vec, weight.vec, n.data,
max.segments (input parameters), cost.mat (optimal Poisson loss),
ends.mat (optimal position of segment ends, 1-indexed), mean.mat
(optimal segment means), intervals.mat (number of intervals stored
by the functional pruning algorithm). To recover the solution in
terms of (M,C) variables, see the example.}

\author{Toby Dylan Hocking}




\examples{

## Use the algo to compute the solution list.
library(coseg)
data("H3K4me3_XJ_immune_chunk1", envir=environment())
by.sample <-
  split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
n.data.vec <- sapply(by.sample, nrow)
one <- by.sample[[1]]
count.vec <- one$coverage
weight.vec <- with(one, chromEnd-chromStart)
max.segments <- 19L

library(data.table)
ic.list <- list()
for(fun.name in c("PeakSegPDPA", "PeakSegPDPAInf")){
  fun <- get(fun.name)
  fit <- fun(count.vec, weight.vec, max.segments)
  ic.list[[fun.name]] <- data.table(
    fun.name,
    segments=as.numeric(row(fit$intervals.mat)),
    data=as.numeric(col(fit$intervals.mat)),
    cost=as.numeric(fit$cost.mat),
    intervals=as.numeric(fit$intervals.mat))
}
ic <- do.call(rbind, ic.list)[0 < intervals]
intervals <- dcast(ic, data + segments ~ fun.name, value.var="intervals")
cost <- dcast(ic, data + segments ~ fun.name, value.var="cost")
not.equal <- cost[PeakSegPDPA != PeakSegPDPAInf]
stopifnot(nrow(not.equal)==0)

intervals[, increase := PeakSegPDPAInf-PeakSegPDPA]
table(intervals$increase)
quantile(intervals$increase)
ic[, list(
  mean=mean(intervals),
  max=max(intervals)
  ), by=list(fun.name)]

}
