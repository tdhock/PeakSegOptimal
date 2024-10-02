\name{PeakSegPDPA}
\alias{PeakSegPDPA}
\title{PeakSegPDPA}
\description{Find the optimal change-points using the Poisson loss and the
PeakSeg constraint. For N data points and S segments, the
functional pruning algorithm is O(S*NlogN) space and O(S*NlogN)
time. It recovers the exact solution to the following optimization
problem. Let Z be an N-vector of count data (\code{count.vec},
non-negative integers) and let W be an N-vector of positive
weights (\code{weight.vec}). Find the N-vector M of real numbers (segment
means) and (N-1)-vector C of change-point indicators in -1,0,1
which minimize the Poisson Loss, sum_[i=1]^N
w_i*[m_i-z_i*log(m_i)], subject to constraints: (1) there are
exactly S-1 non-zero elements of C, and (2) the first change is up
and the next change is down, etc (sum_[i=1]^t c_i in 0,1 for all
t<N), and (3) Every zero-valued change-point variable has an equal
segment mean after: c_i=0 implies m_i=m_[i+1], (4) every
positive-valued change-point variable may have an up change after:
c_i=1 implies m_i<=m_[i+1], (5) every negative-valued change-point
variable may have a down change after: c_i=-1 implies
m_i>=m_[i+1]. Note that when the equality constraints are active
for non-zero change-point variables, the recovered model is not
feasible for the strict inequality constraints of the PeakSeg
problem, and the optimum of the PeakSeg problem is undefined.}
\usage{PeakSegPDPA(count.vec, 
    weight.vec = rep(1, 
        length(count.vec)), 
    max.segments = NULL)}
\arguments{
  \item{count.vec}{integer vector of count data.}
  \item{weight.vec}{numeric vector (same length as \code{count.vec}) of positive weights.}
  \item{max.segments}{integer of length 1: maximum number of segments (must be >= 2).}
}

\value{List of model parameters. \code{count.vec}, \code{weight.vec}, n.data,
\code{max.segments} (input parameters), cost.mat (optimal Poisson loss),
ends.mat (optimal position of segment ends, 1-indexed), mean.mat
(optimal segment means), intervals.mat (number of intervals stored
by the functional pruning algorithm). To recover the solution in
terms of (M,C) variables, see the example.}

\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre]}




\examples{

## Use the algo to compute the solution list.
data("H3K4me3_XJ_immune_chunk1", envir=environment())
by.sample <-
  split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
n.data.vec <- sapply(by.sample, nrow)
one <- by.sample[[1]]
count.vec <- one$coverage
weight.vec <- with(one, chromEnd-chromStart)
max.segments <- 19L
fit <- PeakSegPDPA(count.vec, weight.vec, max.segments)

## Recover the solution in terms of (M,C) variables.
n.segs <- 11L
change.vec <- fit$ends.mat[n.segs, 2:n.segs]
change.sign.vec <- rep(c(1, -1), length(change.vec)/2)
end.vec <- c(change.vec, fit$n.data)
start.vec <- c(1, change.vec+1)
length.vec <- end.vec-start.vec+1
mean.vec <- fit$mean.mat[n.segs, 1:n.segs]
M.vec <- rep(mean.vec, length.vec)
C.vec <- rep(0, fit$n.data-1)
C.vec[change.vec] <- change.sign.vec
diff.vec <- diff(M.vec)
data.frame(
  change=c(C.vec, NA),
  mean=M.vec,
  equality.constraint.active=c(sign(diff.vec) != C.vec, NA))
stopifnot(cumsum(sign(C.vec)) \%in\% c(0, 1))

## Compute Poisson loss of M.vec and compare to the value reported
## in the fit solution list.
rbind(
  PoissonLoss(count.vec, M.vec, weight.vec),
  fit$cost.mat[n.segs, fit$n.data])

## Plot the number of intervals stored by the algorithm.
PDPA.intervals <- data.frame(
  segments=as.numeric(row(fit$intervals.mat)),
  data=as.numeric(col(fit$intervals.mat)),
  intervals=as.numeric(fit$intervals.mat))
some.intervals <- subset(PDPA.intervals, segments<data & 1<segments)
library(ggplot2)
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(segments ~ .)+
  geom_line(aes(data, intervals), data=some.intervals)+
  scale_y_continuous(
    "intervals stored by the\\nconstrained optimal segmentation algorithm",
    breaks=c(20, 40))

}
