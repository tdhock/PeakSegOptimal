\name{PeakSegFPOP}
\alias{PeakSegFPOP}
\title{PeakSegFPOP}
\description{Find the optimal change-points using the Poisson loss and the
PeakSeg constraint. For N data points, the functional pruning
algorithm is O(N log N) time and memory. It recovers the exact
solution to the following optimization problem. Let Z be an
N-vector of count data (\code{count.vec}, non-negative integers), let W
be an N-vector of positive weights (\code{weight.vec}), and let \code{penalty}
be a non-negative real number. Find the N-vector M of real numbers
(segment means) and (N-1)-vector C of change-point indicators in
{-1,0,1} which minimize the penalized Poisson Loss,
penalty*sum_{i=1}^{N_1} I(c_i=1) + sum_{i=1}^N
w_i*[m_i-z_i*log(m_i)], subject to constraints: (1) the first
change is up and the next change is down, etc (sum_{i=1}^t c_i in
{0,1} for all t<N-1), and (2) the last change is down
0=sum_{i=1}^{N-1}c_i, and (3) Every zero-valued change-point
variable has an equal segment mean after: c_i=0 implies
m_i=m_{i+1}, (4) every positive-valued change-point variable may
have an up change after: c_i=1 implies m_i<=m_{i+1}, (5) every
negative-valued change-point variable may have a down change
after: c_i=-1 implies m_i>=m_{i+1}. Note that when the equality
constraints are active for non-zero change-point variables, the
recovered model is not feasible for the strict inequality
constraints of the PeakSeg problem, and the optimum of the PeakSeg
problem is undefined.}
\usage{PeakSegFPOP(count.vec, weight.vec = rep(1, length(count.vec)), 
    penalty = NULL)}
\arguments{
  \item{count.vec}{integer vector of length >= 3: non-negative count data to segment.}
  \item{weight.vec}{numeric vector (same length as \code{count.vec}) of positive weights.}
  \item{penalty}{non-negative numeric scalar: \code{penalty} parameter (smaller for more
peaks, larger for fewer peaks).}
}

\value{List of model parameters. \code{count.vec}, \code{weight.vec}, n.data, \code{penalty}
(input parameters), cost.mat (optimal Poisson loss), ends.vec
(optimal position of segment ends, 1-indexed), mean.vec (optimal
segment means), intervals.mat (number of intervals stored by the
functional pruning algorithm). To recover the solution in terms of
(M,C) variables, see the example.}

\author{Toby Dylan Hocking}




\examples{

## Use the algo to compute the solution list.
library(PeakSegOptimal)
data("H3K4me3_XJ_immune_chunk1", envir=environment())
by.sample <-
  split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
n.data.vec <- sapply(by.sample, nrow)
one <- by.sample[[1]]
count.vec <- one$coverage
weight.vec <- with(one, chromEnd-chromStart)
penalty <- 1000
fit <- PeakSegFPOP(count.vec, weight.vec, penalty)

## Recover the solution in terms of (M,C) variables.
change.vec <- with(fit, rev(ends.vec[ends.vec>0]))
change.sign.vec <- rep(c(1, -1), length(change.vec)/2)
end.vec <- c(change.vec, fit$n.data)
start.vec <- c(1, change.vec+1)
length.vec <- end.vec-start.vec+1
mean.vec <- rev(fit$mean.vec[1:(length(change.vec)+1)])
M.vec <- rep(mean.vec, length.vec)
C.vec <- rep(0, fit$n.data-1)
C.vec[change.vec] <- change.sign.vec
diff.vec <- diff(M.vec)
data.frame(
  change=c(C.vec, NA),
  mean=M.vec,
  equality.constraint.active=c(sign(diff.vec) != C.vec, NA))
stopifnot(cumsum(sign(C.vec)) \%in\% c(0, 1))

## Compute penalized Poisson loss of M.vec and compare to the value reported
## in the fit solution list.
n.peaks <- sum(C.vec==1)
rbind(
  n.peaks*penalty + PoissonLoss(count.vec, M.vec, weight.vec),
  fit$cost.mat[2, fit$n.data])

## Plot the number of intervals stored by the algorithm.
FPOP.intervals <- data.frame(
  label=ifelse(as.numeric(row(fit$intervals.mat))==1, "up", "down"),
  data=as.numeric(col(fit$intervals.mat)),
  intervals=as.numeric(fit$intervals.mat))
library(ggplot2)
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(label ~ .)+
  geom_line(aes(data, intervals), data=FPOP.intervals)+
  scale_y_continuous(
    "intervals stored by the\\nconstrained optimal segmentation algorithm")

}
