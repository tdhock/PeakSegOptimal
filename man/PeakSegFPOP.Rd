\name{PeakSegFPOP}
\alias{PeakSegFPOP}
\title{PeakSegFPOP}
\description{Find the optimal change-points using the Poisson loss and the
PeakSeg constraint. For N data points, the functional pruning
algorithm is O(N) space and O(NlogN) time. It recovers the exact
solution to the following optimization problem. Let Z be an
N-vector of count data (non-negative integers). Find the N-vector
M of real numbers (segment means) which minimize the penalized
Poisson Loss, sum_{i=2}^N I(m_i != m_{i-1})*penalty + sum_{i=1}^N
m_i - z_i * log(m_i), subject to constraint: up changes are
followed by down changes, and vice versa. Note that the segment
means can be equal, in which case the recovered model is not
feasible for the PeakSeg problem. Unlike PeakSegPDPA which forces
the first segment mean to be down (mu1 <= mu2), PeakSegFPOP may
recover a model with the first segment mean up (mu1 >= mu2).}
\usage{PeakSegFPOP(count.vec, weight.vec = rep(1, length(count.vec)), 
    penalty = NULL)}
\arguments{
  \item{count.vec}{integer vector of count data.}
  \item{weight.vec}{numeric vector (same length as count.vec) of positive weights.}
  \item{penalty}{numeric of length 1: penalty parameter (smaller for more peaks,
larger for fewer peaks).}
}

\value{List of model parameters. count.vec, weight.vec, n.data, penalty
(input parameters), cost.mat (optimal Poisson loss), ends.vec
(optimal position of segment ends, 1-indexed), mean.vec (optimal
segment means), intervals.mat (number of intervals stored by the
functional pruning algorithm), label.vec (1=up or 0=down).}

\author{Toby Dylan Hocking}




\examples{
library(coseg)
data("H3K4me3_XJ_immune_chunk1")
by.sample <-
  split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
n.data.vec <- sapply(by.sample, nrow)
one <- by.sample[[1]]
count.vec <- one$coverage
weight.vec <- with(one, chromEnd-chromStart)
fit <- PeakSegFPOP(count.vec, weight.vec, 1000)
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
    "intervals stored by the\nconstrained optimal segmentation algorithm")
}
