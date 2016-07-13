PeakSegFPOP <- structure(function
### Find the optimal change-points using the Poisson loss and the
### PeakSeg constraint. For N data points and S segments, the
### functional pruning algorithm is O(S*N) space and O(S*NlogN)
### time. It recovers the exact solution to the following optimization
### problem. Let Z be an N-vector of count data (non-negative
### integers). Find the N-vector M of real numbers (segment means)
### which minimize the Poisson Loss, sum_i m_i - z_i * log(m_i),
### subject to constraints: (1) there are at most S changes in M, and
### (2) up changes are followed by down changes, and vice versa (mu1
### <= mu2 >= mu3 <= mu4 >= mu5, etc). Note that the segment means can
### be equal, in which case the recovered model is not feasible for
### the PeakSeg problem.
(count.vec,
### integer vector of count data.
 weight.vec=rep(1, length(count.vec)),
### numeric vector (same length as count.vec) of positive weights.
 penalty=NULL
### numeric of length 1: penalty parameter (smaller for more peaks,
### larger for fewer peaks).
){
  n.data <- length(count.vec)
  stopifnot(is.integer(count.vec))
  stopifnot(is.numeric(weight.vec))
  stopifnot(n.data==length(weight.vec))
  stopifnot(0 < weight.vec)
  stopifnot(is.numeric(penalty))
  stopifnot(length(penalty)==1)
  stopifnot(0 <= penalty)
  cost.mat <- double(n.data*2)
  ends.vec <- integer(n.data)
  mean.vec <- double(n.data)
  intervals.mat <- integer(n.data*2)
  result.list <- .C(
    "PeakSegFPOPLog_interface",
    count.vec=as.integer(count.vec),
    weight.vec=as.numeric(weight.vec),
    n.data=as.integer(n.data),
    penalty=as.numeric(penalty),
    cost.mat=as.double(cost.mat),
    ends.vec=as.integer(ends.vec),
    mean.vec=as.double(mean.vec),
    intervals.mat=as.integer(intervals.mat),
    PACKAGE="coseg")
  result.list$cost.mat <- matrix(
    result.list$cost.mat, 2, n.data, byrow=TRUE)
  result.list$intervals.mat <- matrix(
    result.list$intervals.mat, 2, n.data, byrow=TRUE)
  result.list
### List of model parameters. count.vec, weight.vec, n.data,
### max.segments (input parameters), cost.mat (optimal Poisson loss),
### ends.mat (optimal position of segment ends, 1-indexed), mean.mat
### (optimal segment means), intervals.mat (number of intervals stored
### by the functional pruning algorithm).
}, ex=function(){
  library(coseg)
  data("H3K4me3_XJ_immune_chunk1")
  by.sample <-
    split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
  n.data.vec <- sapply(by.sample, nrow)
  one <- by.sample[[1]]
  count.vec <- one$coverage
  weight.vec <- with(one, chromEnd-chromStart)
  count.vec <- as.integer(c(1, 10, 14, 13))
  weight.vec <- c(1,1,1,1)
  fit <- PeakSegFPOP(count.vec, weight.vec, 1)
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
      "intervals stored by the\nconstrained optimal segmentation algorithm",
      breaks=c(20, 40))
})

