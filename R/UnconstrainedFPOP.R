UnconstrainedFPOP <- structure(function
### Find the optimal change-points using the Poisson loss and no
### constraints. For N data points, the functional pruning algorithm
### is O(N log N) time and memory. It recovers the exact solution to
### the following optimization problem. Let Z be an N-vector of count
### data (count.vec, non-negative integers), let W be an N-vector of
### positive weights (weight.vec), and let penalty be a non-negative
### real number. Find the N-vector M of real numbers (segment means)
### and (N-1)-vector C of change-point indicators in 0,1 which
### minimize the penalized Poisson Loss, penalty*sum_[i=1]^[N_1]
### I(c_i=1) + sum_[i=1]^N w_i*[m_i-z_i*log(m_i)].
(count.vec,
### integer vector of length >= 3: non-negative count data to segment.
 weight.vec=rep(1, length(count.vec)),
### numeric vector (same length as count.vec) of positive weights.
 penalty=NULL,
### non-negative numeric scalar: penalty parameter (smaller for more
### peaks, larger for fewer peaks).
  verbose_file=""
### String: file name to output candidate indices (slow but useful for
### analysis of the algorithm). Default "" means to not output indices
### (fast).
){
  n.data <- length(count.vec)
  stopifnot(1 <= n.data)
  stopifnot(is.integer(count.vec))
  stopifnot(0 <= count.vec)
  stopifnot(is.numeric(weight.vec))
  stopifnot(n.data==length(weight.vec))
  stopifnot(0 < weight.vec)
  stopifnot(is.numeric(penalty))
  stopifnot(length(penalty)==1)
  stopifnot(0 <= penalty)
  cost.vec <- double(n.data)
  ends.vec <- integer(n.data)
  mean.vec <- double(n.data)
  intervals.vec <- integer(n.data)
  result.list <- .C(
    "UnconstrainedFPOP_interface",
    count.vec=as.integer(count.vec),
    weight.vec=as.numeric(weight.vec),
    n.data=as.integer(n.data),
    penalty=as.numeric(penalty),
    verbose_file=as.character(verbose_file),
    cost.vec=as.double(cost.vec),
    ends.vec=as.integer(ends.vec),
    mean.vec=as.double(mean.vec),
    intervals.vec=as.integer(intervals.vec),
    PACKAGE="PeakSegOptimal")
  ## 1-indexed segment ends!
  result.list$ends.vec <- result.list$ends.vec+1L
  result.list
### List of model parameters. count.vec, weight.vec, n.data, penalty
### (input parameters), cost.vec (optimal Poisson loss), ends.vec
### (optimal position of segment ends, 1-indexed), mean.vec (optimal
### segment means), intervals.vec (number of intervals stored by the
### functional pruning algorithm). To recover the solution in terms of
### (M,C) variables, see the example.
}, ex=function(){

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
  stopifnot(cumsum(sign(C.vec)) %in% c(0, 1))

  ## Compute penalized Poisson loss of M.vec and compare to the value reported
  ## in the fit solution list.
  n.peaks <- sum(C.vec==1)
  rbind(
    n.peaks*penalty + PoissonLoss(count.vec, M.vec, weight.vec),
    fit$cost.vec[2, fit$n.data])

  ## Plot the number of intervals stored by the algorithm.
  FPOP.intervals <- data.frame(
    label=ifelse(as.numeric(row(fit$intervals.vec))==1, "up", "down"),
    data=as.numeric(col(fit$intervals.vec)),
    intervals=as.numeric(fit$intervals.vec))
  library(ggplot2)
  ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(label ~ .)+
    geom_line(aes(data, intervals), data=FPOP.intervals)+
    scale_y_continuous(
      "intervals stored by the\nconstrained optimal segmentation algorithm")

})

