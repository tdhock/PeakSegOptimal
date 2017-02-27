PeakSegFPOP <- structure(function
### Find the optimal change-points using the Poisson loss and the
### PeakSeg constraint. For N data points, the functional pruning
### algorithm is O(N log N) time and memory. It recovers the exact
### solution to the following optimization problem. Let Z be an
### N-vector of count data (count.vec, non-negative integers), let W
### be an N-vector of positive weights (weight.vec), and let penalty
### be a non-negative real number. Find the N-vector M of real numbers
### (segment means) and (N-1)-vector C of change-point indicators in
### {-1,0,1} which minimize the penalized Poisson Loss,
### penalty*sum_{i=1}^{N_1} I(c_i=1) + sum_{i=1}^N
### w_i*[m_i-z_i*log(m_i)], subject to constraints: (1) the first
### change is up and the next change is down, etc (sum_{i=1}^t c_i in
### {0,1} for all t<N-1), and (2) the last change is down
### 0=sum_{i=1}^{N-1}c_i, and (3) Every zero-valued change-point
### variable has an equal segment mean after: c_i=0 implies
### m_i=m_{i+1}, (4) every positive-valued change-point variable may
### have an up change after: c_i=1 implies m_i<=m_{i+1}, (5) every
### negative-valued change-point variable may have a down change
### after: c_i=-1 implies m_i>=m_{i+1}. Note that when the equality
### constraints are active for non-zero change-point variables, the
### recovered model is not feasible for the strict inequality
### constraints of the PeakSeg problem, and the optimum of the PeakSeg
### problem is undefined.
(count.vec,
### integer vector of length >= 3: non-negative count data to segment.
 weight.vec=rep(1, length(count.vec)),
### numeric vector (same length as count.vec) of positive weights.
 penalty=NULL
### non-negative numeric scalar: penalty parameter (smaller for more
### peaks, larger for fewer peaks).
){
  n.data <- length(count.vec)
  stopifnot(3 <= n.data)
  stopifnot(is.integer(count.vec))
  stopifnot(0 <= count.vec)
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
    ##label.vec=as.integer(label.vec),
    PACKAGE="coseg")
  ## 1-indexed segment ends!
  result.list$ends.vec <- result.list$ends.vec+1L
  result.list$cost.mat <- matrix(
    result.list$cost.mat*cumsum(weight.vec), 2, n.data, byrow=TRUE)
  result.list$intervals.mat <- matrix(
    result.list$intervals.mat, 2, n.data, byrow=TRUE)
  result.list
### List of model parameters. count.vec, weight.vec, n.data, penalty
### (input parameters), cost.mat (optimal Poisson loss), ends.vec
### (optimal position of segment ends, 1-indexed), mean.vec (optimal
### segment means), intervals.mat (number of intervals stored by the
### functional pruning algorithm). To recover the solution in terms of
### (M,C) variables, see the example.
}, ex=function(){

  ## Use the algo to compute the solution list.
  library(coseg)
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
      "intervals stored by the\nconstrained optimal segmentation algorithm")

})

PeakSegFPOPchrom <- structure(function
### Find the optimal change-points using the Poisson loss and the
### PeakSeg constraint. This function is a user-friendly interface to
### the PeakSegFPOP function.
(count.df,
### data.frame with columns count, chromStart, chromEnd.
 penalty=NULL
### non-negative numeric scalar: penalty parameter (smaller for more
### peaks, larger for fewer peaks).
){
  stopifnot(is.data.frame(count.df))
  n.data <- nrow(count.df)
  stopifnot(3 <= n.data)
  stopifnot(is.integer(count.df$chromStart))
  stopifnot(is.integer(count.df$chromEnd))
  stopifnot(is.integer(count.df$count))
  stopifnot(count.df$chromStart < count.df$chromEnd)
  stopifnot(0 <= count.df$chromStart)
  stopifnot(0 <= count.df$count)
  weight.vec <- with(count.df, chromEnd - chromStart)
  stopifnot(is.numeric(penalty))
  stopifnot(length(penalty)==1)
  stopifnot(0 <= penalty)
  fit <- PeakSegFPOP(count.df$count, weight.vec, penalty)
  break.vec <- rev(fit$ends.vec[0<fit$ends.vec])
  first <- c(1, break.vec+1)
  last <- c(break.vec, nrow(count.df))
  ##label.vec <- rev(fit$label.vec[0 <= fit$label.vec])
  label.vec <- rep(c(0, 1), l=sum(is.finite(fit$mean.vec)))
  status.str <- ifelse(label.vec==0, "background", "peak")
  peaks <- sum(label.vec==1)
  mean.vec <- rev(fit$mean.vec[is.finite(fit$mean.vec)])
  list(
    segments=data.frame(
      mean=mean.vec,
      first,
      last,
      chromStart=count.df$chromStart[first],
      chromEnd=count.df$chromEnd[last],
      status=factor(status.str, c("background", "peak")),
      peaks,
      segments=length(first)),
    loss=data.frame(
      segments=length(first),
      peaks,
      penalized.loss=min(fit$cost.mat[, n.data]),
      feasible=all(diff(mean.vec)!=0)
      )
    )
### List of data.frames: segments can be used for plotting the
### segmentation model, loss summarizes the penalized PoissonLoss and
### feasibilty of the computed model.
}, ex=function(){

  library(coseg)
  data("H3K4me3_XJ_immune_chunk1", envir=environment())
  sample.id <- "McGill0106"
  H3K4me3_XJ_immune_chunk1$count <- H3K4me3_XJ_immune_chunk1$coverage
  by.sample <-
    split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
  one.sample <- by.sample[[sample.id]]

  penalty.constant <- 1100
  fpop.fit <- PeakSegFPOPchrom(one.sample, penalty.constant)
  fpop.breaks <- subset(fpop.fit$segments, 1 < first)
  library(ggplot2)
  ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    geom_step(aes(chromStart/1e3, coverage),
              data=one.sample, color="grey")+
    geom_segment(aes(chromStart/1e3, mean,
                     xend=chromEnd/1e3, yend=mean),
                 color="green",
                 data=fpop.fit$segments)+
    geom_vline(aes(xintercept=chromStart/1e3),
               color="green",
               linetype="dashed",
               data=fpop.breaks)

  max.peaks <- as.integer(fpop.fit$segments$peaks[1]+1)
  pdpa.fit <- PeakSegPDPAchrom(one.sample, max.peaks)
  models <- pdpa.fit$modelSelection.decreasing
  models$PoissonLoss <- pdpa.fit$loss[paste(models$peaks), "PoissonLoss"]
  models$algorithm <- "PDPA"
  fpop.fit$loss$algorithm <- "FPOP"
  ggplot()+
    geom_abline(aes(slope=peaks, intercept=PoissonLoss, color=peaks),
                data=pdpa.fit$loss)+
    geom_text(aes(0, PoissonLoss, label=peaks),
              hjust=1,
              data=pdpa.fit$loss)+
    geom_point(aes(penalty.constant, penalized.loss, fill=algorithm),
               shape=21,
               data=fpop.fit$loss)+
    geom_point(aes(min.lambda, min.lambda*model.complexity + PoissonLoss,
                   fill=algorithm),
               shape=21,
               data=models)
  
})

