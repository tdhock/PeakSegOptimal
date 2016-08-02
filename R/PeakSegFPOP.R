PeakSegFPOP <- structure(function
### Find the optimal change-points using the Poisson loss and the
### PeakSeg constraint. For N data points, the functional pruning
### algorithm is O(N) space and O(NlogN) time. It recovers the exact
### solution to the following optimization problem. Let Z be an
### N-vector of count data (non-negative integers). Find the N-vector
### M of real numbers (segment means) which minimize the penalized
### Poisson Loss, sum_{i=2}^N I(m_i != m_{i-1})*penalty + sum_{i=1}^N
### m_i - z_i * log(m_i), subject to constraint: up changes are
### followed by down changes, and vice versa. Note that the segment
### means can be equal, in which case the recovered model is not
### feasible for the PeakSeg problem. Unlike PeakSegPDPA which forces
### the first segment mean to be down (mu1 <= mu2), PeakSegFPOP may
### recover a model with the first segment mean up (mu1 >= mu2).
(count.vec,
### integer vector of count data.
 weight.vec=rep(1, length(count.vec)),
### numeric vector (same length as count.vec) of positive weights.
 penalty=NULL
### numeric of length 1: penalty parameter (smaller for more peaks,
### larger for fewer peaks).
){
  n.data <- length(count.vec)
  stopifnot(1 < n.data)
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
    ##label.vec=as.integer(label.vec),
    PACKAGE="coseg")
  ## 1-indexed segment ends!
  result.list$ends.vec <- result.list$ends.vec+1L
  result.list$cost.mat <- matrix(
    result.list$cost.mat, 2, n.data, byrow=TRUE)
  result.list$intervals.mat <- matrix(
    result.list$intervals.mat, 2, n.data, byrow=TRUE)
  result.list
### List of model parameters. count.vec, weight.vec, n.data, penalty
### (input parameters), cost.mat (optimal Poisson loss), ends.vec
### (optimal position of segment ends, 1-indexed), mean.vec (optimal
### segment means), intervals.mat (number of intervals stored by the
### functional pruning algorithm), label.vec (1=up or 0=down).
}, ex=function(){

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

})

PeakSegFPOPchrom <- structure(function
### Find the optimal change-points using the Poisson loss and the
### PeakSeg constraint. This function is a user-friendly interface to
### the PeakSegFPOP function.
(count.df,
### data.frame with columns count, chromStart, chromEnd.
 penalty=NULL
### integer > 0: maximum number of peaks.
){
  stopifnot(is.data.frame(count.df))
  n.data <- nrow(count.df)
  stopifnot(1 < n.data)
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
  data("H3K4me3_XJ_immune_chunk1")
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

