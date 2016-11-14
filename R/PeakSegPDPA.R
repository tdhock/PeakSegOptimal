PeakSegPDPA <- structure(function
### Find the optimal change-points using the Poisson loss and the
### PeakSeg constraint. For N data points and S segments, the
### functional pruning algorithm is O(S*NlogN) space and O(S*NlogN)
### time. It recovers the exact solution to the following optimization
### problem. Let Z be an N-vector of count data (count.vec,
### non-negative integers) and let W be an N-vector of positive
### weights (weight.vec). Find the N-vector M of real numbers (segment
### means) and (N-1)-vector C of change-point indicators in {-1,0,1}
### which minimize the Poisson Loss, sum_{i=1}^N
### w_i*[m_i-z_i*log(m_i)], subject to constraints: (1) there are
### exactly S-1 non-zero elements of C, and (2) the first change is up
### and the next change is down, etc (sum_{i=1}^t c_i in {0,1} for all
### t<N), and (3) Every zero-valued change-point variable has an equal
### segment mean after: c_i=0 implies m_i=m_{i+1}, (4) every
### positive-valued change-point variable may have an up change after:
### c_i=1 implies m_i<=m_{i+1}, (5) every negative-valued change-point
### variable may have a down change after: c_i=-1 implies
### m_i>=m_{i+1}. Note that when the equality constraints are active
### for non-zero change-point variables, the recovered model is not
### feasible for the strict inequality constraints of the PeakSeg
### problem, and the optimum of the PeakSeg problem is undefined.
(count.vec,
### integer vector of count data.
 weight.vec=rep(1, length(count.vec)),
### numeric vector (same length as count.vec) of positive weights.
 max.segments=NULL
### integer of length 1: maximum number of segments (must be >= 2).
){
  n.data <- length(count.vec)
  stopifnot(is.integer(count.vec))
  stopifnot(0 <= count.vec)
  stopifnot(is.numeric(weight.vec))
  stopifnot(n.data==length(weight.vec))
  stopifnot(0 < weight.vec)
  stopifnot(is.integer(max.segments))
  stopifnot(length(max.segments)==1)
  stopifnot(2 <= max.segments && max.segments <= n.data)
  cost.mat <- double(n.data*max.segments)
  ends.mat <- integer(max.segments*max.segments)
  mean.mat <- double(max.segments*max.segments)
  intervals.mat <- integer(n.data*max.segments)
  result.list <- .C(
    "PeakSegPDPALog_interface",
    count.vec=as.integer(count.vec),
    weight.vec=as.numeric(weight.vec),
    n.data=as.integer(n.data),
    max.segments=as.integer(max.segments),
    cost.mat=as.double(cost.mat),
    ends.mat=as.integer(ends.mat),
    mean.mat=as.double(mean.mat),
    intervals.mat=as.integer(intervals.mat),
    PACKAGE="coseg")
  result.list$cost.mat <- matrix(
    result.list$cost.mat*cumsum(weight.vec), max.segments, n.data, byrow=TRUE)
  result.list$ends.mat <- matrix(
    result.list$ends.mat+1L, max.segments, max.segments, byrow=TRUE)
  result.list$mean.mat <- matrix(
    result.list$mean.mat, max.segments, max.segments, byrow=TRUE)
  result.list$intervals.mat <- matrix(
    result.list$intervals.mat, max.segments, n.data, byrow=TRUE)
  result.list
### List of model parameters. count.vec, weight.vec, n.data,
### max.segments (input parameters), cost.mat (optimal Poisson loss),
### ends.mat (optimal position of segment ends, 1-indexed), mean.mat
### (optimal segment means), intervals.mat (number of intervals stored
### by the functional pruning algorithm). To recover the solution in
### terms of (M,C) variables, see the example.
}, ex=function(){

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
  stopifnot(cumsum(sign(C.vec)) %in% c(0, 1))

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
      "intervals stored by the\nconstrained optimal segmentation algorithm",
      breaks=c(20, 40))
  
})

PeakSegPDPAchrom <- structure(function
### Find the optimal change-points using the Poisson loss and the
### PeakSeg constraint. This function is a user-friendly interface to
### the PeakSegPDPA function.
(count.df,
### data.frame with columns count, chromStart, chromEnd.
 max.peaks=NULL
### integer > 0: maximum number of peaks.
){
  stopifnot(is.data.frame(count.df))
  stopifnot(is.integer(count.df$chromStart))
  stopifnot(is.integer(count.df$chromEnd))
  stopifnot(is.integer(count.df$count))
  stopifnot(is.integer(max.peaks))
  stopifnot(0 < max.peaks)
  stopifnot(count.df$chromStart < count.df$chromEnd)
  stopifnot(0 <= count.df$chromStart)
  stopifnot(0 <= count.df$count)
  weight.vec <- with(count.df, chromEnd - chromStart)
  data.vec <- count.df$count
  max.segments <- as.integer(max.peaks*2+1)
  stopifnot(max.segments <= nrow(count.df))
  fit <- PeakSegPDPA(data.vec, weight.vec, max.segments)
  segments.list <- list()
  seg.vec <- seq(1, max.segments, by=2)
  for(n.segments in seg.vec){
    break.vec <- if(n.segments==1){
      c()
    }else{
      fit$ends.mat[n.segments, 2:n.segments]
    }
    first <- c(1, break.vec+1)
    last <- c(break.vec, nrow(count.df))
    status.str <- rep(c("background", "peak"), l=n.segments)
    segments.list[[paste(n.segments)]] <- data.frame(
      mean=fit$mean.mat[n.segments, 1:n.segments],
      first,
      last,
      chromStart=count.df$chromStart[first],
      chromEnd=count.df$chromEnd[last],
      status=factor(status.str, c("background", "peak")),
      peaks=(n.segments-1)/2,
      segments=n.segments)
  }
  is.feasible <- function(loss.vec){
    !any(diff(loss.vec) == 0, na.rm=TRUE)
  }
  loss.df <- data.frame(
    segments=seg.vec,
    peaks=(seg.vec-1)/2,
    PoissonLoss=fit$cost.mat[seg.vec, length(data.vec)],
    feasible=apply(fit$mean.mat[seg.vec,], 1, is.feasible))
  rownames(loss.df) <- loss.df$peaks
  seg.df <- do.call(rbind, segments.list)
  only.feasible <- loss.df[loss.df$feasible,]
  rownames(seg.df) <- NULL
  cum.vec <- cummin(loss.df$PoissonLoss)
  min.loss <- loss.df[cum.vec==loss.df$PoissonLoss,]
  is.dec <- c(TRUE, diff(min.loss$PoissonLoss) < 0)
  dec.loss <- min.loss[is.dec,]
  dec.models <- with(dec.loss, exactModelSelection(PoissonLoss, peaks, peaks))
  list(
    segments=seg.df,
    loss=loss.df,
    modelSelection.feasible=with(only.feasible, {
      exactModelSelection(PoissonLoss, peaks, peaks)
    }),
    modelSelection.decreasing=dec.models)
### List of data.frames: segments can be used for plotting the
### segmentation model, loss describes model loss and feasibility,
### modelSelection.feasible describes the set of all linear penalty
### (lambda) values which can be used to select the feasible models,
### modelSelection.decreasing selects from all models that decrease
### the Poisson loss relative to simpler models (same as PeakSegFPOP).
}, ex=function(){

  ## samples for which pdpa recovers a more likely model, but it is
  ## not feasible for the PeakSeg problem (some segment means are
  ## equal).
  sample.id <- "McGill0322"
  sample.id <- "McGill0079"
  sample.id <- "McGill0106"
  n.peaks <- 3
  data("H3K4me3_XJ_immune_chunk1", envir=environment())
  H3K4me3_XJ_immune_chunk1$count <- H3K4me3_XJ_immune_chunk1$coverage
  by.sample <-
    split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
  one.sample <- by.sample[[sample.id]]
  pdpa.fit <- PeakSegPDPAchrom(one.sample, 9L)
  pdpa.segs <- subset(pdpa.fit$segments, n.peaks == peaks)
  both.segs.list <- list(pdpa=data.frame(pdpa.segs, algorithm="PDPA"))
  if(require(PeakSegDP)){
    dp.fit <- PeakSegDP(one.sample, 9L)
    dp.segs <- subset(dp.fit$segments, n.peaks == peaks)
    both.segs.list$dp <- data.frame(dp.segs, algorithm="cDPA")
  }
  both.segs <- do.call(rbind, both.segs.list)
  both.breaks <- subset(both.segs, 1 < first)
  library(ggplot2)
  ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(algorithm ~ ., scales="free")+
    geom_step(aes(chromStart/1e3, coverage),
              data=one.sample, color="grey")+
    geom_segment(aes(chromStart/1e3, mean,
                     xend=chromEnd/1e3, yend=mean),
                 color="green",
                 data=both.segs)+
    geom_vline(aes(xintercept=chromStart/1e3),
               color="green",
               linetype="dashed",
               data=both.breaks)

  ## samples for which pdpa recovers some feasible models that the
  ## heuristic dp does not.
  sample.id.vec <- c(
    "McGill0091", "McGill0107", "McGill0095",
    "McGill0059", "McGill0029", "McGill0010")
  sample.id <- sample.id.vec[1]
  one.sample <- by.sample[[sample.id]]
  pdpa.fit <- PeakSegPDPAchrom(one.sample, 9L)
  gg.loss <- ggplot()+
    scale_size_manual(values=c(cDPA=2, PDPA=3))+
    scale_fill_manual(values=c(cDPA="white", PDPA="black"))+
    geom_point(aes(peaks, PoissonLoss,
                   size=algorithm, fill=algorithm),
               shape=21,
               data=data.frame(pdpa.fit$loss, algorithm="PDPA"))
  if(require(PeakSegDP)){
    dp.fit <- PeakSegDP(one.sample, 9L)
    gg.loss <- gg.loss+
      geom_point(aes(peaks, error,
                     size=algorithm, fill=algorithm),
                 shape=21,
                 data=data.frame(dp.fit$error, algorithm="cDPA"))
  }
  gg.loss
  
})
