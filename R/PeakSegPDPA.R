PeakSegPDPA <- structure(function
### Compute the PeakSeg constrained, Poisson loss, Segment Neighborhood
### model using a constrained version of the Pruned Dynamic
### Programming Algorithm. 
(count.vec,
### integer vector of count data.
 weight.vec=rep(1, length(count.vec)),
### integer vector (same length as count.vec) of positive weights.
 max.segments=NULL
### integer of length 1: maximum number of segments (must be >= 2).
){
  n.data <- length(count.vec)
  stopifnot(is.integer(count.vec))
  stopifnot(is.integer(weight.vec))
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
    "PeakSegPDPA_interface",
    count.vec=as.integer(count.vec),
    weight.vec=as.integer(weight.vec),
    n.data=as.integer(n.data),
    max.segments=as.integer(max.segments),
    cost.mat=as.double(cost.mat),
    ends.mat=as.integer(ends.mat),
    mean.mat=as.double(mean.mat),
    intervals.mat=as.integer(intervals.mat),
    PACKAGE="coseg")
  result.list$cost.mat <- matrix(
    result.list$cost.mat, max.segments, n.data, byrow=TRUE)
  result.list$ends.mat <- matrix(
    result.list$ends.mat+1L, max.segments, max.segments, byrow=TRUE)
  result.list$mean.mat <- matrix(
    result.list$mean.mat, max.segments, max.segments, byrow=TRUE)
  result.list$intervals.mat <- matrix(
    result.list$intervals.mat, max.segments, n.data, byrow=TRUE)
  result.list
}, ex=function(){
  data("H3K4me3_XJ_immune_chunk1")
  by.sample <-
    split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
  n.data.vec <- sapply(by.sample, nrow)
  one <- by.sample[[1]]
  count.vec <- one$coverage
  weight.vec <- with(one, chromEnd-chromStart)
  max.segments <- 19L
  fit <- PeakSegPDPA(count.vec, weight.vec, max.segments)
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

PeakSegPDPALog <- structure(function
### Compute the PeakSeg constrained, Poisson loss, Segment Neighborhood
### model using a constrained version of the Pruned Dynamic
### Programming Algorithm. 
(count.vec,
### integer vector of count data.
 weight.vec=rep(1, length(count.vec)),
### integer vector (same length as count.vec) of positive weights.
 max.segments=NULL
### integer of length 1: maximum number of segments (must be >= 2).
){
  n.data <- length(count.vec)
  stopifnot(is.integer(count.vec))
  stopifnot(is.integer(weight.vec))
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
    weight.vec=as.integer(weight.vec),
    n.data=as.integer(n.data),
    max.segments=as.integer(max.segments),
    cost.mat=as.double(cost.mat),
    ends.mat=as.integer(ends.mat),
    mean.mat=as.double(mean.mat),
    intervals.mat=as.integer(intervals.mat),
    PACKAGE="coseg")
  result.list$cost.mat <- matrix(
    result.list$cost.mat, max.segments, n.data, byrow=TRUE)
  result.list$ends.mat <- matrix(
    result.list$ends.mat+1L, max.segments, max.segments, byrow=TRUE)
  result.list$mean.mat <- matrix(
    result.list$mean.mat, max.segments, max.segments, byrow=TRUE)
  result.list$intervals.mat <- matrix(
    result.list$intervals.mat, max.segments, n.data, byrow=TRUE)
  result.list
}, ex=function(){
  data("H3K4me3_XJ_immune_chunk1")
  by.sample <-
    split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
  n.data.vec <- sapply(by.sample, nrow)
  one <- by.sample[[1]]
  count.vec <- one$coverage
  weight.vec <- with(one, chromEnd-chromStart)
  max.segments <- 19L
  fit <- PeakSegPDPA(count.vec, weight.vec, max.segments)
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
### Compute the PeakSeg constrained, Poisson loss, Segment Neighborhood
### model using a constrained version of the Pruned Dynamic
### Programming Algorithm.
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
    break.vec <- fit$ends.mat[n.segments, 2:n.segments]
    first <- c(1, break.vec+1)
    last <- c(break.vec, nrow(count.df))
    segments.list[[paste(n.segments)]] <- data.frame(
      mean=fit$mean.mat[n.segments, 1:n.segments],
      first,
      last,
      chromStart=count.df$chromStart[first],
      chromEnd=count.df$chromEnd[last],
      status=rep(c("background", "peak"), l=n.segments),
      peaks=(n.segments-1)/2,
      segments=n.segments)
  }
  is.feasible <- function(loss.vec){
    !any(diff(loss.vec) == 0, na.rm=TRUE)
  }
  list(
    segments=do.call(rbind, segments.list),
    loss=data.frame(
      segments=seg.vec,
      peaks=(seg.vec-1)/2,
      PoissonLoss=fit$cost.mat[seg.vec, length(data.vec)],
      feasible=apply(fit$mean.mat[seg.vec,], 1, is.feasible)))
### List of data.frames. 
}, ex=function(){

  ## samples for which pdpa recovers a more likely model, but it is
  ## not feasible for the PeakSeg problem (some segment means are
  ## equal).
  sample.id <- "McGill0322"
  sample.id <- "McGill0079"
  sample.id <- "McGill0106"
  n.peaks <- 3
  data("H3K4me3_XJ_immune_chunk1")
  H3K4me3_XJ_immune_chunk1$count <- H3K4me3_XJ_immune_chunk1$coverage
  by.sample <-
    split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
  one.sample <- by.sample[[sample.id]]
  pdpa.fit <- PeakSegPDPAchrom(one.sample, 9L)
  pdpa.segs <- subset(pdpa.fit$segments, n.peaks == peaks)
  dp.fit <- PeakSegDP(one.sample, 9L)
  dp.segs <- subset(dp.fit$segments, n.peaks == peaks)
  both.segs <- rbind(
    data.frame(dp.segs, algorithm="cDPA"),
    data.frame(pdpa.segs, algorithm="PDPA"))
  both.breaks <- subset(both.segs, 1 < first)

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
  dp.fit <- PeakSegDP(one.sample, 9L)

  ggplot()+
    scale_size_manual(values=c(cDPA=3, PDPA=1))+
    geom_point(aes(peaks, error,
                   size=algorithm, color=algorithm),
               data=data.frame(dp.fit$error, algorithm="cDPA"))+
    geom_point(aes(peaks, PoissonLoss,
                   size=algorithm, color=algorithm),
               data=data.frame(pdpa.fit$loss, algorithm="PDPA"))
  
})
