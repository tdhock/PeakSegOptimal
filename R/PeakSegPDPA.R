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

PeakSegPDPAchrom <- function
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
  for(n.segments in seq(1, max.segments, by=2)){
    break.vec <- fit$ends.mat[n.segments, 2:n.segments]
    first <- c(1, break.vec+1)
    last <- c(break.vec, n.data)
    pdpa.segments.list[[paste(n.segments)]] <- data.table(
      mean=fit$mean.mat[n.segments, 1:n.segments],
      first,
      last,
      chromStart=count.df$chromStart[first],
      chromEnd=count.df$chromEnd[last],
      status=rep(c("background", "peak"), l=n.segments),
      peaks=(n.segments-1)/2,
      segments=n.segments)
  }
  list(segments=do.call(rbind, segments.list))
### List of data.frames.
}
