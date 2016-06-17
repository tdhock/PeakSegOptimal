PeakSegPDPA <- function
### Compute the PeakSeg constrained, Poisson loss, Segment Neighborhood
### model using a constrained version of the Pruned Dynamic
### Programming Algorithm.
(count.vec,
### integer vector of count data.
 weight.vec=rep(1, length(count.vec)),
### numeric vector (same length as count.vec) of positive weights.
 max.segments=NULL
### integer of length 1: maximum number of segments.
){
  n.data <- length(count.vec)
  stopifnot(is.integer(count.vec))
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
    "PeakSegPDPA_interface",
    count.vec=as.double(count.vec),
    weight.vec=as.double(weight.vec),
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
}

PeakSegPDPAchrom <- function
### Helper function which returns a list of data.frames
(input.dt,
### data.table with columns count, chromStart, chromEnd.
 max.segments=NULL
### maximum number of segments, or NULL which means to keep going
### until an active equality constraint is found.
){
  stop("TODO")
}
