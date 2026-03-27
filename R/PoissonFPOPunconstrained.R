PoissonFPOPunconstrained <- structure(function
### Unconstrained optimal changepoint detection using the Poisson loss
### and the FPOP algorithm. Unlike PeakSegFPOP, there is no up-down
### constraint on segment means.
(count.vec,
### integer vector of non-negative count data to segment.
 weight.vec=rep(1, length(count.vec)),
### numeric vector (same length as count.vec) of positive weights.
 penalty=NULL
### non-negative numeric scalar: penalty for each changepoint.
){
  n.data <- length(count.vec)
  stopifnot(2 <= n.data)
  stopifnot(is.integer(count.vec))
  stopifnot(0 <= count.vec)
  stopifnot(is.numeric(weight.vec))
  stopifnot(n.data == length(weight.vec))
  stopifnot(0 < weight.vec)
  stopifnot(is.numeric(penalty))
  stopifnot(length(penalty) == 1)
  stopifnot(0 <= penalty)
  cost.vec <- double(n.data)
  end.vec <- integer(n.data)
  mean.vec <- double(n.data)
  intervals.vec <- integer(n.data)
  result.list <- .C(
    "PoissonFPOPunconstrained_interface",
    count.vec=as.integer(count.vec),
    weight.vec=as.numeric(weight.vec),
    n.data=as.integer(n.data),
    penalty=as.numeric(penalty),
    cost.vec=as.double(cost.vec),
    end.vec=as.integer(end.vec),
    mean.vec=as.double(mean.vec),
    intervals.vec=as.integer(intervals.vec),
    PACKAGE="PeakSegOptimal")
  n.segments <- sum(result.list$end.vec != -2L)
  if(n.segments == 0) n.segments <- 1L
  seg.mean <- rev(result.list$mean.vec[1:n.segments])
  seg.end <- rev(result.list$end.vec[1:n.segments])
  result.list$seg.mean <- seg.mean
  result.list$seg.end <- seg.end + 1L
  result.list$n.segments <- n.segments
  result.list
### List with seg.mean (chronological segment means), seg.end
### (1-indexed segment ends), n.segments, cost.vec, end.vec,
### mean.vec, intervals.vec.
})
