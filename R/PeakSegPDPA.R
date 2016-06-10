PeakSegPDPA <- function
### Compute the PeakSeg constrained, Poisson loss, Segment Neighborhood
### model using a constrained version of the Pruned Dynamic
### Programming Algorithm.
(count.vec,
### integer vector of count data.
 weight.vec,
### numeric vector of positive weights.
 max.segments=NULL
### maximum number of segments, or NULL which means to keep going
### until an active equality constraint is found.
){
  stop("TODO write C++ code")
  return("full cost matrix for comparison with cDPA function")
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
