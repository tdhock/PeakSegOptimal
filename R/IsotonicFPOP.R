IsotonicFPOP <- structure(function
                      (count.vec,
                          ### double vector of length >= 3: real data to segment.
                          weight.vec=rep(1, length(count.vec)),
                          ### numeric vector (same length as count.vec) of positive weights.
                          penalty=NULL
                          ### non-negative numeric scalar: penalty parameter (smaller for more
                          ### peaks, larger for fewer peaks).
                         ){
                           n.data <- length(count.vec)
                           stopifnot(3 <= n.data)
                           stopifnot(0 <= count.vec)
                           stopifnot(is.numeric(weight.vec))
                           stopifnot(n.data==length(weight.vec))
                           stopifnot(0 < weight.vec)
                           stopifnot(is.numeric(penalty))
                           stopifnot(length(penalty)==1)
                           stopifnot(0 <= penalty)
                           cost.mat <- double(n.data)
                           ends.vec <- integer(n.data)
                           mean.vec <- double(n.data)
                           intervals.mat <- integer(n.data)
                           
                           # cat(paste0("n.data = \t", n.data, "\n"))
                           # cat(paste0("count.vec = \t", count.vec, "\n"))
                           # cat(paste0("weight.vec = \t", weight.vec, "\n"))
                           # cat(paste0("penalty = \t", penalty, "\n"))
                           # cat(paste0("cost.mat = \t", cost.mat, "\n"))
                           # cat(paste0("ends.vec = \t", ends.vec, "\n"))
                           # cat(paste0("intervals.mat = \t", intervals.mat, "\n"))
                           
                           result.list <- .C(
                             "IsotonicFPOP_interface",
                             count.vec=as.numeric(count.vec),
                             n.data=as.integer(n.data),
                             penalty=as.numeric(penalty),
                             cost.mat=as.double(cost.mat),
                             ends.vec=as.integer(ends.vec),
                             mean.vec=as.double(mean.vec),
                             intervals.mat=as.integer(intervals.mat),
                             PACKAGE="coseg")
                           ## 1-indexed segment ends!
                           result.list$ends.vec <- result.list$ends.vec+1L
                           result.list$cost.mat <- matrix(
                             result.list$cost.mat*cumsum(weight.vec), 1, n.data, byrow=TRUE)
                           result.list$intervals.mat <- matrix(
                             result.list$intervals.mat, 1, n.data, byrow=TRUE)
                           result.list
                           ### List of model parameters. count.vec, weight.vec, n.data, penalty
                           ### (input parameters), cost.mat (optimal Poisson loss), ends.vec
                           ### (optimal position of segment ends, 1-indexed), mean.vec (optimal
                           ### segment means), intervals.mat (number of intervals stored by the
                           ### functional pruning algorithm). To recover the solution in terms of
                           ### (M,C) variables, see the example.
                         })