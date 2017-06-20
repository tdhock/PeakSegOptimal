IsotonicFPOP <- structure(function
                      (dat.vec,
                          ### double vector of length >= 3: real data to segment.
                          weight.vec=NULL,
                          ### numeric vector (same length as count.vec) of positive weights.
                          penalty=NULL,
                          ### non-negative numeric scalar: penalty
                          constraint=TRUE
                          ## positivity constraint?
                         ){
                           
                          if (!is.null(weight.vec)) {
                            warning("Weights not currently implemented") 
                          }
                           weight.vec <- rep(1, length(dat.vec))
                           n.data <- length(dat.vec)
                           stopifnot(3 <= n.data)
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
                           constraint <- constraint
                           
                           result.list <- .C(
                             "IsotonicFPOP_interface",
                             dat.vec=as.numeric(dat.vec),
                             n.data=as.integer(n.data),
                             penalty=as.numeric(penalty),
                             cost.mat=as.double(cost.mat),
                             ends.vec=as.integer(ends.vec),
                             mean.vec=as.double(mean.vec),
                             intervals.mat=as.integer(intervals.mat),
                             constraint = constraint,
                             PACKAGE="coseg")
                           
                           ## 1-indexed segment ends!
                           result.list$ends.vec <- result.list$ends.vec+1L
                           result.list$cost.mat <- matrix(
                             result.list$cost.mat*cumsum(weight.vec), 1, n.data, byrow=TRUE)
                           result.list$intervals.mat <- matrix(
                             result.list$intervals.mat, 1, n.data, byrow=TRUE)
                           
                           change.vec <- with(result.list, rev(ends.vec[ends.vec>0]))
                           change.sign.vec <- rep(c(1, -1), length(change.vec)/2)
                           end.vec <- c(change.vec, n.data)
                           start.vec <- c(1, change.vec+1)
                           length.vec <- end.vec-start.vec+1
                           mean.vec <- rev(result.list$mean.vec[1:(length(change.vec)+1)])
                           result.list$fitted.values <- rep(mean.vec, length.vec)
                           
                           return(result.list)
                           ### List of model parameters. count.vec, weight.vec, n.data, penalty
                           ### (input parameters), cost.mat (optimal Poisson loss), ends.vec
                           ### (optimal position of segment ends, 1-indexed), mean.vec (optimal
                           ### segment means), intervals.mat (number of intervals stored by the
                           ### functional pruning algorithm). 
                           ## fitted.values is the M vector 
                         })