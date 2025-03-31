UnconstrainedFPOP <- structure(function
### Find the optimal change-points using the Poisson loss and no
### constraints. For N data points, the functional pruning algorithm
### is O(N log N) time and memory. It recovers the exact solution to
### the following optimization problem. Let Z be an N-vector of count
### data (count.vec, non-negative integers), let W be an N-vector of
### positive weights (weight.vec), and let penalty be a non-negative
### real number. Find the N-vector M of real numbers (segment means)
### and (N-1)-vector C of change-point indicators in 0,1 which
### minimize the penalized Poisson Loss, penalty*sum_[i=1]^[N_1]
### I(c_i=1) + sum_[i=1]^N w_i*[m_i-z_i*log(m_i)].
(count.vec,
### integer vector of length >= 3: non-negative count data to segment.
 weight.vec=rep(1, length(count.vec)),
### numeric vector (same length as count.vec) of positive weights.
 penalty=NULL,
### non-negative numeric scalar: penalty parameter (smaller for more
### peaks, larger for fewer peaks).
  verbose_file=""
### String: file name to output candidate indices (slow but useful for
### analysis of the algorithm). Default "" means to not output indices
### (fast).
){
  n.data <- length(count.vec)
  stopifnot(1 <= n.data)
  stopifnot(is.integer(count.vec))
  stopifnot(0 <= count.vec)
  stopifnot(is.numeric(weight.vec))
  stopifnot(n.data==length(weight.vec))
  stopifnot(0 < weight.vec)
  stopifnot(is.numeric(penalty))
  stopifnot(length(penalty)==1)
  stopifnot(0 <= penalty)
  cost.vec <- double(n.data)
  ends.vec <- integer(n.data)
  mean.vec <- double(n.data)
  intervals.vec <- integer(n.data)
  result.list <- .C(
    "UnconstrainedFPOP_interface",
    count.vec=as.integer(count.vec),
    weight.vec=as.numeric(weight.vec),
    n.data=as.integer(n.data),
    penalty=as.numeric(penalty),
    verbose_file=as.character(verbose_file),
    cost.vec=as.double(cost.vec),
    ends.vec=as.integer(ends.vec),
    mean.vec=as.double(mean.vec),
    intervals.vec=as.integer(intervals.vec),
    PACKAGE="PeakSegOptimal")
  ## 1-indexed segment ends!
  result.list$ends.vec <- result.list$ends.vec+1L
  if(file.exists(verbose_file)){
    result.list$index_dt <- data.table::fread(verbose_file)
  }
  result.list
### List of model parameters. count.vec, weight.vec, n.data, penalty
### (input parameters), cost.vec (optimal Poisson loss), ends.vec
### (optimal position of segment ends, 1-indexed), mean.vec (optimal
### segment means), intervals.vec (number of intervals stored by the
### functional pruning algorithm). To recover the solution in terms of
### (M,C) variables, see the example.
}, ex=function(){

  mean_vec <- c(10,20,5,25)
  N_per_seg <- 10
  data_mean_vec <- rep(mean_vec, each=N_per_seg)
  N_data <- length(data_mean_vec)
  set.seed(1)
  data_value <- rpois(N_data, data_mean_vec)
  fit <- UnconstrainedFPOP(data_value, penalty=10, verbose_file = tempfile())

  if(require(ggplot2)){
    ggplot()+
      geom_tile(aes(
        data_index, change, fill=cost),
        data=fit$index_dt)+
      scale_fill_gradient(low="white", high="red")
  }

})

