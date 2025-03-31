\name{UnconstrainedFPOP}
\alias{UnconstrainedFPOP}
\title{UnconstrainedFPOP}
\description{Find the optimal change-points using the Poisson loss and no
constraints. For N data points, the functional pruning algorithm
is O(N log N) time and memory. It recovers the exact solution to
the following optimization problem. Let Z be an N-vector of count
data (\code{count.vec}, non-negative integers), let W be an N-vector of
positive weights (\code{weight.vec}), and let \code{penalty} be a non-negative
real number. Find the N-vector M of real numbers (segment means)
and (N-1)-vector C of change-point indicators in 0,1 which
minimize the penalized Poisson Loss, penalty*sum_[i=1]^[N_1]
I(c_i=1) + sum_[i=1]^N w_i*[m_i-z_i*log(m_i)].}
\usage{UnconstrainedFPOP(count.vec, 
    weight.vec = rep(1, 
        length(count.vec)), 
    penalty = NULL, verbose_file = "")}
\arguments{
  \item{count.vec}{integer vector of length >= 3: non-negative count data to segment.}
  \item{weight.vec}{numeric vector (same length as \code{count.vec}) of positive weights.}
  \item{penalty}{non-negative numeric scalar: \code{penalty} parameter (smaller for more
peaks, larger for fewer peaks).}
  \item{verbose_file}{String: file name to output candidate indices (slow but useful for
analysis of the algorithm). Default "" means to not output indices
(fast).}
}

\value{List of model parameters. \code{count.vec}, \code{weight.vec}, n.data, \code{penalty}
(input parameters), cost.vec (optimal Poisson loss), ends.vec
(optimal position of segment ends, 1-indexed), mean.vec (optimal
segment means), intervals.vec (number of intervals stored by the
functional pruning algorithm). To recover the solution in terms of
(M,C) variables, see the example.}

\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre]}




\examples{

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

}
