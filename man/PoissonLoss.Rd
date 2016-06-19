\name{PoissonLoss}
\alias{PoissonLoss}
\title{PoissonLoss}
\description{Compute the weighted Poisson loss function, which is seg.mean -
count * log(seg.mean). The edge case is when the mean is zero, in
which case the probability mass function takes a value of 1 when
the data is 0 (and 0 otherwise). Thus the log-likelihood of a
maximum likelihood segment with mean zero must be zero.}
\usage{PoissonLoss(count, seg.mean, weight = 1)}
\arguments{
  \item{count}{
}
  \item{seg.mean}{
}
  \item{weight}{
}
}



\author{Toby Dylan Hocking}




\examples{
PoissonLoss(1, 1)
PoissonLoss(0, 0)
PoissonLoss(1, 0)
PoissonLoss(0, 1)
}
