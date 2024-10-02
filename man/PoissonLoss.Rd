\name{PoissonLoss}
\alias{PoissonLoss}
\title{PoissonLoss}
\description{Compute the weighted Poisson loss function, which is \code{seg.mean} -
\code{count} * log(\code{seg.mean}). The edge case is when the mean is zero, in
which case the probability mass function takes a value of 1 when
the data is 0 (and 0 otherwise). Thus the log-likelihood of a
maximum likelihood segment with mean zero must be zero.}
\usage{PoissonLoss(count, seg.mean, 
    weight = 1)}
\arguments{
  \item{count}{count }
  \item{seg.mean}{seg.mean }
  \item{weight}{weight }
}



\author{Toby Dylan Hocking <toby.hocking@r-project.org> [aut, cre]}




\examples{
PoissonLoss(1, 1)
PoissonLoss(0, 0)
PoissonLoss(1, 0)
PoissonLoss(0, 1)
}
