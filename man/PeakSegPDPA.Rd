\name{PeakSegPDPA}
\alias{PeakSegPDPA}
\title{PeakSegPDPA}
\description{Compute the PeakSeg constrained, Poisson loss, Segment Neighborhood
model using a constrained version of the Pruned Dynamic
Programming Algorithm.}
\usage{PeakSegPDPA(count.vec, weight.vec, max.segments = NULL)}
\arguments{
  \item{count.vec}{integer vector of count data.}
  \item{weight.vec}{numeric vector of positive weights.}
  \item{max.segments}{maximum number of segments, or NULL which means to keep going
until an active equality constraint is found.}
}



\author{Toby Dylan Hocking}





