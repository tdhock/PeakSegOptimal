\name{PeakSegPDPA}
\alias{PeakSegPDPA}
\title{PeakSegPDPA}
\description{Compute the PeakSeg constrained, Poisson loss, Segment Neighborhood
model using a constrained version of the Pruned Dynamic
Programming Algorithm.}
\usage{PeakSegPDPA(count.vec, weight.vec = rep(1, length(count.vec)), 
    max.segments = NULL)}
\arguments{
  \item{count.vec}{integer vector of count data.}
  \item{weight.vec}{numeric vector (same length as count.vec) of positive weights.}
  \item{max.segments}{integer of length 1: maximum number of segments.}
}



\author{Toby Dylan Hocking}




\examples{
data("H3K4me3_XJ_immune_chunk1")
by.sample <-
  split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
n.data.vec <- sapply(by.sample, nrow)
one <- by.sample[[1]]
count.vec <- one$coverage
weight.vec <- with(one, chromEnd-chromStart)
max.segments <- 19L
fit <- PeakSegPDPA(count.vec, weight.vec, max.segments)
PDPA.intervals <- data.frame(
  segments=as.numeric(row(fit$intervals.mat)),
  data=as.numeric(col(fit$intervals.mat)),
  intervals=as.numeric(fit$intervals.mat))
some.intervals <- subset(PDPA.intervals, segments<data & 1<segments)
library(ggplot2)
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(segments ~ .)+
  geom_line(aes(data, intervals), data=some.intervals)+
  scale_y_continuous(
    "intervals stored by the\nconstrained optimal segmentation algorithm",
    breaks=c(20, 40))
}
