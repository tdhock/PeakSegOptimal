\name{largestContinuousMinimum}
\alias{largestContinuousMinimum}
\title{largestContinuousMinimum}
\description{Find the run of minimum cost with the largest size.}
\usage{largestContinuousMinimum(cost, size)}
\arguments{
  \item{cost}{
}
  \item{size}{
}
}



\author{Toby Dylan Hocking}




\examples{

data("H3K4me3_XJ_immune_chunk1")
sample.id <- "McGill0106"
by.sample <-
  split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
one.sample <- by.sample[[sample.id]]
count.vec <- one.sample$coverage
weight.vec <- with(one.sample, chromEnd-chromStart)
max.segments <- 19L
fit <- PeakSegPDPA(count.vec, weight.vec, max.segments)
segs.vec <- seq(1, max.segments, by=2)
cost.vec <- with(fit, cost.mat[segs.vec, n.data])
peaks.vec <- seq(0, by=1, l=length(cost.vec))
## Typically we take model.complexity to be the number of changes,
## so that the penalized cost is the same as in FPOP.
model.complexity <- segs.vec-1
## Calculate the exact path of breakpoints in the optimal number of
## peaks function.
exact.df <- exactModelSelection(cost.vec, model.complexity, peaks.vec)
## Say that we have computed an error function which takes a minimum
## for the models with 1 or 3 peaks.
exact.df$errors <- c(3, 2, 2, 1, 1, 1, 0, 0, 1)
indices <- with(exact.df, {
  largestContinuousMinimum(errors, max.log.lambda-min.log.lambda)
})
target.interval <- data.frame(
  min.log.lambda=exact.df$min.log.lambda[indices$start],
  max.log.lambda=exact.df$max.log.lambda[indices$end],
  errors=exact.df$errors[indices$start])
library(ggplot2)
ggplot()+
  ggtitle(
    "target interval (red) is the set of penalties with min error (black)")+
  geom_segment(aes(min.log.lambda, errors,
                   xend=max.log.lambda, yend=errors),
               data=target.interval,
               color="red",
               size=2)+
  geom_segment(aes(min.log.lambda, errors,
                   xend=max.log.lambda, yend=errors),
               data=exact.df)+
  ylab("errors of selected model")+
  xlab("penalty constant log(lambda)")

}
