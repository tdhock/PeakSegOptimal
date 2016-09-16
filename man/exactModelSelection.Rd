\name{exactModelSelection}
\alias{exactModelSelection}
\title{Exact model selection function}
\description{Given a set of optimal costs C_i, and model complexity values K_i,
and a model selection function i*(lambda) = argmin_i C_i +
lambda*K_i, compute a set of consecutive (K_i, min.lambda,
max.lambda) with i being the solution for every lambda in
(min.lambda, max.lambda).}
\usage{exactModelSelection(cost, model.complexity, peaks)}
\arguments{
  \item{cost}{numeric vector: optimal costs C_i.}
  \item{model.complexity}{numeric vector: model complexity K_i.}
  \item{peaks}{
}
}



\author{Toby Dylan Hocking}




\examples{

data("H3K4me3_XJ_immune_chunk1", envir=environment())
sample.id <- "McGill0106"
by.sample <-
  split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
one.sample <- by.sample[[sample.id]]
count.vec <- one.sample$coverage
weight.vec <- with(one.sample, chromEnd-chromStart)
max.segments <- 5L
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
exact.df$cost <- rev(cost.vec) + exact.df$min.lambda * exact.df$model.complexity
exact.df$next.cost <- c(exact.df$cost[-1], NA)
exact.df$PoissonLoss <- rev(cost.vec)
library(ggplot2)
ggplot()+
  geom_point(aes(min.lambda, cost),
             data=exact.df, pch=1, color="red")+
  geom_segment(aes(min.lambda, cost,
                   xend=max.lambda, yend=next.cost),
               data=exact.df, color="red", size=1.5)+
  geom_text(aes((min.lambda+max.lambda)/2, (cost+next.cost)/2,
                label=sprintf("\%d peak\%s optimal", peaks,
                  ifelse(peaks==1, "", "s"))),
            data=exact.df, color="red", hjust=0, vjust=1.5)+
  geom_abline(aes(slope=model.complexity, intercept=PoissonLoss),
              data=exact.df)+
  geom_text(aes(0, PoissonLoss, label=peaks),
            data=exact.df, hjust=1.5, color="red")+
  ggtitle("model selection: cost = PoissonLoss_k + lambda*changes_k")

## Solve the optimization using grid search.
L.grid <- with(exact.df,{
  seq(min(max.log.lambda)-1,
      max(min.log.lambda)+1,
      l=100)
})
lambda.grid <- exp(L.grid)
kstar.grid <- sapply(lambda.grid, function(lambda){
  crit <- with(exact.df, model.complexity * lambda + PoissonLoss)
  picked <- which.min(crit)
  exact.df$peaks[picked]
})
grid.df <- data.frame(log.lambda=L.grid, peaks=kstar.grid)
## Compare the results.
ggplot()+
  ggtitle("grid search (red) agrees with exact path computation (black)")+
  geom_segment(aes(min.log.lambda, peaks,
                   xend=max.log.lambda, yend=peaks),
               data=exact.df)+
  geom_point(aes(log.lambda, peaks),
             data=grid.df, color="red", pch=1)+
  ylab("optimal model complexity (peaks)")+
  xlab("log(lambda)")

}
