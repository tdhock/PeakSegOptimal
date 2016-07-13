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
data(chr11ChIPseq)
one <- subset(chr11ChIPseq$coverage, sample.id=="McGill0002")
fit <- PeakSegDP(one, 5L)
library(ggplot2)
ggplot()+
  geom_step(aes(chromStart/1e3, count), data=one)+
  geom_segment(aes(chromStart/1e3, mean,
                   xend=chromEnd/1e3, yend=mean),
               data=fit$segments, color="green")+
  geom_segment(aes(chromStart/1e3, 0,
                   xend=chromEnd/1e3, yend=0),
               data=subset(fit$segments, status=="peak"),
               size=3, color="deepskyblue")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(peaks ~ ., scales="free", labeller=function(var, val){
    s <- ifelse(val==1, "", "s")
    paste0(val, " peak", s)
  })
## Calculate the exact path of breakpoints in the optimal number of
## peaks function.
rownames(fit$error) <- fit$error$peaks
exact.df <- with(fit$error, exactModelSelection(error, segments, peaks))
intercept <- fit$error[as.character(exact.df$peaks), "error"]
exact.df$cost <- intercept + exact.df$min.lambda * exact.df$model.complexity
exact.df$next.cost <- c(exact.df$cost[-1], NA)
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
  geom_abline(aes(slope=segments, intercept=error), data=fit$error)+
  geom_text(aes(0, error, label=peaks),
            data=fit$error, hjust=1.5, color="red")+
  ggtitle("model selection: cost = loss_k + lambda*segments_k")
## Solve the optimization using grid search.
L.grid <- with(exact.df,{
  seq(min(max.log.lambda)-1,
      max(min.log.lambda)+1,
      l=100)
})
lambda.grid <- exp(L.grid)
kstar.grid <- sapply(lambda.grid,function(lambda){
  crit <- fit$error$segments * lambda + fit$error$error
  picked <- which.min(crit)
  fit$error$peaks[picked]
})
grid.df <- data.frame(log.lambda=L.grid, peaks=kstar.grid)
## Compare the results.
ggplot()+
  geom_segment(aes(min.log.lambda, peaks,
                   xend=max.log.lambda, yend=peaks),
               data=exact.df)+
  geom_point(aes(log.lambda, peaks),
             data=grid.df, color="red", pch=1)+
  ylab("optimal model complexity (peaks)")+
  xlab("log(lambda)")
}
