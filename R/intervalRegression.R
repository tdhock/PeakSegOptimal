exactModelSelection <- structure(function # Exact model selection function
### Given a set of optimal costs C_i, and model complexity values K_i,
### and a model selection function i*(lambda) = argmin_i C_i +
### lambda*K_i, compute a set of consecutive (K_i, min.lambda,
### max.lambda) with i being the solution for every lambda in
### (min.lambda, max.lambda).
(cost,
### numeric vector: optimal costs C_i.
 model.complexity,
### numeric vector: model complexity K_i.
 peaks){
  stopifnot(is.numeric(cost))
  stopifnot(is.numeric(model.complexity))
  stopifnot(diff(model.complexity) > 0)
  stopifnot(diff(cost) < 0)
  stopifnot(length(cost) == length(model.complexity))
  n.models <- length(cost)
  Kmax <- model.complexity[n.models]
  Kcurrent <- Kmax
  Lcurrent <- 0
  vK <- Kmax
  vL <- 0
  vP <- peaks[n.models]
  i <- 2
  min.complexity <- model.complexity[1]
  while(Kcurrent > min.complexity) {
    is.smaller <- model.complexity < Kcurrent
    is.current <- model.complexity == Kcurrent
    smallerK <- model.complexity[is.smaller]
    smallerPeaks <- peaks[is.smaller]
    cost.term <- cost[is.current] - cost[is.smaller]
    complexity.term <- smallerK - model.complexity[is.current]
    lambdaTransition <- cost.term/complexity.term
    next.i <- which.min(lambdaTransition)
    Kcurrent <- smallerK[next.i]
    Lcurrent <- min(lambdaTransition)
    vL[i] <- Lcurrent
    vK[i] <- Kcurrent
    vP[i] <- smallerPeaks[next.i]
    i <- i + 1
  }
  L <- log(vL)
  data.frame(min.log.lambda = L,
             max.log.lambda = c(L[-1], Inf),
             model.complexity = vK,
             peaks=vP,
             min.lambda = vL,
             max.lambda = c(vL[-1], Inf))
},ex=function(){

  data("H3K4me3_XJ_immune_chunk1")
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
                  label=sprintf("%d peak%s optimal", peaks,
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
  
})

largestContinuousMinimum <- structure(function
### Find the run of minimum cost with the largest size.
(cost,
 size
 ){
  m <- min(cost)
  is.min <- cost == m
  d <- c(diff(c(FALSE,is.min,FALSE)))
  ##print(data.frame(cost=c(cost,NA),size=c(size,NA),diff=d))
  starts <- which(d==1)
  ends <- which(d==-1)-1
  runs <- data.frame(starts,ends)
  ##print(runs)
  runs$size <- sapply(seq_along(starts),function(i){
    sum(size[ starts[i]:ends[i] ])
  })
  ##print(runs)
  largest <- which.max(runs$size)
  list(start=starts[largest],end=ends[largest])
}, ex=function(){

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

})

