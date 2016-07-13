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
                  label=sprintf("%d peak%s optimal", peaks,
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
  data(chr11ChIPseq)
  one <- subset(chr11ChIPseq$coverage, sample.id=="McGill0322")
  fit <- PeakSegDP(one, 5L)
  ## First compute the optimal number of peaks function.
  exact.df <- with(fit$error, exactModelSelection(error, segments, peaks))
  ## Then compute the PeakError of these models with respect to the
  ## annotated regions.
  regions <- subset(chr11ChIPseq$regions, sample.id=="McGill0322")
  peak.list <- split(fit$segments, fit$segments$peaks)
  require(PeakError)
  all.error <- NULL
  error.regions <- NULL
  for(peaks.chr in names(peak.list)){
    peaks <- subset(peak.list[[peaks.chr]], status=="peak")
    error <- PeakErrorChrom(peaks, regions)
    error.regions <- rbind(error.regions, data.frame(peaks=peaks.chr, error))
    all.error <- rbind(all.error, {
      data.frame(peaks=as.integer(peaks.chr),
                 errors=with(error, sum(fp+fn)))
    })
  }
  ## plot the annotation error of the 6 models.
  ann.colors <- c(noPeaks = "#f6f4bf", peakStart = "#ffafaf", 
                  peakEnd = "#ff4c4c", peaks = "#a445ee")
  library(ggplot2)
  ggplot()+
    geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                      fill=annotation),
                  data=error.regions, alpha=1/2)+
    geom_step(aes(chromStart/1e3, count), data=one, color="grey40")+
    geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                      linetype=status),
                  data=error.regions, fill=NA, color="black", size=2)+
    scale_fill_manual(values=ann.colors)+
    scale_linetype_manual(values=c("false negative"=3,
                          "false positive"=1,
                          "correct"=0))+
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
  rownames(all.error) <- all.error$peaks
  exact.df$errors <-
    all.error[as.character(exact.df$model.complexity), "errors"]
  indices <- with(exact.df, {
    largestContinuousMinimum(errors, max.log.lambda-min.log.lambda)
  })
  ## The target interval (min.log.lambda, max.log.lambda) is the
  ## largest continuous interval such that the error is minimal.
  target.interval <- with(indices, {
    data.frame(min.log.lambda=exact.df$min.log.lambda[start],
               max.log.lambda=exact.df$max.log.lambda[end])
  })
  print(target.interval)
})

