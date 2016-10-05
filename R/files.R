problem.coverage <- function
### Ensure that coverage.bedGraph has been correctly computed for a
### particular genomic segmentation problem.
(problem.dir
### Path to a directory like sampleID/problems/problemID where
### sampleID/coverage.bedGraph contains counts of aligned reads in the
### entire genome, and sampleID/problems/problemID/problem.bed
### contains one line that indicates the genomic coordinates of a
### particular segmentation problem.
 ){
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  problems.dir <- dirname(problem.dir)
  sample.dir <- dirname(problems.dir)
  problem.bed <- file.path(problem.dir, "problem.bed")
  problem <- fread(problem.bed)
  setnames(problem, c("chrom", "problemStart", "problemEnd"))
  ## First check if problem/coverage.bedGraph has been created.
  prob.cov.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
  coverage.ok <- tryCatch({
    tail.cmd <- paste("tail -1", prob.cov.bedGraph)
    last.cov <- fread(tail.cmd)
    setnames(last.cov, c("chrom", "chromStart", "chromEnd", "coverage"))
    last.cov$chromEnd == problem$problemEnd
  }, error=function(e){
    FALSE
  })
  ## Create problemID/coverage.bedGraph if it is not present.
  if(!coverage.ok){
    ## Use intersectBed -sorted to avoid memory
    ## problems. coverage.bedGraph needs to be -a since that is
    ## reported in the output.
    coverage.bedGraph <- file.path(sample.dir, "coverage.bedGraph")
    intersect.cmd <- paste(
      "intersectBed -sorted",
      "-a", coverage.bedGraph,
      "-b", problem.bed,
      ">", prob.cov.bedGraph)
    cat(intersect.cmd, "\n")
    system(intersect.cmd)
  }
### Nothing. If necessary, the intersectBed command line program is
### used to create problemID/coverage.bedGraph, but it is not read
### into memory.
}

problem.PeakSegFPOP <- function
### Run PeakSegFPOP on one genomic segmentation problem directory.
(problem.dir,
### Path to a directory like sampleID/problems/problemID which
### contains a coverage.bedGraph file with the aligned read counts for
### one genomic segmentation problem.
 penalty.str
### Penalty parameter to pass to the PeakSegFPOP command line program.
 ){
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  stopifnot(is.character(penalty.str))
  stopifnot(length(penalty.str)==1)
  prob.cov.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
  pre <- paste0(prob.cov.bedGraph, "_penalty=", penalty.str)
  penalty_segments.bed <- paste0(pre, "_segments.bed")
  penalty_loss.tsv <- paste0(pre, "_loss.tsv")
  penalty_timing.tsv <- paste0(pre, "_timing.tsv")
  already.computed <- tryCatch({
    penalty.segs <- fread(penalty_segments.bed)
    penalty.loss <- fread(penalty_loss.tsv)
    setnames(penalty.loss, c(
      "penalty", "segments", "peaks", "bases",
      "mean.pen.cost", "total.cost", "status",
      "mean.intervals", "max.intervals"))
    if(penalty.loss$segments == nrow(penalty.segs)){
      TRUE
    }else{
      FALSE
    }
  }, error=function(e){
    FALSE
  })
  if(!already.computed){
    penalty.db <- paste0(pre, ".db")
    fpop.cmd <- paste(
      "PeakSegFPOP", prob.cov.bedGraph, penalty.str, penalty.db)
    cat(fpop.cmd, "\n")
    seconds <- system.time({
      system(fpop.cmd)
    })[["elapsed"]]
    megabytes <- if(file.exists(penalty.db)){
      file.size(penalty.db)/1024/1024
    }else{
      0
    }
    timing <- data.table(
      penalty.str,
      megabytes,
      seconds)
    write.table(
      timing,
      penalty_timing.tsv,
      row.names=FALSE, col.names=FALSE,
      quote=FALSE, sep="\t")
    unlink(penalty.db)
  }
  timing <- fread(penalty_timing.tsv)
  setnames(timing, c("penalty", "megabytes", "seconds"))
  penalty.loss <- fread(penalty_loss.tsv)
  setnames(penalty.loss, c(
    "penalty", "segments", "peaks", "bases",
    "mean.pen.cost", "total.cost", "status",
    "mean.intervals", "max.intervals"))
  penalty.segs <- fread(penalty_segments.bed)
  setnames(penalty.segs, c("chrom","chromStart", "chromEnd", "status", "mean"))
  list(
    segments=penalty.segs,
    loss=penalty.loss,
    timing=timing)
### List of data.tables: segments has one row for every segment in the
### optimal model, loss has one row and contains the Poisson loss and
### feasibility, and timing is one row with the time and disk usage.
}  

problem.features <- function
### Make sure features are computed for one segmentation problem.
(problem.dir
### problemID directory with problemID/coverage.bedGraph.
 ){
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  coverage <- fread(file.path(problem.dir, "coverage.bedGraph"))
  setnames(coverage, c("chrom", "chromStart", "chromEnd", "count"))
  bases <- with(coverage, chromEnd-chromStart)
  long <- rep(coverage$count, bases)
  diff.vec <- abs(diff(long))
  feature.vec <- c(
    quartile=quantile(long),
    mean=mean(long),
    sd=sd(long),
    bases=sum(bases),
    data=nrow(coverage))
  log.features <- suppressWarnings(c(
    feature.vec,
    `log+1`=log(feature.vec+1),
    log=log(feature.vec),
    log.log=log(log(feature.vec))))
  feature.dt <- data.table(t(log.features))
  write.table(
    feature.dt,
    file.path(problem.dir, "features.tsv"),
    quote=FALSE,
    row.names=FALSE,
    col.names=TRUE,
    sep="\t")
### Nothing, but writes problemID/features.tsv if it does not exist
### already.
}

problem.target <- function
### Compute target interval for a segmentation problem.
(problem.dir
### problemID directory in which coverage.bedGraph has already been
### computed. If there is a labels.bed file then the number of
### incorrect labels will be computed in order to find the target
### interval of minimal error penalty values.
 ){
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  ## Check if problem/labels.bed exists.
  problem.labels <- tryCatch({
    prob.lab.bed <- file.path(problem.dir, "labels.bed")
    problem.labels <- fread(prob.lab.bed)
    setnames(problem.labels, c("chrom", "chromStart", "chromEnd", "annotation"))
    problem.labels
  }, error=function(e){
    data.frame(
      chrom=character(),
      chromStart=integer(),
      chromEnd=integer(),
      annotation=character())
  })

  ## Compute the target interval. Even if there are no labels, we can
  ## still compute an upper bound on the number of peaks (lower limit of
  ## target interval of penalty values).
  error.list <- list()
  getError <- function(penalty.str){
    stopifnot(is.character(penalty.str))
    result <- problem.PeakSegFPOP(problem.dir, penalty.str)
    penalty.peaks <- result$segments[status=="peak",]
    penalty.error <- PeakErrorChrom(penalty.peaks, problem.labels)
    error.list[[penalty.str]] <<- with(penalty.error, data.table(
      result$loss,
      fn=sum(fn),
      fp=sum(fp)))
    do.call(rbind, error.list)[order(penalty),]
  }

  error.dt <- getError("0")
  error.dt <- getError("Inf")
  min.fn <- error.list[["0"]]$fn
  max.fp <- error.list[["0"]]$fp
  max.fn <- error.list[["Inf"]]$fn

  ## mx+b = lossInf => x = (lossInf-b)/m
  lossInf <- error.list[["Inf"]]$total.cost
  next.pen <- with(error.list[["0"]], (lossInf-total.cost)/peaks)
  next.side <- "upper"
  while(!is.null(next.pen)){
    if(interactive()){
      gg <- ggplot()+
        geom_abline(aes(slope=peaks, intercept=total.cost),
                    data=error.dt)+
        geom_vline(aes(xintercept=penalty),
                   data=data.table(penalty=next.pen))+
        geom_point(aes(penalty, mean.pen.cost*bases),
                   data=error.dt)
      print(gg)
    }
    next.str <- paste(next.pen)
    error.dt <- getError(next.str)
    error.dt[, errors := ifelse(status=="feasible", fp, Inf)+fn]
    print(error.dt[,.(penalty, peaks, status, fp, fn, errors)])
    peaks.tab <- table(error.dt$peaks)
    check <- function(is.vec){
      stopifnot(is.logical(is.vec))
      stopifnot(nrow(error.dt) == length(is.vec))
      i.vec <- which(is.vec)
      if(length(i.vec)==0)return(NULL)
      two.i <- if(1 %in% i.vec){
        m <- max(i.vec)
        c(m, m+1)
      }else{
        m <- min(i.vec)
        c(m-1, m)
      }
      if(any(!two.i %in% seq_along(error.dt$peaks)))return(NULL)
      two.error <- error.dt[two.i,]
      adjacent.models <- diff(two.error$peaks) == -1
      two.lambda <- any(1 < peaks.tab[paste(two.error$peaks)])
      data.table(
        penalty=(two.error[1, total.cost]-two.error[2, total.cost])/
          (two.error[2, peaks]-two.error[1, peaks]),
        found=two.lambda||adjacent.models)
    }
    feasible <- check(error.dt$status=="infeasible")
    fp <- check(error.dt$fp==0)
    fn <- check(error.dt$fn==min.fn)
    min.errors <- min(error.dt$errors)
    min.errors.i.vec <- which(error.dt$errors==min.errors)
    min.errors.first <- min(min.errors.i.vec)
    min.errors.last <- max(min.errors.i.vec)
    before.min.error <- check(seq_along(error.dt$errors) < min.errors.first)
    after.min.error <- check(min.errors.last < seq_along(error.dt$errors))
    lower.candidates <- rbind(feasible, fp, before.min.error)[order(penalty),]
    lower <- lower.candidates[.N,]
    upper <- if(!is.null(after.min.error)){
      after.min.error
    }else{
      fn
    }
    next.pen <- if(0 == max.fn){
      ## Do not search for the upper limit when there are no
      ## positive labels.
      if(!lower$found){
        lower$penalty
      }
    }else{
      if((!lower$found) && (!upper$found)){
        ## Neither the upper limit nor the lower limit have been
        ## found, so we can search for one or the other. One search
        ## strategy is to always search for the upper limit until it
        ## is found, and then search for the lower limit. That could
        ## be bad if we have an error profile as follows,
        ## H3K4me3_TDH_immune/McGill0322/problems/chr3:93504854-194041961

        ##  1:     0.0000 2495509 infeasible 23  0    Inf
        ##  2:    46.2596  297639 infeasible 23  0    Inf
        ##  3:   347.5916   65093 infeasible 20  0    Inf
        ##  4:  1132.9245    9402 infeasible 14  0    Inf
        ##  5:  4317.8561    1460 infeasible  4  2    Inf
        ##  6: 18045.5544     384 infeasible  0  3    Inf
        ##  7: 20763.5231     345 infeasible  0  3    Inf
        ##  8: 20877.9948     344 infeasible  0  3    Inf
        ##  9: 20904.2069     343   feasible  0  3      3
        ## 10: 20990.6296     342   feasible  0  3      3
        ## 11: 21165.0400     340   feasible  0  3      3
        ## 12: 21669.8758     337   feasible  0  3      3
        ## 13: 22530.0711     331   feasible  0  3      3
        ## 14: 23795.3361     314   feasible  0  3      3
        ## 15: 29561.2146     247   feasible  0  3      3
        ## 16: 30935.5267     230   feasible  0  2      2
        ## 17: 32698.7993     215   feasible  0  4      4
        ## 18: 36089.5234     181   feasible  0  4      4
        ## 19: 45600.3192     126   feasible  0  5      5
        ## 20:        Inf       0   feasible  0 11     11

        ## Using the old/bad search method, after finding the upper
        ## limit between 3 and Inf errors, the lower limit search
        ## discovered a model with only 2 errors, so had to re-start
        ## the upper limit search. Using the new search method below,
        ## we alternate between looking for the upper and lower
        ## limits.
        if(next.side == "upper"){
          next.side <- "lower"
          upper$penalty
        }else{
          next.side <- "upper"
          lower$penalty
        }
      }else if(!lower$found){
        lower$penalty
      }else if(!upper$found){
        upper$penalty
      }
    }
  }#while(!is.null(pen))

  write.table(
    error.dt,
    file.path(problem.dir, "target_models.tsv"),
    quote=FALSE,
    row.names=FALSE,
    col.names=TRUE)

  error.sorted <- error.dt[order(peaks), ][c(TRUE, diff(peaks) != 0),]

  ## Oracle model complexity?

  ## error.sorted[, `:=`(complexity={
  ##   in.sqrt <- 1.1 + log(bases / segments)
  ##   in.square <- 1 + 4 * sqrt(in.sqrt)
  ##   in.square * in.square * segments
  ## })]

  path <- error.sorted[, exactModelSelection(total.cost, peaks, peaks)]
  setkey(error.sorted, peaks)
  path$errors <- error.sorted[J(path$peaks), errors]
  indices <- with(path, largestContinuousMinimum(
    errors, max.log.lambda-min.log.lambda))
  target <- with(path, data.table(
    min.log.lambda=min.log.lambda[indices$start],
    max.log.lambda=max.log.lambda[indices$end]))
  write.table(
    target,
    file.path(problem.dir, "target.tsv"),
    quote=FALSE,
    col.names=FALSE,
    row.names=FALSE)
  list(
    target=as.numeric(unlist(target)),
    models=error.dt,
    selected=path)
### List of info related to target interval computation: target is the
### interval of log(penalty) values that achieve minimum incorrect
### labels (numeric vector of length 2), models is a data.table with
### one row per model for which the label error was computed, selected
### is a data.frame that describes which penalty values will select
### which models.
}
