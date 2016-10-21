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
  ## If problemID/coverage.bedGraph has already been computed, than we
  ## have nothing to do.
  if(coverage.ok){
    return(NULL)
  }
  ## Create problemID/coverage.bedGraph from either
  ## sampleID/coverage.bigWig or sampleID/coverage.bedGraph.
  coverage.bedGraph <- file.path(sample.dir, "coverage.bedGraph")
  coverage.bigWig <- file.path(sample.dir, "coverage.bigWig")
  cov.cmd <- if(file.exists(coverage.bigWig)){
    problem[, sprintf(
      "bigWigToBedGraph -chrom=%s -start=%d -end=%d %s %s",
      chrom, problemStart, problemEnd,
      coverage.bigWig, prob.cov.bedGraph)]
  }else if(file.exists(coverage.bedGraph)){
    ## Use intersectBed -sorted to avoid memory
    ## problems. coverage.bedGraph needs to be -a since that is
    ## reported in the output.
    paste(
      "intersectBed -sorted",
      "-a", coverage.bedGraph,
      "-b", problem.bed,
      ">", prob.cov.bedGraph)
  }else{
    stop("To compute ", prob.cov.bedGraph,
         " need either ", coverage.bigWig,
         " or ", coverage.bedGraph,
         " which do not exist.")
  }
  cat(cov.cmd, "\n")
  status <- system(cov.cmd)
  if(status != 0){
    stop("non-zero status code ", status)
  }
### Nothing. If necessary, the intersectBed or bigWigToBedGraph
### command line program is used to create
### problemID/coverage.bedGraph, but it is not read into memory.
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
    penalty.segs <- fread(penalty_segments.bed)
  }
  timing <- fread(penalty_timing.tsv)
  setnames(timing, c("penalty", "megabytes", "seconds"))
  penalty.loss <- fread(penalty_loss.tsv)
  setnames(penalty.loss, c(
    "penalty", "segments", "peaks", "bases",
    "mean.pen.cost", "total.cost", "status",
    "mean.intervals", "max.intervals"))
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
### Compute features for one segmentation problem.
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

  ## Compute the label error for one penalty parameter.
  getError <- function(penalty.str){
    stopifnot(is.character(penalty.str))
    stopifnot(length(penalty.str) == 1)
    result <- problem.PeakSegFPOP(problem.dir, penalty.str)
    penalty.peaks <- result$segments[status=="peak",]
    penalty.error <- PeakErrorChrom(penalty.peaks, problem.labels)
    with(penalty.error, data.table(
      result$loss,
      fn=sum(fn),
      fp=sum(fp)))
  }
  
  ## Compute the target interval given the errors computed in dt.
  getTarget <- function(dt){
    peaks.tab <- table(dt$peaks)
    error.sorted <- dt[order(peaks), ][c(TRUE, diff(peaks) != 0),]
    error.sorted[, n.infeasible := cumsum(status=="infeasible")]
    error.sorted[, errors := fp + fn]
    setkey(error.sorted, peaks)
    path <- error.sorted[, exactModelSelection(total.cost, peaks, peaks)]
    path.dt <- data.table(path)
    setkey(path.dt, peaks)
    join.dt <- error.sorted[path.dt][order(penalty),]
    direction.list <- list(start=1, end=-1)
    side.vec.list <- list(fn="end", fp="start", errors=c("start", "end"))
    result <- list(models=path, candidates=list())
    for(error.col in c("fp", "fn", "errors")){
      incorrect.or.Inf <- ifelse(
        join.dt$n.infeasible==0, join.dt[[error.col]], Inf)
      indices <- join.dt[, largestContinuousMinimum(
        incorrect.or.Inf, max.log.lambda-min.log.lambda)]
      side.vec <- side.vec.list[[error.col]]
      for(side in side.vec){
        direction <- direction.list[[side]]
        index <- indices[[side]]
        model <- join.dt[index,]
        index.outside <- index - direction
        neighboring.peaks <- model$peaks + direction
        found.neighbor <- neighboring.peaks %in% join.dt$peaks
        multiple.penalties <- if(index.outside %in% seq_along(join.dt$peaks)){
          model.outside <- join.dt[index.outside,]
          peaks.num <- c(model.outside$peaks, model$peaks)
          peaks.str <- paste(peaks.num)
          peaks.counts <- peaks.tab[peaks.str]
          any(1 < peaks.counts)
        }else{
          FALSE
        }
        next.pen <- ifelse(side=="start", model$min.lambda, model$max.lambda)
        already.computed <- paste(next.pen) %in% names(error.list)
        done <- found.neighbor | multiple.penalties | already.computed
        result$candidates[[paste(error.col, side)]] <- data.table(
          model, found.neighbor, multiple.penalties, already.computed,
          done, next.pen)
      }
    }
    result
  }

  error.list <- list()
  next.pen <- c(0, Inf)
  while(length(next.pen)){
    cat("Next =", paste(next.pen, collapse=", "), "\n")
    next.str <- paste(next.pen)
    error.list[next.str] <- mclapply.or.stop(next.str, getError)
    error.dt <- do.call(rbind, error.list)[order(-penalty),]
    print(error.dt[,.(penalty, peaks, status, fp, fn)])
    target.list <- getTarget(error.dt)
    target.vec <- c(
      target.list$candidates[["errors start"]]$min.log.lambda,
      target.list$candidates[["errors end"]]$max.log.lambda)
    is.error <- grepl("error", names(target.list$candidates))
    error.candidates <- do.call(rbind, target.list$candidates[is.error])
    other.candidates <- do.call(rbind, target.list$candidates[!is.error])
    other.in.target <- other.candidates[done==FALSE &
        target.vec[1] < log(next.pen) & log(next.pen) < target.vec[2],]
    next.pen <- if(nrow(other.in.target)){
      other.in.target[, unique(next.pen)]
    }else{
      error.candidates[done==FALSE, unique(next.pen)]
    }
    if(interactive() && length(next.pen)){
      gg <- ggplot()+
        geom_abline(aes(slope=peaks, intercept=total.cost),
                    data=error.dt)+
        geom_vline(aes(xintercept=penalty),
                   color="red",
                   data=data.table(penalty=next.pen))+
        geom_point(aes(penalty, mean.pen.cost*bases),
                   data=error.dt)
      print(gg)
    }
  }#while(!is.null(pen))

  write.table(
    error.dt,
    file.path(problem.dir, "target_models.tsv"),
    sep="\t",
    quote=FALSE,
    row.names=FALSE,
    col.names=TRUE)

  write(target.vec, file.path(problem.dir, "target.tsv"), sep="\t")
  
  list(
    target=target.vec,
    models=error.dt)
### List of info related to target interval computation: target is the
### interval of log(penalty) values that achieve minimum incorrect
### labels (numeric vector of length 2), models is a data.table with
### one row per model for which the label error was computed, selected
### is a data.frame that describes which penalty values will select
### which models.
}

problem.predict <- function
### Predict peaks for a genomic segmentation problem.
(problem.dir,
### Problem directory.
 model.RData
### Model file created via train_model.R
 ){
  stopifnot(is.character(problem.dir))
  stopifnot(length(problem.dir)==1)
  stopifnot(is.character(model.RData))
  stopifnot(length(model.RData)==1)
  
  load(model.RData)

  problem.coverage(problem.dir)

  features.tsv <- file.path(problem.dir, "features.tsv")
  is.computed <- if(file.exists(features.tsv)){
    TRUE
  }else{
    tryCatch({
      problem.features(problem.dir)
      cat(sprintf("Computed %s\n", features.tsv))
      TRUE
    }, error=function(e){
      FALSE
    })
  }

  if(!is.computed){
    cat("Unable to compute", features.tsv, "so not predicting.\n")
    return(NULL)
  }

  features <- fread(features.tsv)
  feature.mat <- as.matrix(features)
  pred.penalty <- as.numeric(exp(model$predict(feature.mat)))
  n.features <- length(model$pred.feature.names)
  cat(paste0(
    "Predicting penalty=", pred.penalty,
    " log(penalty)=", log(pred.penalty),
    " based on ", n.features,
    " feature", ifelse(n.features==1, "", "s"),
    ".\n"))

  loss.glob <- file.path(problem.dir, "*_loss.tsv")
  loss.file.vec <- Sys.glob(loss.glob)
  loss.ord <- if(length(loss.file.vec)){
    loss <- fread(paste("cat", loss.glob))
    setnames(loss, c("penalty", "segments", "peaks", "bases", "mean.pen.cost", "total.cost", "status", "mean.intervals", "max.intervals"))
    loss[, log.penalty := log(penalty)]
    loss[order(-penalty),]
  }else{
    data.table()
  }

  ## TODO: If we have already computed the target interval and the
  ## prediction is outside, then we should choose the minimal error
  ## model which is closest to the predicted penalty.
  target.vec <- tryCatch({
    suppressWarnings(scan(file.path(problem.dir, "target.tsv"), quiet=TRUE))
  }, error=function(e){
    NULL
  })
  if(nrow(loss.ord) && length(target.vec)==2){
    cat(sprintf(
      "Target interval %f < log(penalty) < %f\n",
      target.vec[1], target.vec[2]))
    if(log(pred.penalty) < target.vec[1]){
      pred.penalty <- loss.ord[target.vec[1] < log.penalty, penalty[.N]]
      cat(sprintf(
        "Closest inside target penalty=%f log(penalty)=%f\n",
        pred.penalty, log(pred.penalty)))
    }
    if(target.vec[2] < log(pred.penalty)){
      pred.penalty <- loss.ord[log.penalty < target.vec[2], penalty[1]]
      cat(sprintf(
        "Closest inside target penalty=%f log(penalty)=%f\n",
        pred.penalty, log(pred.penalty)))
    }
  }  

  ## This will be NULL until we find or compute a model that can be used
  ## for predicted peaks.
  pen.str <- NULL

  ## If two neighboring penalties have already been computed, then we do
  ## not have to re-run PeakSegFPOP.
  if(is.null(pen.str) && 2 <= nrow(loss.ord)){
    is.after <- loss.ord[, penalty < pred.penalty]
    first.after <- which(is.after)[1]
    last.before <- first.after - 1
    smaller.peaks <- loss.ord[last.before, peaks]
    bigger.peaks <- loss.ord[first.after, peaks]
    if(any(is.after) && 0 < last.before && bigger.peaks - smaller.peaks <= 1){
      loss.unique <- loss.ord[c(TRUE, diff(peaks) != 0), ]
      exact <- loss.unique[, exactModelSelection(total.cost, peaks, peaks)]
      selected <- subset(
        exact, min.lambda < pred.penalty & pred.penalty < max.lambda)
      same.peaks <- loss.ord[peaks==selected$peaks, ]
      pen.num <- same.peaks$penalty[1]
      cat(
        "Based on previous computations, penalty of ",
        pred.penalty, " and ", 
        pen.num, " both recover ",
        selected$peaks, " peak",
        ifelse(selected$peaks==1, "", "s"), ".\n",
        sep="")
      pen.str <- paste(pen.num)
    }
  }

  ## TODO: If we have not already computed the target interval, then we
  ## can run PeakSegFPOP at the predicted penalty value. If the
  ## resulting model is feasible then we are done. Otherwise, we need to
  ## compute the target interval to find the biggest feasible model,
  ## which we return.
  if(is.null(pen.str)){
    pen.str <- paste(pred.penalty)
    result <- problem.PeakSegFPOP(problem.dir, pen.str)
    if(result$loss$status=="infeasible"){
      t.info <- problem.target(problem.dir)
      models.in.target <- with(t.info, {
        models[target[1] <= log(penalty) & log(penalty) <= target[2],]
      })
      biggest.feasible <- models.in.target[which.max(peaks),]
      pen.str <- paste(biggest.feasible$penalty)
    }
  }

  ## compute peaks.
  prob.cov.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
  pre <- paste0(prob.cov.bedGraph, "_penalty=", pen.str)
  penalty_segments.bed <- paste0(pre, "_segments.bed")
  penalty.segs <- fread(penalty_segments.bed)
  setnames(penalty.segs, c("chrom","chromStart", "chromEnd", "status", "mean"))
  peaks <- penalty.segs[status=="peak", ]
  peaks.bed <- file.path(problem.dir, "peaks.bed")
  cat(
    "Writing ", peaks.bed,
    " with ", nrow(peaks),
    " peak", ifelse(nrow(peaks)==1, "", "s"),
    " based on ", penalty_segments.bed,
    ".\n", sep="")

  write.table(
    peaks,
    peaks.bed,
    quote=FALSE,
    sep="\t",
    col.names=FALSE,
    row.names=FALSE)

  peaks
### data.table of peak predictions.
}
