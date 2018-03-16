pattern <- paste0(
  "=(?<fun>.*?)\n",
  "(?<table>",
  "(?:[^=].*?\n)*",
  ")")
library(namedCapture)
library(data.table)
library(ggplot2)
ploss <- function(dt, x){
  ## need to make a new data table, otherwise ifelse may only get one
  ## element, and return only one element.
  new.dt <- data.table(dt, x)
  new.dt[, ifelse(Log==0, 0, Log*log(x)) + Linear*x + Constant]
}
getLines <- function(dt){
  line.list <- list()
  print(dt)
  for(piece.i in 1:nrow(dt)){
    piece <- dt[piece.i,]
    mean.vec <- piece[, {
      min.mean <- exp(min_log_mean)
      max.mean <- if(max_log_mean==Inf){
        min.mean+1
      }else{
        exp(max_log_mean)
      }
      seq(min.mean, max.mean, l=1000)
    }]
    line.list[[piece.i]] <- data.table(
      piece.i,
      piece,
      mean=mean.vec,
      log.mean=log(mean.vec),
      cost=ploss(piece, mean.vec))
  }
  do.call(rbind, line.list)
}
getLines <- function(dt){
  line.list <- list()
  print(dt)
  for(piece.i in 1:nrow(dt)){
    piece <- dt[piece.i,]
    log.mean.vec <- piece[, {
      min.log.mean <- if(min_log_mean==-Inf){
        -600
      }else{
        min_log_mean
      }
      max.log.mean <- if(max_log_mean==Inf){
        min.log.mean+1
      }else{
        max_log_mean
      }
      seq(min.log.mean, max.log.mean, l=100)
    }]
    line.list[[piece.i]] <- data.table(
      piece.i,
      piece,
      log.mean=log.mean.vec,
      ##cost=piece[, ifelse(Log==0, 0, Log*log.mean.vec) + Linear*exp(log.mean.vec) + Constant])
      cost=ploss(piece, exp(log.mean.vec)))
  }
  do.call(rbind, line.list)
}
gdata <- function(txt){
  mat <- str_match_all_named(txt, pattern)[[1]]
  funs.list <- list()
  vlines.list <- list()
  coef.list <- list()
  for(row.i in 1:nrow(mat)){
    r <- mat[row.i,]
    df <- read.table(text=r[["table"]], header=TRUE)
    dt <- data.table(df)
    l <- getLines(dt)
    fun <- r[["fun"]]
    coef.list[[fun]] <- dt
    funs.list[[row.i]] <- data.table(fun, l)
    if(1 < nrow(dt)){
      vlines.list[[row.i]] <- data.table(fun, dt[-1,])
    }
  }
  list(
    funs=do.call(rbind, funs.list),
    vlines=do.call(rbind, vlines.list),
    coefs=coef.list)
}

viz.list <- gdata("
=prev cost
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
2.00000000000000011102e-01 -4.00000000000000000000e+00 -1.41163501263891912885e+01        0.693147        2.570013        2.803360 3
4.00000000000000022204e-01 -7.00000000000000000000e+00 -9.01950994083225410236e+00        2.570013        2.840016        2.890372 2
2.00000000000000011102e-01 -4.00000000000000000000e+00 -1.41163501263891912885e+01        2.840016        2.995732        2.803360 3
")
xi <- NA # the point at which the error was detected.
gg <- ggplot()+
  geom_vline(xintercept=xi, linetype="dashed")+
  geom_point(aes(log.mean, cost, color=fun),
            size=2,
            alpha=0.5,
            data=viz.list$funs)
if(!is.null(viz.list$vlines)){
  gg <- gg+
    geom_vline(aes(xintercept=min_log_mean, color=fun),
               data=viz.list$vlines)
}
print(gg)

sapply(viz.list$coefs, function(dt){
  ploss(dt[min_log_mean < xi & xi < max_log_mean], exp(xi))
})

gg+
  coord_cartesian(xlim=c(-290, -280), ylim=c(16225, 16300))

gg+
  coord_cartesian(ylim=c(-0.966165, -0.966160), xlim=c(3.51898, 3.519))

PeakSegPipeline::PeakSegFPOP_disk("~/coverage.bedGraph", "48402378.5676387")
