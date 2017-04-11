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

C12.221minless <- gdata("
=min-less/more
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.00000000000000000000e+00 -4.84507042253521147401e-01 0.00000000000000000000e+00            -inf       -0.724623             inf 2
0.00000000000000000000e+00 0.00000000000000000000e+00 8.35592140219317158767e-01       -0.724623        3.496508       -0.724623 2
=cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.00000000000000000000e+00 -4.84507042253521147401e-01 0.00000000000000000000e+00            -inf        0.552375             inf 1
7.21126760563380320157e-01 0.00000000000000000000e+00 2.16877645665452228885e-01        0.552375        3.496508        0.552375 1
=new cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.00000000000000000000e+00 -4.84507042253521147401e-01 0.00000000000000000000e+00            -inf       -0.724623             inf 2
0.00000000000000000000e+00 0.00000000000000000000e+00 8.35592140219317158767e-01       -0.724623        3.496508       -0.724623 2
")
C12.221minless <- gdata("
=min-less/more
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.00000000000000000000e+00 -4.84507042253521147401e-01 0.00000000000000000000e+00            -inf             inf             inf 2
=cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.00000000000000000000e+00 -4.84507042253521147401e-01 0.00000000000000000000e+00            -inf             inf             inf 1
=new cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.00000000000000000000e+00 -4.84507042253521147401e-01 0.00000000000000000000e+00            -inf             inf             inf 2
")
xi <- 2.024689
gg <- ggplot()+
  geom_vline(xintercept=xi, linetype="dashed")+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            alpha=0.5,
            data=C12.221minless$funs)
if(!is.null(C12.221minless$vlines)){
  gg <- gg+
    geom_vline(aes(xintercept=min_log_mean, color=fun),
               data=C12.221minless$vlines)
}
print(gg)

gg+
  coord_cartesian(ylim=c(-1.6, -1.4), xlim=c(0, 2.5))

