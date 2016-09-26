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
  for(piece.i in 1:nrow(dt)){
    piece <- dt[piece.i,]
    mean.vec <- piece[, seq(exp(min_log_mean), exp(max_log_mean), l=1000)]
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
  for(row.i in 1:nrow(mat)){
    r <- mat[row.i,]
    df <- read.table(text=r[["table"]], header=TRUE)
    dt <- data.table(df)
    l <- getLines(dt)
    fun <- r[["fun"]]
    funs.list[[row.i]] <- data.table(fun, l)
    if(1 < nrow(dt)){
      vlines.list[[row.i]] <- data.table(fun, dt[-1,])
    }
  }
  list(
    funs=do.call(rbind, funs.list),
    vlines=do.call(rbind, vlines.list))
}

C12.221minless <- gdata("
=prev cost
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.20481927710843379764e-02 -5.18072289156626508699e-01 -6.68495245036515512993e+01        2.639057        3.335874        3.761200 81
7.10843373493975971833e-01 -2.05421686746987894878e+01 -1.96898571312930741328e+01        3.335874        3.374668        3.601868 23
1.20481927710843379764e-02 -5.18072289156626508699e-01 -6.68495245036515512993e+01        3.374668        3.761200        3.761200 81
=min less/more(prev cost)
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.20481927710843379764e-02 -5.18072289156626508699e-01 -6.68495245036515512993e+01        2.639057        3.335874             inf 81
7.10843373493975971833e-01 -2.05421686746987894878e+01 -1.96898571312930741328e+01        3.335874        3.363783             inf 23
0.00000000000000000000e+00 0.00000000000000000000e+00 -6.82470851145374695079e+01        3.363783        3.382081        3.363783 23
1.20481927710843379764e-02 -5.18072289156626508699e-01 -6.68495245036515512993e+01        3.382081        3.761200             inf 81
0.00000000000000000000e+00 0.00000000000000000000e+00 inf        3.761200        3.761200        3.363783 23
")
ggplot()+
  geom_vline(aes(xintercept=min_log_mean, color=fun),
             data=C12.221minless$vlines)+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=C12.221minless$funs)

