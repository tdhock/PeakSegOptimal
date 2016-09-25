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
=prev up cost
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.64473684210526306659e-02 -1.15131578947368418131e-01 -2.51468095932980206797e+00            -inf        1.945910        1.945910 137
9.26535087719298489084e-02 -6.48574561403508886848e-01 -2.01009182828220334116e+00        1.945910        2.257400             inf 137
8.93640350877192984891e-02 -6.32127192982456231896e-01 -2.01577865491268060083e+00        2.257400        2.428837             inf 137
1.80921052631578954673e-02 -1.28289473684210508786e-01 -2.43088980736734638910e+00        2.428837        4.595120             inf 137
")
ggplot()+
  geom_vline(aes(xintercept=min_log_mean, color=fun),
             data=C12.221minless$vlines)+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=C12.221minless$funs)

