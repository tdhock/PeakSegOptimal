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
0.00000000000000000000e+00 0.00000000000000000000e+00 6.31252201333190328114e+02            -inf        1.948792        1.948792 43
2.34819734345351055493e-02 -1.64848197343453517494e-01 6.31408607983290949051e+02        1.948792        2.488522             inf 43
1.80265654648956372141e-02 -1.43026565464895644153e-01 6.31420006356658745972e+02        2.488522        2.794012             inf 43
1.77893738140417444205e-02 -1.41840607210626190593e-01 6.31420570022667789090e+02        2.794012        3.399876             inf 43
4.74383301707779895303e-04 -3.32068311195445921291e-03 6.31468383313439289850e+02        3.399876        4.962845             inf 43
=cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
9.99999999999999555911e-01 -8.47248576850095136415e-01 0.00000000000000000000e+00            -inf        4.962845        0.000000 -1
")
ggplot()+
  geom_vline(aes(xintercept=min_log_mean, color=fun),
             data=C12.221minless$vlines)+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=C12.221minless$funs)

