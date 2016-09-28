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
=prev up cost
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
3.26541274817136850770e-05 -1.67935512763098959569e-04 -6.33013001572306865938e+00            -inf        1.791759             inf 8322
2.79892521271831576694e-05 -1.39946260635915795123e-04 -6.33015217647847716620e+00        1.791759        2.016512        1.791759 8322
1.91259889535751610159e-04 -1.43678160919540249575e-03 -6.32876359403382249269e+00        2.016512        2.016512             inf 8321
1.91259889535751610159e-04 -1.43678160919540249575e-03 -6.32876359403382249269e+00        2.016512        2.072273             inf 8322
2.79892521271831576694e-05 -1.39946260635915795123e-04 -6.33015415563122463283e+00        2.072273        5.416100        2.072273 8322
=min more(prev up cost)
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
0.00000000000000000000e+00 0.00000000000000000000e+00 -6.33023709288205882473e+00            -inf        1.637609        1.637609 8322
3.26541274817136850770e-05 -1.67935512763098959569e-04 -6.33013001572306865938e+00        1.637609        1.791759             inf 8322
0.00000000000000000000e+00 0.00000000000000000000e+00 -6.33022412236162157484e+00        1.791759        2.016512        2.016512 8321
1.91259889535751610159e-04 -1.43678160919540249575e-03 -6.32876359403382249269e+00        2.016512        2.072273             inf 8322
2.79892521271831576694e-05 -1.39946260635915795123e-04 -6.33015415563122463283e+00        2.072273        5.416100             inf 8322
")
ggplot()+
  geom_vline(aes(xintercept=min_log_mean, color=fun),
             data=C12.221minless$vlines)+
  geom_vline(xintercept=1.904136, linetype="dashed")+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=C12.221minless$funs)

ggplot()+
  coord_cartesian(xlim=c(1.5, 2.5), ylim=c(-6.33025, -6.33))+
  geom_vline(aes(xintercept=min_log_mean, color=fun),
             data=C12.221minless$vlines)+
  geom_vline(xintercept=1.904136, linetype="dashed")+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=C12.221minless$funs)

