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

viz.list <- gdata("
=min more(prev up cost)
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
0.00000000000000000000e+00 0.00000000000000000000e+00 -1.50579293620309506707e+00            -inf        1.098612        1.098612 13374
1.47372169993798088653e-04 -4.42116509981394238855e-04 -1.50574933808218780484e+00        1.098612        4.127134             inf 13374
=prev down cost
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.47372169993798088653e-04 -4.42116509981394238855e-04 -1.50582645029657946623e+00            -inf        0.000000        0.000000 13373
2.94744339987596177307e-04 -5.89488679975192354614e-04 -1.50597382246657351956e+00        0.000000        4.127134             inf 13373
=new down cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
0.00000000000000000000e+00 0.00000000000000000000e+00 -1.50579293620309506707e+00            -inf        0.000000        1.098612 13374
2.94744339987596177307e-04 -5.89488679975192354614e-04 -1.50597382246657351956e+00        0.000000        1.098612             inf 13373
1.47372169993798088653e-04 -4.42116509981394238855e-04 -1.50574933808218780484e+00        1.098612        4.127134             inf 13374
")
xi <- 0.549306 # the point at which the error was detected.
=======
=min less(prev down cost) + 866939314852865280.000000
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.00000000000000000000e+00 -1.95652173913043481157e+00 1.88465068446275080000e+16            -inf        0.671168             inf 1
0.00000000000000000000e+00 0.00000000000000000000e+00 1.88465068446275080000e+16        0.671168        0.693147        0.671168 1
=prev up cost
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.00000000000000000000e+00 -1.95652173913043481157e+00 1.88465068446275040000e+16            -inf        0.693147             inf 0
=new up cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.00000000000000000000e+00 -1.95652173913043481157e+00 1.88465068446275040000e+16            -inf        0.693147             inf 0
")
xi <- -0.306853 # the point at which the error was detected.
sapply(viz.list$coefs, function(dt){
  ploss(dt[min_log_mean < xi & xi < max_log_mean], exp(xi))
})

gg+
  coord_cartesian(xlim=c(-1, 0))

gg+
  coord_cartesian(ylim=c(-1.506, -1.505), xlim=c(-1, 2))

