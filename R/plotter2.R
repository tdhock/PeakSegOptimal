pattern <- paste0(
  "=(?<fun>.*?)\n",
  "(?<table>",
  "(?:[^=].*?\n)*",
  ")")
library(namedCapture)
library(data.table)
library(ggplot2)
MAX <- 10
sloss <- function(dt, x){
  ## need to make a new data table, otherwise ifelse may only get one
  ## element, and return only one element.
  new.dt <- data.table(dt, x)
  new.dt[, Square * x ^ 2 + Linear * x + Constant]
}
getLines <- function(dt){
  line.list <- list()
  print(dt)
  for(piece.i in 1:nrow(dt)){
    piece <- dt[piece.i,]
    mean.vec <- piece[, {
      min.mean <- min_mean
      max.mean <- min(max_mean, MAX)
      seq(min.mean, max.mean, l=1000)
    }]
    line.list[[piece.i]] <- data.table(
      piece.i,
      piece,
      mean=mean.vec,
      cost=sloss(piece, mean.vec))
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

viz.list <- gdata( "
=min-less/more
    Square     Linear        Constant        min_mean        max_mean       prev_mean data_i
                   0.00000000000000000000e+00 0.00000000000000000000e+00 2.00000000000000000000e+00        0.000000             inf        2.000000 1
=cost model
    Square     Linear        Constant        min_mean        max_mean       prev_mean data_i
                   1.00000000000000000000e+00 -4.00000000000000000000e+00 5.00000000000000000000e+00        0.000000        3.000000        4.000000 0
                   2.00000000000000000000e+00 -1.20000000000000000000e+01 2.00000000000000000000e+01        3.000000        5.000000        0.000000 -1
                   1.00000000000000000000e+00 -4.00000000000000000000e+00 5.00000000000000000000e+00        5.000000             inf        4.000000 0
=outerer
    Square     Linear        Constant        min_mean        max_mean       prev_mean data_i
0.00000000000000000000e+00 0.00000000000000000000e+00 2.00000000000000000000e+00        0.000000        1.000000        2.000000 1
                   1.00000000000000000000e+00 -4.00000000000000000000e+00 5.00000000000000000000e+00        1.000000        3.000000        4.000000 0
                   0.00000000000000000000e+00 0.00000000000000000000e+00 2.00000000000000000000e+00        3.000000        5.000000        2.000000 1
                  ")
ggplot()+
  # geom_vline(xintercept=xi, linetype="dashed")+
  geom_line(aes(mean, cost, color=fun),
            size=2,
            alpha=0.5,
            data=viz.list$funs) + 
  geom_vline(xintercept=3, linetype="dashed")+ geom_vline(xintercept=5, linetype="dashed")

