pattern <- paste0(
  "=(?<fun>.*?)\n",
  "(?<table>",
  "(?:[^=].*?\n)*",
  ")")
library(namedCapture)
library(data.table)
library(ggplot2)
MAX <- 5
MIN <- -MAX
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
      min.mean <- max(min_mean, MIN)
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
                   0.00000000000000000000e+00 0.00000000000000000000e+00 1.71372993871380030484e+00            -inf             inf       -0.335433 1
=cost model
                   Square     Linear        Constant        min_mean        max_mean       prev_mean data_i
                   4.00000000000000000000e+00 5.70600000000000040501e+00 3.03490225000000046762e+00            -inf       -0.490978       -0.481957 0
                   2.00000000000000000000e+01 1.34173040000000014516e+01 2.96403052156900059799e+00       -0.490978        0.009022        0.000000 -1
                   4.00000000000000000000e+00 5.70600000000000040501e+00 3.03490225000000046762e+00        0.009022             inf       -0.481957 0
=result
    Square     Linear        Constant        min_mean        max_mean       prev_mean data_i
4.00000000000000000000e+00 5.70600000000000040501e+00 3.03490225000000046762e+00            -inf       -0.490978       -0.481957 0
                   2.00000000000000000000e+01 1.34173040000000014516e+01 2.96403052156900059799e+00       -0.490978       -0.111826        0.000000 -1
                   0.00000000000000000000e+00 0.00000000000000000000e+00 1.71372993871380030484e+00       -0.111826             inf       -0.335433 1
                  ")
ggplot()+
  # geom_vline(xintercept=xi, linetype="dashed")+
  geom_line(aes(mean, cost, color=fun),
            size=1,
            alpha=1,
            data=viz.list$funs) #+ geom_vline(xintercept=-1.176176) 


#+ facet_wrap(~fun) +
 #                 geom_vline(xintercept=0.4054054) +
  #                geom_hline(yintercept = 2.284542)
  
  
  #+ geom_vline(xintercept=y[2], linetype="dashed")

# + geom_vline(xintercept=y[2], linetype="dashed")+ geom_vline(xintercept=y[3], linetype="dashed")

