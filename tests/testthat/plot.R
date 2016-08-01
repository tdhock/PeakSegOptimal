subject <- "
=prev up cost
    Linear        Log   Constant min_log_mean max_log_mean     data_i
         8         -6 229.114503       -inf   0.000000 7
         2          0 235.114503   0.000000   0.595818 7
        66       -122 191.675361   0.595818   0.693661 5
         2          0 235.114503   0.693661   3.496508 7
=min more(prev up cost)
    Linear        Log   Constant min_log_mean max_log_mean     data_i
         0          0 238.722672       -inf   0.614366 5
        66       -122 191.675361   0.614366   0.693661 5
         2          0 235.114503   0.693661   3.496508 7
"
pattern <- paste0(
  "=(?<fun>.*?)\n",
  "(?<table>",
  "(?:[^=].*?\n)*",
  ")")
library(namedCapture)
(mat <- str_match_all_named(subject, pattern)[[1]])
funs.list <- list()
library(data.table)
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
    mean.vec <- piece[, seq(exp(min_log_mean), exp(max_log_mean), l=100)]
    line.list[[piece.i]] <- data.table(
      piece.i,
      piece,
      mean=mean.vec,
      log.mean=log(mean.vec),
      cost=ploss(piece, mean.vec))
  }
  do.call(rbind, line.list)
}
vlines.list <- list()
for(row.i in 1:nrow(mat)){
  r <- mat[row.i,]
  df <- read.table(text=r[["table"]], header=TRUE)
  dt <- data.table(df)
  l <- getLines(dt)
  fun <- r[["fun"]]
  funs.list[[row.i]] <- data.table(fun, l)
  vlines.list[[row.i]] <- data.table(fun, dt[-1,])
}
funs <- do.call(rbind, funs.list)
vlines <- do.call(rbind, vlines.list)
library(ggplot2)
ggplot()+
  geom_vline(aes(xintercept=min_log_mean, color=fun), data=vlines)+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=funs)

ggplot()+
  geom_vline(aes(xintercept=min_log_mean, color=fun), data=vlines)+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=funs)+
  coord_cartesian(xlim=c(0.5,0.75), ylim=c(238, 240))
