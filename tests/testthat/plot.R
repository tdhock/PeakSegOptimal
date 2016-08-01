subject <- "
=prev cost model
    Linear        Log   Constant min_log_mean max_log_mean     data_i
        40          0 -20724.638304       -inf  -1.867945 292
        57        -17 -20759.018862  -1.867945   0.100106 291
       105       -113 -20802.462502   0.100106   0.249963 290
       189       -298 -20864.073454   0.249963   0.456140 286
       208       -345 -20872.616321   0.456140   0.869969 284
       226       -399 -20868.601058   0.869969   0.926387 283
       667      -1793 -20690.905073   0.926387   1.136442 257
       887      -2830 -20197.860455   1.136442   1.397192 236
       888      -2841 -20186.535169   1.397192   1.443866 235
        57        -17 -20743.028552   1.443866   3.496508 292
=min prev cost
    Linear        Log   Constant min_log_mean max_log_mean     data_i
         0          0 -20724.638304       -inf   3.496508 293
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
  if(1 < nrow(dt)){
    vlines.list[[row.i]] <- data.table(fun, dt[-1,])
  }
}
funs <- do.call(rbind, funs.list)
vlines <- do.call(rbind, vlines.list)
library(ggplot2)
ggplot()+
  geom_vline(aes(xintercept=exp(min_log_mean), color=fun), data=vlines)+
  geom_line(aes(mean, cost, color=fun),
            size=2,
            data=funs)

ggplot()+
  geom_vline(aes(xintercept=exp(min_log_mean), color=fun), data=vlines)+
  coord_cartesian(ylim=c(-21000, -20500), xlim=c(0, 3))+
  geom_line(aes(mean, cost, color=fun),
            size=2,
            data=funs)
