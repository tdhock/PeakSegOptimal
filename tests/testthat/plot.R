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

C1.60minless <- gdata("
=prev cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.00000000000000000000e+00 -3.04876615746180945621e+00 0.00000000000000000000e+00            -inf        3.555348        0.000000 -1
=min prev cost
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.00000000000000000000e+00 -3.04876615746180945621e+00 0.00000000000000000000e+00            -inf        1.114737             inf 59
0.00000000000000000000e+00 0.00000000000000000000e+00 -3.49806191860402737603e-01        1.114737        3.555348        1.114737 59
")
ggplot()+
  geom_vline(aes(xintercept=min_log_mean, color=fun),
             data=C1.60minless$vlines)+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=C1.60minless$funs)

C1.60minenv <- gdata("
=min prev cost
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.00000000000000000000e+00 -3.04876615746180945621e+00 0.00000000000000000000e+00            -inf        1.114737             inf 59
0.00000000000000000000e+00 0.00000000000000000000e+00 -3.49806191860402737603e-01        1.114737        3.555348        1.114737 59
=cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
9.99999999999999888978e-01 -3.04876615746180945621e+00 0.00000000000000000000e+00            -inf        1.189116             inf 56
7.16803760282021112937e-02 0.00000000000000000000e+00 -5.76570928152859929483e-01        1.189116        1.204649        1.189116 58
8.96004700352526439744e-01 -3.04876615746180945621e+00 0.00000000000000000000e+00        1.204649        1.493913            -inf 0
8.78965922444183145323e-01 -3.03172737955346693894e+00 5.04446708726733333839e-02        1.493913        1.562189       -1.960580 1
3.17861339600470027555e-01 -1.42361927144535838075e+00 2.14323509143853258019e-01        1.562189        1.602478        0.868121 36
3.05522914218566477018e-01 -1.38014101057579341436e+00 2.05914795383460108580e-01        1.602478        2.358202        0.876596 38
7.16803760282021112937e-02 0.00000000000000000000e+00 -5.76570928152859929483e-01        2.358202        3.555348        1.189116 58
=new cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
9.99999999999999888978e-01 -3.04876615746180945621e+00 0.00000000000000000000e+00            -inf        1.114737             inf 56
0.00000000000000000000e+00 0.00000000000000000000e+00 -3.49806191860402737603e-01        1.114737        1.204649        1.114737 59
8.96004700352526439744e-01 -3.04876615746180945621e+00 0.00000000000000000000e+00        1.204649        1.493913            -inf 0
8.78965922444183145323e-01 -3.03172737955346693894e+00 5.04446708726733333839e-02        1.493913        1.562189       -1.960580 1
3.17861339600470027555e-01 -1.42361927144535838075e+00 2.14323509143853258019e-01        1.562189        1.602478        0.868121 36
3.05522914218566477018e-01 -1.38014101057579341436e+00 2.05914795383460108580e-01        1.602478        1.934180        0.876596 38
0.00000000000000000000e+00 0.00000000000000000000e+00 -3.49806191860402737603e-01        1.934180        3.555348        1.114737 59
")
ggplot()+
  coord_cartesian(xlim=c(1, 2), ylim=c(-1, 0))+
  geom_vline(aes(xintercept=min_log_mean, color=fun),
             data=C1.60minenv$vlines)+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=C1.60minenv$funs)


C1.60minenv <- gdata("
=cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
9.99999999999999888978e-01 -3.04876615746180945621e+00 0.00000000000000000000e+00            -inf        1.189116             inf 56
7.16803760282021112937e-02 0.00000000000000000000e+00 -5.76570928152859929483e-01        1.189116        1.204649        1.189116 58
8.96004700352526439744e-01 -3.04876615746180945621e+00 0.00000000000000000000e+00        1.204649        1.493913            -inf 0
8.78965922444183145323e-01 -3.03172737955346693894e+00 5.04446708726733333839e-02        1.493913        1.562189       -1.960580 1
3.17861339600470027555e-01 -1.42361927144535838075e+00 2.14323509143853258019e-01        1.562189        1.602478        0.868121 36
3.05522914218566477018e-01 -1.38014101057579341436e+00 2.05914795383460108580e-01        1.602478        2.358202        0.876596 38
7.16803760282021112937e-02 0.00000000000000000000e+00 -5.76570928152859929483e-01        2.358202        3.555348        1.189116 58
")
ggplot()+
  coord_cartesian(xlim=c(1, 2), ylim=c(-1, 0))+
  geom_vline(aes(xintercept=min_log_mean, color=fun),
             data=C1.60minenv$vlines)+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=C1.60minenv$funs)

C1.1 <- gdata("
=prev cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.00000000000000000000e+00 -0.00000000000000000000e+00 0.00000000000000000000e+00            -inf        3.555348        0.000000 -1
=new cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
0.00000000000000000000e+00 0.00000000000000000000e+00 0.00000000000000000000e+00            -inf        3.555348            -inf 0
")
ggplot()+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=C1.1$funs)

C1.2prev <- gdata("
=prev cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.00000000000000000000e+00 -1.40776699029126206586e-01 0.00000000000000000000e+00            -inf        3.555348        0.000000 -1
=min prev cost
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.00000000000000000000e+00 -1.40776699029126206586e-01 0.00000000000000000000e+00            -inf       -1.960580             inf 1
0.00000000000000000000e+00 0.00000000000000000000e+00 4.16780727307233478385e-01       -1.960580        3.555348       -1.960580 1
")
ggplot()+
  geom_vline(aes(xintercept=min_log_mean, color=fun),
             data=C1.2prev$vlines)+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=C1.2prev$funs)

C1.54 <- gdata("
=min prev cost
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
6.78554521066750820912e-03 0.00000000000000000000e+00 1.17256803266681891573e-01            -inf      -44.368374             inf 209
1.11251380779548686428e-01 -6.31213507969070625886e-03 -1.62802368939145597482e-01      -44.368374       -3.025947             inf 209
1.55594129714375889462e-01 -3.77150071011519635866e-02 -2.59976955043114466015e-01       -3.025947       -1.417193             inf 209
0.00000000000000000000e+00 0.00000000000000000000e+00 -1.68812511325519581940e-01       -1.417193        3.555348       -1.417193 209
=cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.30976802903582136006e-02 -6.31213507969070625886e-03 -1.62802368939145597482e-01            -inf       -2.859501            -inf 207
1.41865235916048632037e-01 -3.77150071011519635866e-02 -2.59976955043114466015e-01       -2.859501       -0.992364            -inf 203
1.30976802903582136006e-02 -6.31213507969070625886e-03 -1.81079822644533677822e-01       -0.992364        0.975518       -1.512418 207
6.78554521066750820912e-03 0.00000000000000000000e+00 -1.70494226206477972330e-01        0.975518        3.555348       -1.372603 208
=new cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
1.30976802903582136006e-02 -6.31213507969070625886e-03 -1.62802368939145597482e-01            -inf       -2.859501            -inf 207
1.41865235916048632037e-01 -3.77150071011519635866e-02 -2.59976955043114466015e-01       -2.859501       -0.992364            -inf 203
1.30976802903582136006e-02 -6.31213507969070625886e-03 -1.81079822644533677822e-01       -0.992364       -0.140584       -1.512418 207
0.00000000000000000000e+00 0.00000000000000000000e+00 -1.68812511325519581940e-01       -0.140584        3.555348       -1.417193 209
")
ggplot()+
  geom_vline(aes(xintercept=min_log_mean, color=fun),
             data=C1.54$vlines)+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=C1.54$funs)
