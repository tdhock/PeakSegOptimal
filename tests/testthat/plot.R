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

C11.301minless <- gdata("
=prev cost model
    Linear        Log   Constant min_log_mean max_log_mean     data_i
        22          0 -19706.460691       -inf  -1.917963 299
        71        -49 -19807.639270  -1.917963   0.370000 298
       402       -797 -20010.079492   0.370000   0.445068 286
       421       -844 -20018.812620   0.445068   0.850674 284
       439       -898 -20015.018261   0.850674   0.865637 283
      1097      -3300 -19499.507878   0.865637   1.127393 238
      1098      -3309 -19492.448939   1.127393   1.414200 237
      1100      -3329 -19472.391333   1.414200   1.482339 236
        71        -49 -19803.536430   1.482339   3.496508 299
=min prev cost
    Linear        Log   Constant min_log_mean max_log_mean     data_i
         0          0 -19706.460691       -inf  -1.833168 300
        71        -49 -19807.639270  -1.833168  -0.370860 300
         0          0 -19740.467150  -0.370860   0.461754 300
       421       -844 -20018.812620   0.461754   0.695520 300
         0          0 -19761.831214   0.695520   0.884747 300
      1097      -3300 -19499.507878   0.884747   1.101343 300
         0          0 -19833.940726   1.101343   3.496508 300
")
ggplot()+
  geom_vline(xintercept=log(2.85))+
  geom_vline(xintercept=0.656780)+
  coord_cartesian(xlim=c(0,2), ylim=c(-20000, -19500))+
  geom_vline(aes(xintercept=min_log_mean, color=fun),
             data=C11.301minless$vlines)+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=C11.301minless$funs)

C11.301minenv <- gdata("
=min prev cost
    Linear        Log   Constant min_log_mean max_log_mean     data_i
         0          0 -19706.460691       -inf  -1.833168 300
        71        -49 -19807.639270  -1.833168  -0.370860 300
         0          0 -19740.467150  -0.370860   0.461754 300
       421       -844 -20018.812620   0.461754   0.695520 300
         0          0 -19761.831214   0.695520   0.884747 300
      1097      -3300 -19499.507878   0.884747   1.101343 300
         0          0 -19833.940726   1.101343   3.496508 300
=cost model
    Linear        Log   Constant min_log_mean max_log_mean     data_i
        71        -49 -19807.639270       -inf   0.000000 299
        22          0 -19758.639270   0.000000   0.309734 299
       213       -499 -19864.426823   0.309734   0.930885 293
      1097      -3300 -19499.507878   0.930885   1.121603 299
        22          0 -19900.793841   1.121603   3.496508 299
=new cost model
    Linear        Log   Constant min_log_mean max_log_mean     data_i
         0          0 -19706.460691       -inf  -1.833168 300
        71        -49 -19807.639270  -1.833168  -0.370860 300
         0          0 -19740.467150  -0.370860   0.368836 300
       213       -499 -19864.426823   0.368836   0.930885 293
      1097      -3300 -19499.507878   0.930885   1.101343 300
         0          0 -19833.940726   1.101343   3.496508 300
")
ggplot()+
  geom_vline(xintercept=log(2.85))+
  geom_vline(xintercept=0.656780)+
  coord_cartesian(xlim=c(0,2), ylim=c(-20000, -19500))+
  geom_vline(aes(xintercept=min_log_mean, color=fun),
             data=C11.301minenv$vlines)+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=C11.301minenv$funs)
## The min comes from (should not be here)
##        213       -499 -19864.426823   0.368836   0.930885 293

C11.301 <- gdata("
=new cost model
    Linear        Log   Constant min_log_mean max_log_mean     data_i
        95        -95 -19706.460691       -inf  -1.833168 300
       166       -144 -19807.639270  -1.833168  -0.370860 300
        95        -95 -19740.467150  -0.370860   0.368836 300
       308       -594 -19864.426823   0.368836   0.930885 293
      1192      -3395 -19499.507878   0.930885   1.101343 300
        95        -95 -19833.940726   1.101343   3.496508 300
")
ggplot()+
  geom_vline(xintercept=log(2.85))+
  coord_cartesian(xlim=c(-1,2), ylim=c(-19700, -19600))+
  geom_vline(aes(xintercept=min_log_mean, color=fun),
             data=C11.301$vlines)+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=C11.301$funs)+
  geom_point(aes(x=0.656780, y=-19660.553868))

## prev cost model is C10.293.
C11.294minless <- gdata("
=prev cost model
    Linear        Log   Constant min_log_mean max_log_mean     data_i
        40          0 -19864.426823       -inf  -2.005918 292
        57        -17 -19900.814557  -2.005918   0.042381 291
       105       -113 -19946.823999   0.042381   0.228953 290
       189       -298 -20010.079492   0.228953   0.445068 286
       208       -345 -20018.812620   0.445068   0.850674 284
       226       -399 -20015.018261   0.850674   0.865637 283
       884      -2801 -19499.507878   0.865637   1.127393 238
       885      -2810 -19492.448939   1.127393   1.414200 237
       887      -2830 -19472.391333   1.414200   1.568696 236
        57        -17 -19900.814557   1.568696   3.496508 292
=min prev cost
    Linear        Log   Constant min_log_mean max_log_mean     data_i
         0          0 -19864.426823       -inf   0.930885 293
       884      -2801 -19499.507878   0.930885   1.127393 293
       885      -2810 -19492.448939   1.127393   1.155352 293
         0          0 -19928.988388   1.155352   3.496508 293
")
ggplot()+
  geom_vline(xintercept=0.656780)+
  geom_vline(xintercept=log(2.85))+
  coord_cartesian(ylim=c(-20000, -19800))+
  geom_vline(aes(xintercept=min_log_mean, color=fun),
             data=C11.294minless$vlines)+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=C11.294minless$funs)
## Is there a problem with creating a constant segment on the left?

C153.370minless <- gdata("
=prev cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
         5         -9   -20936.673374            -inf        0.000000             inf 368
        54        -58   -20985.673374        0.000000        0.072759             inf 368
         1         -1   -20932.820658        0.072759        3.496508        0.072759 368
=min prev cost
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
         0          0   -20931.817994            -inf        0.071459        0.071459 369
        54        -58   -20985.673374        0.071459        0.072759             inf 369
         1         -1   -20932.820658        0.072759        3.496508             inf 369
")
ggplot()+
  geom_vline(aes(xintercept=min_log_mean, color=fun),
             data=C153.370minless$vlines)+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=C153.370minless$funs)

C153.370minenv <- gdata("
=min prev cost
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
         0          0   -20931.817994            -inf        0.071459        0.071459 369
        54        -58   -20985.673374        0.071459        0.072759             inf 369
         1         -1   -20932.820658        0.072759        3.496508             inf 369
=cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
         1         -1   -20932.820658            -inf        0.072759        0.072759 368
        54        -58   -20985.673374        0.072759        0.438176             inf 368
         1         -1   -20928.505894        0.438176        0.693147        0.693147 368
         5         -9   -20930.960717        0.693147        3.496508             inf 368
=new cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
         0          0   -20931.817994            -inf       -0.073882        0.071459 369
         1         -1   -20932.820658       -0.073882        3.496508        0.072759 368
")
ggplot()+
  ##coord_cartesian(xlim=c(-1, 1), ylim=c(-20900, -20897.5))+
  geom_vline(aes(xintercept=min_log_mean, color=fun),
             data=C153.370minenv$vlines)+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=C153.370minenv$funs)


C153.370 <- gdata("
=prev down cost
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
       256          0     6834.434385            -inf      -39.735084        0.552375 1
       355       -172        0.000000      -39.735084        3.496508        0.000000 -1
=min less(prev down cost)
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
         0          0     6834.434385            -inf        3.496508            -inf 1
")
ggplot()+
  ##coord_cartesian(xlim=c(-1, 1), ylim=c(-20900, -20897.5))+
  geom_vline(aes(xintercept=min_log_mean, color=fun),
             data=C153.370$vlines)+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            data=C153.370$funs)

ggplot()+
  ##coord_cartesian(xlim=c(-1, 1), ylim=c(-20900, -20897.5))+
  geom_vline(aes(xintercept=exp(min_log_mean), color=fun),
             data=C153.370$vlines)+
  geom_line(aes(mean, cost, color=fun),
            size=2,
            data=C153.370$funs)

