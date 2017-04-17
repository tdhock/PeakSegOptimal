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

C12.221minless <- gdata("
=prev up cost
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
9.99999999999999666933e-01 -8.48470476642162907410e-01 6.30332141714930116905e+02            -inf      -25.095106             inf 43
9.84586198719468397300e-01 -8.48470476642162907410e-01 6.30332141714930230592e+02      -25.095106       -1.635367            -inf 0
6.62319184254209103457e-01 -8.01517666587621580021e-01 6.30471730450081622621e+02       -1.635367        0.560019       -1.972958 4
6.34100071140621057708e-01 -7.82783969646668409403e-01 6.30510642587088227629e+02        0.560019        0.611887       -1.717466 6
2.37135404315864384284e-02 -1.66231918425420926999e-01 6.31258878647748360891e+02        0.611887        2.488522       -0.358377 39
1.82594261323215552306e-02 -1.44415461228361408086e-01 6.31270274318158271853e+02        2.488522        2.564949       -0.332471 40
=min more(prev up cost)
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
0.00000000000000000000e+00 0.00000000000000000000e+00 6.31101400884306713124e+02            -inf        1.947338        1.947338 44
2.37135404315864384284e-02 -1.66231918425420926999e-01 6.31258878647748360891e+02        1.947338        2.488522             inf 44
1.82594261323215552306e-02 -1.44415461228361408086e-01 6.31270274318158271853e+02        2.488522        2.564949             inf 44
=prev down cost
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
2.37135404315864370189e-04 -1.42281242589518616692e-03 6.31102509087201951843e+02            -inf     -745.061389        1.948792 43
9.99999999999999666933e-01 -8.48470476642162907410e-01 0.00000000000000000000e+00     -745.061389        2.564949        0.000000 -1
=new down cost model
    Linear        Log        Constant    min_log_mean    max_log_mean   prev_log_mean data_i
0.00000000000000000000e+00 0.00000000000000000000e+00 6.31101400884306713124e+02            -inf     -745.061389        1.947338 44
9.99999999999999666933e-01 -8.48470476642162907410e-01 0.00000000000000000000e+00     -745.061389        2.564949        0.000000 -1
")
xi <- -745.061389
gg <- ggplot()+
  geom_vline(xintercept=xi, linetype="dashed")+
  geom_line(aes(log.mean, cost, color=fun),
            size=2,
            alpha=0.5,
            data=C12.221minless$funs)
if(!is.null(C12.221minless$vlines)){
  gg <- gg+
    geom_vline(aes(xintercept=min_log_mean, color=fun),
               data=C12.221minless$vlines)
}
print(gg)

gg+
  coord_cartesian(ylim=c(-1.6, -1.4), xlim=c(0, 2.5))

