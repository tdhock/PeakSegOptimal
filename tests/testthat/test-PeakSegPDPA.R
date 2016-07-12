library(PeakSegDP)
library(coseg)
data.vec <- as.integer(c(1, 10, 14, 13))
fit <- PeakSegPDPA(data.vec, rep(1L, 4), 3L)
library(testthat)
test_that("first segment is OK", {
  cumsum.vec <- cumsum(data.vec)
  n.vec <- seq_along(data.vec)
  mean.vec <- cumsum.vec/n.vec
  expect_equal(fit$mean.mat[1, 1], mean.vec[4])
  for(i in n.vec){
    expected.loss <- PoissonLoss(data.vec[1:i], mean.vec[i])
    expect_equal(fit$cost.mat[1,i], expected.loss)
  }
})

test_that("second segment is OK", {
  expected22 <- PoissonLoss(c(1,10), c(1,10))
  expect_equal(fit$cost.mat[2,2], expected22)
  expected23 <- PoissonLoss(c(1,10,14), c(1,12,12))
  expect_equal(fit$cost.mat[2,3], expected23)
  mean3 <- (10+14+13)/3
  mean24.vec <- c(1, rep(mean3, 3))
  expected24 <- PoissonLoss(data.vec, mean24.vec)
  expect_equal(fit$cost.mat[2,4], expected24)
  expect_equal(fit$ends.mat[2,2], 1)
  expect_equal(fit$mean.mat[2,1:2], c(1, mean3))
})

test_that("third segment is OK", {
  expected33 <- PoissonLoss(c(1,10,14), c(1,12,12))
  expect_equal(fit$cost.mat[3,3], expected33)
  mean3 <- (10+14+13)/3
  mean34.vec <- c(1, rep(mean3, 3))
  expected34 <- PoissonLoss(c(1,10,14,13), mean34.vec)
  expect_equal(fit$cost.mat[3,4], expected34)
  expect_equal(fit$mean.mat[3,], c(1, mean3, mean3))
})

test_that("segment mean 0 before is OK", {
  fit <- PeakSegPDPA(as.integer(c(0, 10, 14, 13)), rep(1L, 4), 3L)
  expect_identical(fit$mean.mat[3,], c(0, 37/3, 37/3))
})

test_that("segment mean 0 after is OK", {
  fit <- PeakSegPDPA(as.integer(c(1, 10, 14, 0)), rep(1L, 4), 3L)
  expect_identical(fit$mean.mat[3,], c(1, 12, 0))
})

test_that("weighted loss same as duplicated loss", {
  fit.id <- PeakSegPDPA(
    as.integer(c(1, 10, 14, 0)), as.integer(c(1, 1, 1, 1)), 3L)
  fit.weighted <- PeakSegPDPA(
    as.integer(c(1, 10, 14, 0)), as.integer(c(2, 1, 1, 1)), 3L)
  fit.duplicated <- PeakSegPDPA(
    as.integer(c(1, 1, 10, 14, 0)), as.integer(c(1, 1, 1, 1, 1)), 3L)
  expect_equal(fit.weighted$cost.mat[, 4], fit.duplicated$cost.mat[, 5])
})

data.vec <- as.integer(c(0, 10, 14, 13))
fit <- PeakSegPDPA(data.vec, rep(1L, 4), 3L)
test_that("segment mean 0 is OK", {
  expect_identical(fit$mean.mat, rbind(
    c(9.25, Inf, Inf),
    c(0, 37/3, Inf),
    c(0, 37/3, 37/3)))
})

data(H3K4me3_XJ_immune_chunk1)
H3K4me3_XJ_immune_chunk1$count <- H3K4me3_XJ_immune_chunk1$coverage
by.sample <- split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
sapply(by.sample, nrow)

## load("H3K4me3_XJ_immune_chunk1_McGill0101.RData")
## timestep <- 422
## total.segments <- 9
## one.name <- "McGill0101"

## library(ggplot2)
## ggplot()+
##   theme_bw()+
##   theme(panel.margin=grid::unit(0, "lines"))+
##   facet_grid(sample.id ~ ., scales="free", labeller=function(df){
##     df$sample.id <- sub("McGill0", "", df$sample.id)
##     df
##   })+
##   geom_step(aes(chromStart/1e3, coverage),
##             data=H3K4me3_XJ_immune_chunk1, color="grey")

max.segments <- 19L
one.name <- "McGill0004"
for(one.name in names(by.sample)){
  one <- by.sample[[one.name]]
  count.vec <- one$coverage
  weight.vec <- with(one, chromEnd-chromStart)
  pdpa <- PeakSegPDPALog(count.vec, weight.vec, max.segments)
  ## Code to compare Log space computation with positive space.
  ## pdpa.orig <- PeakSegPDPA(count.vec, weight.vec, max.segments)
  ## s <- 5
  ## plot(pdpa$cost.mat[s,], pdpa.orig$cost.mat[s,])
  ## plot(pdpa$intervals.mat[s,], pdpa.orig$intervals.mat[s,])
  ## rbind(pdpa$intervals.mat[s,], pdpa.orig$intervals.mat[s,])
  ## library(data.table)
  ## intervals.dt <- data.table(
  ##   col=as.integer(col(pdpa$cost.mat)),
  ##   row=as.integer(row(pdpa$cost.mat)),
  ##   pdpa.cost=as.numeric(pdpa$cost.mat),
  ##   pdpa.orig.cost=as.numeric(pdpa.orig$cost.mat),
  ##   pdpa=as.integer(pdpa$intervals.mat),
  ##   pdpa.orig=as.integer(pdpa.orig$intervals.mat))
  peakseg <- PeakSegDP(one, 9L)
  seg.vec <- seq(1,19,by=2)
  all.loss <- data.frame(
    pdpa=pdpa$cost.mat[seg.vec,length(count.vec)],
    dp=NA,
    row.names=seg.vec)
  is.feasible <- function(loss.vec){
    !any(diff(loss.vec) == 0, na.rm=TRUE)
  }
  all.loss$pdpa.feasible <- apply(pdpa$mean.mat[seg.vec,], 1, is.feasible)
  all.loss[paste(peakseg$error$segments), "dp"] <- peakseg$error$error
  all.loss$pdpa.better <- with(all.loss, dp - pdpa)
  ##print(subset(all.loss, pdpa.feasible & is.na(dp)))
  cost.mat <- rbind(
    pdpa=pdpa$cost.mat[peakseg$error$segments, length(count.vec)],
    cdpa=peakseg$error$error)
  ## cdpa <- cDPA(count.vec, weight.vec, max.segments)
  ## prob.segs <- 9
  ## cost.prob <- data.frame(
  ##   data.i=seq_along(count.vec),
  ##   pdpa=pdpa$cost.mat[prob.segs,],
  ##   cdpa=cdpa$loss[prob.segs,])
  ## cost.prob$should.be.positive <- with(cost.prob, cdpa-pdpa)
  ## prob.labels <- subset(cost.prob, cdpa-pdpa < -1e-10)
  ## gg <- ggplot()+
  ##   geom_line(aes(data.i,cdpa-pdpa),data=cost.prob)+
  ##   geom_point(aes(data.i,cdpa-pdpa),data=prob.labels,shape=1,color="red")
  ## print(gg)
  ##print(one.name)
  ## print(min.diff)
  diff.vec <- apply(cost.mat, 2, diff)
  ##print(diff.vec)
  min.diff <- min(diff.vec)
  expect_gt(min.diff, -1e-10)
}

##library(data.table)
## plot(pdpa$intervals.mat[9,], Rintervals.mat[9,])
## plot(seq_along(count.vec), pdpa$cost.mat[9,]-Rcost.mat[9,])
## plot(seq_along(count.vec), pdpa$intervals.mat[9,]-Rintervals.mat[9,])
## int.dt <- data.table(
##   data.i=seq_along(count.vec),
##   Cintervals=pdpa$intervals.mat[9,],
##   Rintervals=Rintervals.mat[9,])
## int.dt[Cintervals!=Rintervals,]

## ploss <- function(dt, x){
##   ## need to make a new data table, otherwise ifelse may only get one
##   ## element, and return only one element.
##   new.dt <- data.table(dt, x)
##   new.dt[, ifelse(Log==0, 0, Log*log(x)) + Linear*x + Constant]
## }

## plossC <- function(dt, x){
##   ## need to make a new data table, otherwise ifelse may only get one
##   ## element, and return only one element.
##   new.dt <- data.table(dt, x)
##   loss <- new.dt[, Linear*x + Constant]
##   cat(sprintf("Rloss before adding log term=%a\n", loss))
##   log.mean.only <- log(x)
##   cat(sprintf("log.mean.only=%a\n", log.mean.only))
##   log.coef.only <- as.numeric(new.dt$Log)
##   cat(sprintf("log.coef.only=%a\n", log.coef.only))
##   product <- log.mean.only * log.coef.only
##   cat(sprintf("product=%a\n", product))
##   loss + ifelse(new.dt$Log==0, 0, product)
## }

## pderiv <- function(dt, x){
##   dt[, Linear+Log/x]
## }

## getLines <- function(dt, n.data=100){
##   line.list <- list()
##   for(piece.i in 1:nrow(dt)){
##     piece <- dt[piece.i,]
##     mean.vec <- piece[, seq(min.mean, max.mean, l=n.data)]
##     line.list[[piece.i]] <- data.table(
##       piece.i,
##       piece,
##       mean=mean.vec,
##       cost=ploss(piece, mean.vec))
##   }
##   do.call(rbind, line.list)
## }

## getMinMean <- function(dt){
##   dt[, -Log/Linear]
## }

## AddFuns <- function(dt1, dt2){
##   if(nrow(dt1)==0)return(dt1)
##   if(nrow(dt2)==0)return(dt2)
##   ## NOTE asymmetrey of dt1 and dt2 -- dt1 is used for data.i
##   i1 <- 1
##   i2 <- 1
##   new.dt.list <- list()
##   while(i1 <= nrow(dt1) && i2 <= nrow(dt2)){
##     row1 <- dt1[i1,]
##     row2 <- dt2[i2,]
##     this.min <- max(row1$min.mean, row2$min.mean)
##     this.max <- if(row1$max.mean < row2$max.mean){
##       i1 <- i1+1
##       row1$max.mean
##     }else{
##       i2 <- i2+1
##       row2$max.mean
##     }
##     new.dt.list[[paste(i1, i2)]] <- data.table(
##       Linear=row1$Linear+row2$Linear,
##       Log=row1$Log+row2$Log,
##       Constant=row1$Constant+row2$Constant,
##       min.mean=this.min,
##       max.mean=this.max,
##       data.i=row1$data.i)
##     this.min <- this.max
##   }
##   do.call(rbind, new.dt.list)
## }

## less.more.min.list <- list(
##   less=function(dt){
##     new.dt.list <- list()
##     prev.min.cost <- NULL
##     prev.data.i <- NULL
##     row.i <- 1
##     prev.min.mean <- dt$min.mean[1]
##     while(row.i <= nrow(dt)){
##       this.row <- dt[row.i,]
##       if(is.null(prev.min.cost)){
##         ## Look for min achieved in this interval.
##         mu <- getMinMean(this.row)
##         cat(sprintf("RgetMinMean=%a\n", mu))
##         if(mu <= this.row$min.mean){
##           ## The minimum is achieved before this interval, so this
##           ## function is always increasing in this interval. We don't
##           ## need to store it.
##           prev.min.cost <- ploss(this.row, this.row$min.mean)
##           prev.data.i <- this.row$data.i
##         }else if(mu < this.row$max.mean){
##           ## Minimum in this interval.
##           new.row <- this.row
##           new.row$min.mean <- prev.min.mean
##           new.row$max.mean <- mu
##           new.dt.list[[paste(row.i)]] <- new.row
##           prev.min.mean <- mu
##           prev.min.cost <- plossC(this.row, mu)
##           cat(sprintf("prev.min.cost=%a\n", prev.min.cost))
##           prev.data.i <- this.row$data.i
##         }else{
##           ## Minimum after this interval, so this function is
##           ## decreasing on this entire interval, and so we can just
##           ## store it as is.
##           new.row <- this.row
##           new.row$min.mean <- prev.min.mean
##           new.dt.list[[paste(row.i)]] <- new.row
##           prev.min.mean <- this.row$max.mean
##         }
##       }else{
##         ## Look for a function with prev.min.cost.
##         if(this.row$Log==0){
##           ## degenerate linear case.
##           if(this.row$Linear < 0){
##             ## decreasing linear function
##             stop("this should never happen")
##           }else{
##             ##increasing linear function, so will not intersect the
##             ##constant below.
##           }
##         }else{
##           discriminant <-
##             this.row[, Linear/Log*exp((prev.min.cost-Constant)/Log)]
##           if(-1/exp(1) < discriminant){
##             ## Since the Log constant is negative, the principal branch
##             ## W(,0) results in the smaller of the two mean values.
##             mu <- this.row[, Log/Linear*LambertW::W(discriminant, 0)]
##             if(this.row[, min.mean < mu & mu < max.mean]){
##               new.dt.list[[paste(row.i, "constant")]] <- data.table(
##                 Linear=0,
##                 Log=0,
##                 Constant=prev.min.cost,
##                 min.mean=prev.min.mean,
##                 max.mean=mu,
##                 data.i=prev.data.i)
##               prev.min.cost <- NULL
##               prev.min.mean <- mu
##               row.i <- row.i-1
##             }
##           }#if(there are two roots
##         }#if(degenerate linear) else
##       }#if(looking for a min)else
##       row.i <- row.i+1
##     }
##     if(!is.null(prev.data.i)){
##       new.dt.list[["last"]] <- data.table(
##         Linear=0,
##         Log=0,
##         Constant=prev.min.cost,
##         min.mean=prev.min.mean,
##         max.mean=this.row$max.mean,
##         data.i=prev.data.i)
##     }
##     do.call(rbind, new.dt.list)
##   }, more=function(dt){
##     new.dt.list <- list()
##     prev.min.cost <- NULL
##     prev.data.i <- NULL
##     row.i <- nrow(dt)
##     prev.max.mean <- dt$max.mean[row.i]
##     while(1 <= row.i){
##       this.row <- dt[row.i,]
##       if(is.null(prev.min.cost)){
##         ## Look for min achieved in this interval.
##         mu <- getMinMean(this.row)
##         if(this.row$max.mean <= mu){
##           ## The minimum is achieved after this interval, so this
##           ## function is always decreasing in this interval. We don't
##           ## need to store it.
##           prev.min.cost <- ploss(this.row, this.row$max.mean)
##           prev.data.i <- this.row$data.i
##         }else if(this.row$min.mean < mu){
##           ## Minimum in this interval.
##           new.row <- this.row
##           new.row$max.mean <- prev.max.mean
##           new.row$min.mean <- mu
##           new.dt.list[[paste(row.i)]] <- new.row
##           prev.max.mean <- mu
##           prev.min.cost <- ploss(this.row, mu)
##           prev.data.i <- this.row$data.i
##         }else{
##           ## Minimum before this interval, so this function is
##           ## increasing on this entire interval, and so we can just
##           ## store it as is.
##           new.row <- this.row
##           new.row$max.mean <- prev.max.mean
##           new.dt.list[[paste(row.i)]] <- new.row
##           prev.max.mean <- this.row$min.mean
##         }
##       }else{
##         ## Look for a function with prev.min.cost.
##         mu <- if(this.row$Log==0){
##           ## degenerate case where the function is linear.
##           this.row[, (prev.min.cost - Constant)/Linear]
##         }else{
##           discriminant <-
##             this.row[, Linear/Log*exp((prev.min.cost-Constant)/Log)]
##           if(-1/exp(1) < discriminant){
##             ## Since the Log constant is negative, the non-principal
##             ## branch W(,-1) results in the larger of the two mean
##             ## values.
##             browser(expr=!is.finite(discriminant))
##             this.row[, Log/Linear*LambertW::W(discriminant, -1)]
##           }#if(-1/e < discriminant
##         }
##         if(is.numeric(mu) && this.row[, min.mean < mu & mu < max.mean]){
##           new.dt.list[[paste(row.i, "constant")]] <- data.table(
##             Linear=0,
##             Log=0,
##             Constant=prev.min.cost,
##             min.mean=mu,
##             max.mean=prev.max.mean,
##             data.i=prev.data.i)
##           prev.min.cost <- NULL
##           prev.max.mean <- mu
##           row.i <- row.i+1
##         }#if(mu in interval
##       }#if(is.null(prev.min.cost)else
##       row.i <- row.i-1
##     }
##     if(!is.null(prev.data.i)){
##       new.dt.list[["last"]] <- data.table(
##         Linear=0,
##         Log=0,
##         Constant=prev.min.cost,
##         min.mean=this.row$min.mean,
##         max.mean=prev.max.mean,
##         data.i=prev.data.i)
##     }
##     do.call(rbind, rev(new.dt.list))
##   })

## Minimize <- function(dt, from=min(dt$min.mean), to=max(dt$max.mean)){
##   stopifnot(from < to)
##   is.before <- dt$max.mean < from
##   is.after <- to < dt$min.mean
##   feasible <- dt[!(is.before | is.after),]
##   feasible$min.mean[1] <- from
##   feasible$max.mean[nrow(feasible)] <- to
##   feasible$fun.min.mean <- getMinMean(feasible)
##   feasible[, min.cost.mean := ifelse(
##     fun.min.mean < min.mean, min.mean,
##     ifelse(max.mean < fun.min.mean,
##            max.mean,
##            fun.min.mean))]
##   feasible[, min.cost := ploss(feasible, min.cost.mean)]
##   feasible[which.min(min.cost),]
## }

## sameFuns <- function(row1, row2){
##   row1$Linear==row2$Linear &&
##     row1$Log==row2$Log &&
##     row1$Constant==row2$Constant
## }

## CompareRows <- function(dt1, dt2, i1, i2){
##   row1 <- dt1[i1,]
##   row2 <- dt2[i2,]
##   if(row1$min.mean < row2$min.mean){
##     prev2 <- dt2[i2-1, ]
##     same.at.left <- sameFuns(prev2, row1)
##     last.min.mean <- row2$min.mean
##   }else{
##     last.min.mean <- row1$min.mean
##     prev1 <- dt1[i1-1, ]
##     same.at.left <- if(row2$min.mean < row1$min.mean){
##       sameFuns(prev1, row2)
##     }else{
##       if(i1==1 && i2==1){
##         FALSE
##       }else{
##         prev2 <- dt2[i2-1, ]
##         sameFuns(prev1, prev2)
##       }
##     }
##   }
##   if(row1$max.mean < row2$max.mean){
##     next1 <- dt1[i1+1,]
##     same.at.right <- sameFuns(next1, row2)
##     first.max.mean <- row1$max.mean
##   }else{
##     first.max.mean <- row2$max.mean
##     same.at.right <- if(row2$max.mean < row1$max.mean){
##       next2 <- dt2[i2+1,]
##       sameFuns(row1, next2)
##     }else{
##       if(i1==nrow(dt1) && i2==nrow(dt2)){
##         FALSE
##       }else{
##         next1 <- dt1[i1+1,]
##         next2 <- dt2[i2+1,]
##         sameFuns(next1, next2)
##       }
##     }
##   }
##   stopifnot(last.min.mean < first.max.mean)
##   is.same <- sameFuns(row1, row2)
##   if(is.same){
##     ## The functions are exactly equal over the entire interval so we
##     ## can return either one of them.
##     row1$min.mean <- last.min.mean
##     row1$max.mean <- first.max.mean
##     return(row1)
##   }
##   cat(sprintf("row1$Constant=%a\nrow2$Constant=%a\n",
##               row1$Constant,              row2$Constant))
##   row.diff <- row1-row2
##   row.diff$min.mean <- last.min.mean
##   row.diff$max.mean <- first.max.mean
##   if(row.diff$Log==0){
##     if(row.diff$Linear==0){
##       ## They are offset by a constant.
##       new.row <- if(row.diff$Constant < 0)row1 else row2
##       new.row$min.mean <- last.min.mean
##       new.row$max.mean <- first.max.mean
##       return(new.row)
##     }
##     if(row.diff$Constant==0){
##       ## The only difference is the Linear coef.
##       new.row <- if(row.diff$Linear < 0)row1 else row2
##       new.row$min.mean <- last.min.mean
##       new.row$max.mean <- first.max.mean
##       return(new.row)
##     }
##     mean.at.equal.cost <- row.diff[, -Constant/Linear]
##     root.in.interval <-
##       last.min.mean < mean.at.equal.cost &&
##       mean.at.equal.cost < first.max.mean
##     if(root.in.interval){
##       new.rows <- if(0 < row.diff$Linear){
##         rbind(row1, row2)
##       }else{
##         rbind(row2, row1)
##       }
##       new.rows$min.mean <- c(last.min.mean, mean.at.equal.cost)
##       new.rows$max.mean <- c(mean.at.equal.cost, first.max.mean)
##       return(new.rows)
##     }else{
##       new.row <- if(mean.at.equal.cost < last.min.mean)row1 else row2
##       new.row$min.mean <- last.min.mean
##       new.row$max.mean <- first.max.mean
##       return(new.row)
##     }
##   }
##   ## cost2.left <- ploss(row2, last.min.mean)
##   ## cost1.left <- ploss(row1, last.min.mean)
##   ## cost1.right <- ploss(row1, first.max.mean)
##   ## cost2.right <- ploss(row2, first.max.mean)
##   cost.diff.right <- ploss(row.diff, first.max.mean)
##   cost.diff.left <- ploss(row.diff, last.min.mean)
##   discriminant <- row.diff[, Linear/Log*exp(-Constant/Log)]
##   str(row.diff)
##   cat(sprintf("Rdiscriminant=%a\nConstant=%a\n", discriminant, row.diff$Constant))
##   ## The discriminant could be -Inf, if the exp argument is larger
##   ## than about 700.
##   two.roots <- -1/exp(1) < discriminant
##   if(!same.at.right){
##     row1.min.on.right <- cost.diff.right < 0
##   }else{
##     ## They are equal on the right limit, so use the first and second
##     ## derivatives to see which is minimal just before the right
##     ## limit. Do we need to check if they intersect before the right
##     ## limit? Only if the sign of the slope of the more curved
##     ## function is positive. And in that case we just need to check
##     ## the smaller root.
##     deriv1.right <- pderiv(row1, first.max.mean)
##     deriv2.right <- pderiv(row2, first.max.mean)
##     sign1 <- sign(deriv1.right)
##     sign2 <- sign(deriv2.right)
##     maybe.cross <-
##       (row2$Log < row1$Log && 0 < sign2) ||
##       (row1$Log < row2$Log && 0 < sign1)
##     if(two.roots && maybe.cross){
##       ## There could be a crossing point to the left.
##       root.left <- row.diff[, Log/Linear*LambertW::W(discriminant, 0)]
##       cat(sprintf("smaller_mean=%a\n", root.left))
##       mean.at.equal.cost <- root.left
##       in.interval <-
##         last.min.mean < mean.at.equal.cost &
##         mean.at.equal.cost < first.max.mean
##       if(in.interval){
##         new.rows <- if(cost.diff.left < 0){
##           rbind(row1, row2)
##         }else{
##           rbind(row2, row1)
##         }
##         new.rows$min.mean <- c(last.min.mean, mean.at.equal.cost)
##         new.rows$max.mean <- c(mean.at.equal.cost, first.max.mean)
##         return(new.rows)
##       }
##     }
##     row1.min.before.right <- cost.diff.left < 0
##     this.row <- if(row1.min.before.right)row1 else row2
##     this.row$min.mean <- last.min.mean
##     this.row$max.mean <- first.max.mean
##     return(this.row)
##   }
##   if(!same.at.left){
##     row1.min.on.left <- cost.diff.left < 0
##   }else{
##     ## Equal on the left.
##     deriv1.left <- pderiv(row1, last.min.mean)
##     deriv2.left <- pderiv(row2, last.min.mean)
##     sign1 <- sign(deriv1.left)
##     sign2 <- sign(deriv2.left)
##     maybe.cross <-
##       (row2$Log < row1$Log && sign2 < 0) ||
##       (row1$Log < row2$Log && sign1 < 0)
##     if(two.roots && maybe.cross){
##       ## There could be a crossing point to the right.
##       root.right <- row.diff[, Log/Linear*LambertW::W(discriminant, -1)]
##       cat(sprintf("larger_mean=%a\n", root.right))
##       mean.at.equal.cost <- root.right
##       in.interval <-
##         last.min.mean < mean.at.equal.cost &
##         mean.at.equal.cost < first.max.mean
##       if(in.interval){
##         new.rows <- if(cost.diff.right < 0){
##           rbind(row2, row1)
##         }else{
##           rbind(row1, row2)
##         }
##         new.rows$min.mean <- c(last.min.mean, mean.at.equal.cost)
##         new.rows$max.mean <- c(mean.at.equal.cost, first.max.mean)
##         return(new.rows)
##       }
##     }
##     row1.min.after.left <- cost.diff.right < 0
##     this.row <- if(row1.min.after.left)row1 else row2
##     this.row$min.mean <- last.min.mean
##     this.row$max.mean <- first.max.mean
##     return(this.row)
##   }
##   ## The only remaining case is that the curves are equal neither on
##   ## the left nor on the right of the interval. However they may be
##   ## equal inside the interval, so let's check for that.
##   mean.in.interval <- if(two.roots){
##     root.right <- row.diff[, Log/Linear*LambertW::W(discriminant, -1)]
##     root.left <- row.diff[, Log/Linear*LambertW::W(discriminant, 0)]
##     cat(sprintf("smaller_mean=%a\nlarger_mean=%a\n", root.left, root.right))
##     mean.at.equal.cost <- sort(c(root.left, root.right))
##     in.interval <-
##       last.min.mean < mean.at.equal.cost &
##       mean.at.equal.cost < first.max.mean
##     mean.at.equal.cost[in.interval]
##   }
##   new.intervals <- if(length(mean.in.interval)==1){
##     if(row1.min.on.left)rbind(row1,row2) else rbind(row2, row1)
##   }else if(length(mean.in.interval)==2){
##     if(row1.min.on.left){
##       rbind(row1, row2, row1)
##     }else{
##       rbind(row2, row1, row2)
##     }
##   }else{
##     ## functions do not cross in this interval.
##     if(row1.min.on.right)row1 else row2
##   }
##   new.intervals$min.mean <- c(last.min.mean, mean.in.interval)
##   new.intervals$max.mean <- c(mean.in.interval, first.max.mean)
##   new.intervals
## }

## MinEnvelope <- function(dt1, dt2){
##   stopifnot(dt1[, min.mean < max.mean])
##   stopifnot(dt2[, min.mean < max.mean])
##   stopifnot(dt1[1, min.mean]==dt2[1, min.mean])
##   stopifnot(dt1[.N, max.mean]==dt2[.N, max.mean])
##   i1 <- 1
##   i2 <- 1
##   new.dt.list <- list()
##   row.to.add <- NULL
##   while(i1 <= nrow(dt1) && i2 <= nrow(dt2)){
##     row1 <- dt1[i1,]
##     row2 <- dt2[i2,]
##     new.rows <- CompareRows(dt1, dt2, i1, i2)
##     for(row.i in 1:nrow(new.rows)){
##       new.row <- new.rows[row.i,]
##       if(is.null(row.to.add)){
##         row.to.add <- new.row
##       }else{
##         if(sameFuns(new.row, row.to.add)){
##           ## this is the same function as the previous one, so just make
##           ## it extend further to the right.
##           row.to.add$max.mean <- new.row$max.mean
##         }else{
##           ## new.row is a different function than the previous
##           ## row.to.add, so now is the time to store it and move on.
##           new.dt.list[[paste(i1, i2, row.i)]] <- row.to.add
##           row.to.add <- new.row
##         }
##       }
##     }
##     if(row1$max.mean == new.row$max.mean){
##       i1 <- i1+1
##     }
##     if(row2$max.mean == new.row$max.mean){
##       i2 <- i2+1
##     }
##   }
##   new.dt.list[["last"]] <- row.to.add
##   do.call(rbind, new.dt.list)
## }

## count.vec <- data.vec
## weight.vec <- 1
## input.dt <- data.table(count=count.vec, weight=weight.vec)
## max.segments <- 3L
## stopifnot(is.data.table(input.dt))
## stopifnot(2 <= max.segments && max.segments <= nrow(input.dt))
## stopifnot(c("weight", "count") %in% names(input.dt))
## min.mean <- min(input.dt$count)
## max.mean <- max(input.dt$count)
## gamma.dt <- input.dt[, data.table(
##   Linear=weight,
##   Log=-count*weight,
##   Constant=0)]
## C1.dt <- cumsum(gamma.dt)
## gamma.dt$min.mean <- C1.dt$min.mean <- min.mean
## gamma.dt$max.mean <- C1.dt$max.mean <- max.mean
## gamma.dt$data.i <- C1.dt$data.i <- 0
## cost.models.list <- list()
## for(data.i in 1:nrow(C1.dt)){
##   cost.models.list[[paste(1, data.i)]] <- C1.dt[data.i,]
## }
## stopifnot(max.segments <= nrow(input.dt))
## intervals.list <- list()
## for(total.segments in 2:max.segments){
##   prev.cost.model <-
##     cost.models.list[[paste(total.segments-1, total.segments-1)]]
##   if(total.segments %% 2){
##     min.fun.name <- "more"
##   }else{
##     min.fun.name <- "less"
##   }
##   min.fun <- less.more.min.list[[min.fun.name]]
##   first.min <- min.fun(prev.cost.model)
##   first.data <- gamma.dt[total.segments,]
##   first.data$data.i <- total.segments-1
##   cost.model <- AddFuns(first.data, first.min)
##   intervals.list[[paste(total.segments, total.segments)]] <- data.table(
##     total.segments, timestep=total.segments, intervals=nrow(cost.model))
##   cost.models.list[[paste(total.segments, total.segments)]] <- cost.model
##   for(timestep in (total.segments+1):length(input.dt$count)){
##     ## problem at timestep=423? infinite cost on left?
##     cat(sprintf("%4d / %4d segments %4d / %4d data points %d intervals\n",
##                 total.segments, max.segments, timestep, length(input.dt$count),
##                 nrow(cost.model)))
##     prev.cost.model <- cost.models.list[[paste(total.segments-1, timestep-1)]]
##     compare.cost <- min.fun(prev.cost.model)
##     compare.cost$data.i <- timestep-1
##     cost.model <- cost.models.list[[paste(total.segments, timestep-1)]]
##     one.env <- MinEnvelope(compare.cost, cost.model)
##     gg <- ggplot()+
##       ggtitle(paste(total.segments, "segments,", timestep, "data points"))+
##       geom_line(aes(mean, cost),
##                 data=getLines(one.env),
##                 color="grey",
##                 size=2,
##                 alpha=0.5)+
##       geom_line(aes(mean, cost,
##                     color="compare.cost",
##                     group=piece.i),
##                 data=getLines(compare.cost))+
##       geom_line(aes(mean, cost,
##                     group=piece.i,
##                     color="cost.model"),
##                 data=getLines(cost.model))
##     print(gg)
##     ## browser()
##     stopifnot(one.env[, min.mean < max.mean])
##     ## Now that we are done with this step, we can perform the
##     ## recursion by setting the new model of the cost to the min
##     ## envelope, plus a new data point.
##     new.cost.model <- AddFuns(one.env, gamma.dt[timestep,])
##     intervals.list[[paste(total.segments, timestep)]] <- data.table(
##       total.segments, timestep, intervals=nrow(new.cost.model))
##     cost.models.list[[paste(total.segments, timestep)]] <-
##       new.cost.model
##   }#for(timestep
## }#for(total.segments
## minima.list <- list()
## cost.list <- list()
## for(total.segments in seq(1, max.segments, by=2)){
##   peaks <- (total.segments-1)/2
##   timestep <- length(input.dt$count)
##   cat(sprintf(
##     "decoding %4d / %4d segments\n",
##     total.segments, max.segments))
##   data.i <- timestep
##   seg.i <- total.segments
##   no.constraint <- data.table(
##     min.mean,
##     max.mean,
##     data.i=NA)
##   constraint <- no.constraint
##   segment.end <- timestep
##   while(0 < seg.i && length(data.i)==1){
##     unconstrained.fun <- cost.models.list[[paste(seg.i, data.i)]]
##     min.dt <- Minimize(
##       unconstrained.fun,
##       constraint$min.mean,
##       constraint$max.mean)
##     min.dt$segment.end <- segment.end
##     min.dt$peaks <- peaks
##     min.dt$total.segments <- total.segments
##     min.dt$seg.i <- seg.i
##     min.dt[, segment.start := ifelse(seg.i==1, 1, 1+data.i)]
##     segment.end <- min.dt$data.i
##     if(seg.i==total.segments){
##       cost.row <- min.dt
##       cost.row$constraint <- "inactive"
##     }
##     min.dt$constraint <- if(min.dt[, fun.min.mean == min.cost.mean]){
##       "inactive"
##     }else{
##       cost.row$constraint <- "active"
##       "active"
##     }
##     minima.list[[paste(total.segments, seg.i)]] <- min.dt
##     constraint <- no.constraint
##     constraint.side <- if(seg.i %% 2){
##       constraint$min.mean <- min.dt$min.cost.mean
##     }else{
##       constraint$max.mean <- min.dt$min.cost.mean
##     }
##     data.i <- min.dt$data.i
##     seg.i <- seg.i-1
##   }#while(...
##   cost.list[[paste(total.segments)]] <- cost.row
## }#for(total.segments

## Rcost.mat <- matrix(NA, max.segments, length(count.vec))
## Rintervals.mat <- matrix(NA, max.segments, length(count.vec))
## for(total.segments in 1:max.segments){
##   for(data.i in total.segments:length(count.vec)){
##     cat(sprintf("segments=%d data=%d\n", total.segments, data.i))
##     fun.dt <- cost.models.list[[paste(total.segments, data.i)]]
##     min.dt <- Minimize(fun.dt)
##     Rcost.mat[total.segments, data.i] <- min.dt$min.cost
##     Rintervals.mat[total.segments, data.i] <- nrow(fun.dt)
##   }
## }

## rbind(
##   pdpa=pdpa$cost.mat[9, length(count.vec)],
##   cdpa=peakseg$error$error[5],
##   Rpdpa=Rcost.mat[9, length(count.vec)])
