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
  exp.mat <- rbind(
    c(9.25, Inf, Inf),
    c(0, 37/3, Inf),
    c(0, 37/3, 37/3))
  expect_equal(fit$mean.mat, exp.mat)
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
  pdpa <- PeakSegPDPA(count.vec, weight.vec, max.segments)
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

