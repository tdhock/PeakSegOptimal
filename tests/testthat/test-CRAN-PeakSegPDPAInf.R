library(testthat)
context("PeakSegPDPAInf")
library(PeakSegOptimal)
data("H3K4me3_XJ_immune_chunk1", envir=environment())
by.sample <-
  split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
n.data.vec <- sapply(by.sample, nrow)
one <- by.sample[[1]]
count.vec <- one$coverage
weight.vec <- with(one, chromEnd-chromStart)
max.segments <- 19L
library(data.table)
ic.list <- list()
for(fun.name in c("PeakSegPDPA", "PeakSegPDPAInf")){
  fun <- get(fun.name)
  fit <- fun(count.vec, weight.vec, max.segments)
  ic.list[[fun.name]] <- data.table(
    fun.name,
    segments=as.numeric(row(fit$intervals.mat)),
    data=as.numeric(col(fit$intervals.mat)),
    cost=as.numeric(fit$cost.mat),
    intervals=as.numeric(fit$intervals.mat))
}
all.ic <- do.call(rbind, ic.list)
ic <- subset(all.ic, 0 < intervals)
intervals <- dcast(ic, data + segments ~ fun.name, value.var="intervals")
cost <- dcast(ic, data + segments ~ fun.name, value.var="cost")
not.equal <- subset(cost, PeakSegPDPA != PeakSegPDPAInf)
test_that("same cost for both PDPA algos", {
  expect_equal(nrow(not.equal), 0)
})

pos <- 1:3
rep.df <- data.frame(
  count=5L,
  chromStart=pos,
  chromEnd=pos+1L)
test_that("PeakSegPDPAInf is fine for all same data", {
  L <- with(rep.df, PeakSegPDPAInf(count, chromEnd-chromStart, 3L))
  expect_equal(L$mean.mat[1,], c(5, Inf, Inf))
  expect_equal(L$mean.mat[2,], c(5, 5, Inf))
  expect_equal(L$mean.mat[3,], c(5, 5, 5))
  expect_equal(L$ends.mat[1,], c(0, 0, 0))
  expect_equal(L$ends.mat[2,], c(0, 2, 0))
  expect_equal(L$ends.mat[3,], c(0, 1, 2))
})
