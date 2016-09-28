library(testthat)
context("cosegData")
library(coseg)

obj.name <- data(H3K36me3_AM_immune_McGill0079_chr3_60000_66170270, package="cosegData")
test_that("big data and penalty do not crash", {
  data.list <- get(obj.name)
  fit <- with(data.list, PeakSegFPOPchrom(coverage[1:100,], as.numeric(penalty)))
})

obj.name <- data(H3K4me3_TDH_other_McGill0014_chr5_10000_17530657, package="cosegData")
test_that("bigger data and penalty do not crash", {
  data.list <- get(obj.name)
  fit <- with(data.list, PeakSegFPOPchrom(coverage[1:10000,], as.numeric(penalty)))
})
