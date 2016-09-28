library(testthat)
context("cosegData")
library(coseg)

data(H3K36me3_AM_immune_McGill0079_chr3_60000_66170270, package="cosegData")

test_that("big data and penalty do not crash", {
  data.list <- H3K36me3_AM_immune_McGill0079_chr3_60000_66170270
  fit <- with(data.list, PeakSegFPOPchrom(coverage[1:100,], as.numeric(penalty)))
})

