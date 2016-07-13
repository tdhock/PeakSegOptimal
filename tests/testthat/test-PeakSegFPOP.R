library(PeakSegDP)
library(coseg)
data.vec <- as.integer(c(1, 10, 14, 13))
(fit <- PeakSegFPOP(data.vec, rep(1L, 4), 0))
(fit <- PeakSegFPOP(data.vec, rep(1L, 4), 1e6))
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
