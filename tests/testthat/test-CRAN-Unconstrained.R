library(testthat)
test_that("data_index is integer", {
  data.vec <- as.integer(c(1, 10, 14, 13))
  fit <- PeakSegOptimal::UnconstrainedFPOP(data.vec, penalty=0, verbose_file = tempfile())
  expect_identical(fit$index_dt$data_index, 0:3)
})
