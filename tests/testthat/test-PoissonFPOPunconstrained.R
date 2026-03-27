library(testthat)
context("PoissonFPOPunconstrained")
library(PeakSegOptimal)

fpop_to_breaks <- function(result, n) {
  K <- result$n.segments
  if(K == 1) return(n)
  c(result$seg.end[2:K], n)
}

expand_means <- function(result, n) {
  fitted <- numeric(n)
  seg.start <- 1L
  for(k in seq_len(result$n.segments)) {
    if(k < result$n.segments) {
      seg.stop <- result$seg.end[k + 1L]
    } else {
      seg.stop <- n
    }
    fitted[seg.start:seg.stop] <- result$seg.mean[k]
    seg.start <- seg.stop + 1L
  }
  fitted
}

test_that("penalty=0 gives every point its own segment", {
  data_vec <- as.integer(c(3, 10, 3, 10, 3))
  fpop <- PoissonFPOPunconstrained(data_vec, penalty = 0.0)
  expect_identical(as.integer(fpop$n.segments), length(data_vec))
})

test_that("n=2 works for both 1-segment and 2-segment outcomes", {
  data_vec <- as.integer(c(1, 100))
  fpop2 <- PoissonFPOPunconstrained(data_vec, penalty = 0.0)
  expect_identical(as.integer(fpop2$n.segments), 2L)
  fpop1 <- PoissonFPOPunconstrained(data_vec, penalty = 1e10)
  expect_identical(as.integer(fpop1$n.segments), 1L)
})

test_that("rejects negative data", {
  expect_error(PoissonFPOPunconstrained(as.integer(c(-1, 2, 3)), penalty = 5))
})
