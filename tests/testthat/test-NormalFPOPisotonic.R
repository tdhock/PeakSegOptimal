library(testthat)
context("NormalFPOPisotonic")
library(PeakSegOptimal)

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

test_that("strictly decreasing data gives flat mean (PAVA pooling)", {
  y <- c(5, 4, 3, 2, 1)
  result <- NormalFPOPisotonic(y, penalty = 0)
  fitted <- expand_means(result, length(y))
  expect_equal(fitted, rep(mean(y), length(y)), tolerance = 1e-6)
})

test_that("penalty=0 matches isoreg on random data", {
  set.seed(42)
  y <- c(5, 3, 4, 1, 2, 8, 7, 9)
  result <- NormalFPOPisotonic(y, penalty = 0)
  fitted <- expand_means(result, length(y))
  iso_fitted <- as.numeric(isoreg(y)$yf)
  expect_equal(fitted, iso_fitted, tolerance = 1e-5,
               info = "random 8-element vector")
})

test_that("penalty=0 matches isoreg on already sorted data", {
  y <- c(1, 2, 3, 4, 5)
  result <- NormalFPOPisotonic(y, penalty = 0)
  fitted <- expand_means(result, length(y))
  iso_fitted <- as.numeric(isoreg(y)$yf)
  expect_equal(fitted, iso_fitted, tolerance = 1e-5,
               info = "already non-decreasing")
})

test_that("penalty=0 matches isoreg on strictly decreasing data", {
  y <- c(10, 8, 6, 4, 2)
  result <- NormalFPOPisotonic(y, penalty = 0)
  fitted <- expand_means(result, length(y))
  iso_fitted <- as.numeric(isoreg(y)$yf)
  expect_equal(fitted, iso_fitted, tolerance = 1e-5,
               info = "strictly decreasing")
})

test_that("penalty=0 matches isoreg on longer random data", {
  set.seed(2026)
  y <- rnorm(50, mean = cumsum(rnorm(50, 0, 0.3)))
  result <- NormalFPOPisotonic(y, penalty = 0)
  fitted <- expand_means(result, length(y))
  iso_fitted <- as.numeric(isoreg(y)$yf)
  expect_equal(fitted, iso_fitted, tolerance = 1e-5,
               info = "50-point random walk")
})

test_that("n=2 works correctly", {
  y <- c(1.0, 5.0)
  r_low <- NormalFPOPisotonic(y, penalty = 0)
  expect_equal(r_low$n.segments, 2L)
  expect_equal(as.numeric(r_low$seg.mean), c(1, 5), tolerance = 1e-6)
  r_high <- NormalFPOPisotonic(y, penalty = 1e6)
  expect_equal(r_high$n.segments, 1L)
  expect_equal(r_high$seg.mean[1], mean(y), tolerance = 1e-6)
})

test_that("means are always non-decreasing", {
  set.seed(7)
  y <- rnorm(30, mean = 0, sd = 5)
  for (pen in c(0, 1, 5, 20, 100)) {
    result <- NormalFPOPisotonic(y, penalty = pen)
    if (result$n.segments > 1) {
      expect_true(all(diff(result$seg.mean) >= -1e-8),
                  info = paste("non-decreasing at penalty =", pen))
    }
  }
})

test_that("high penalty gives 1 segment at overall mean", {
  set.seed(99)
  y <- rnorm(40, mean = 5)
  result <- NormalFPOPisotonic(y, penalty = 1e8)
  expect_equal(result$n.segments, 1L)
  expect_equal(result$seg.mean[1], mean(y), tolerance = 1e-6)
})

test_that("already isotonic data is not altered", {
  y <- c(1.0, 2.0, 3.0, 4.0, 5.0)
  result <- NormalFPOPisotonic(y, penalty = 0)
  fitted <- expand_means(result, length(y))
  iso_fitted <- as.numeric(isoreg(y)$yf)
  expect_equal(fitted, iso_fitted, tolerance = 1e-6,
               info = "isotonic data, penalty=0")
})

if(requireNamespace("fpop", quietly = TRUE)) {

  test_that("matches Fpop when unconstrained solution is non-decreasing", {
    cases <- list(
      list(y = c(rep(2, 20), rep(5, 20), rep(10, 20)), pens = c(1, 5, 10)),
      list(y = c(rep(0, 10), rep(3, 10), rep(7, 10), rep(15, 10)), pens = c(1, 5, 10)),
      list(y = c(rep(1, 30), rep(8, 30)), pens = c(1, 5, 20))
    )
    for (cc in cases) {
      for (pen in cc$pens) {
        fpop_fit <- fpop::Fpop(cc$y, lambda = pen)
        fpop_seg_means <- unique(fpop_fit$signal)
        if (all(diff(fpop_seg_means) >= -1e-8)) {
          iso_fit <- NormalFPOPisotonic(cc$y, penalty = pen)
          iso_fitted <- expand_means(iso_fit, length(cc$y))
          expect_equal(iso_fitted, as.numeric(fpop_fit$signal),
                       tolerance = 1e-5,
                       info = paste("n =", length(cc$y), "pen =", pen))
        }
      }
    }
  })

} else {
  message("fpop not installed, skipping Fpop comparison tests")
}

test_that("penalty=0 matches isoreg on 200 random trials (N=10)", {
  set.seed(20260222)
  n_fail <- 0L
  for (trial in seq_len(200)) {
    y <- switch(((trial - 1L) %% 3L) + 1L,
      rnorm(10), round(runif(10, 0, 10)), rnorm(10, sd = 1e3))
    if (length(unique(y)) < 2L) next
    result <- NormalFPOPisotonic(y, penalty = 0)
    fitted <- expand_means(result, length(y))
    iso <- as.numeric(isoreg(y)$yf)
    if (max(abs(fitted - iso)) > 1e-5) n_fail <- n_fail + 1L
  }
  expect_equal(n_fail, 0L)
})
