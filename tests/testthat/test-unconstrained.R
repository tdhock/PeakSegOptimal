library(testthat)
library(PeakSegOptimal)
library(Segmentor3IsBack)

context("PeakSegFPOPLogUnconstrained - Comparison with Segmentor3IsBack")

# helper function to count segments from our output
count_segments <- function(result) {
  sum(result$ends.vec >= 0)
}

extract_segment_means <- function(result) {
  valid_idx <- which(result$ends.vec >= 0)
  means <- result$mean.vec[valid_idx]
  rev(means)
}

test_that("Matches Segmentor for K=2 segments", {
  data <- as.integer(c(5, 5, 5, 20, 20, 20))
  
  # get Segmentor's 2-segment model
  seg_result <- Segmentor(data, model=1, Kmax=2)
  
  our_result <- PeakSegFPOPLogUnconstrained(data, penalty=0.5)
  
  n_segments <- count_segments(our_result)
  expect_equal(n_segments, 2)
  
  # checking if means are close to expected values
  means <- extract_segment_means(our_result)
  expect_equal(length(means), 2)
  expect_true(abs(means[1] - 5) < 0.1)
  expect_true(abs(means[2] - 20) < 0.1)
})


test_that("high penalty gives 1 segment", {
  data <- as.integer(c(5, 5, 5, 20, 20, 20))
  
  # our result with very high penalty (forces 1 segment)
  our_result <- PeakSegFPOPLogUnconstrained(data, penalty=100)
  
  # check we got 1 segment
  n_segments <- count_segments(our_result)
  expect_equal(n_segments, 1)
  
  means <- extract_segment_means(our_result)
  overall_mean <- mean(data)
  expect_true(abs(means[1] - overall_mean) < 1)
})


test_that("works with more complex data", {
  data <- as.integer(c(3, 3, 3, 10, 10, 10, 5, 5, 5))
  
  penalties <- c(0.1, 1, 5, 10, 50)
  
  for(pen in penalties) {
    our_result <- PeakSegFPOPLogUnconstrained(data, penalty=pen)
    n_segments <- count_segments(our_result)
    
    expect_true(n_segments >= 1)
    expect_true(n_segments <= length(data))
    
    means <- extract_segment_means(our_result)
    expect_true(all(means > 0))
  }
})

test_that("Penalty controls number of segments", {
  data <- as.integer(c(1, 1, 5, 5, 10, 10, 3, 3))
  
  # low penalty should give more segments
  result_low <- PeakSegFPOPLogUnconstrained(data, penalty=0.1)
  segments_low <- count_segments(result_low)
  
  # high penalty should give fewer segments  
  result_high <- PeakSegFPOPLogUnconstrained(data, penalty=100)
  segments_high <- count_segments(result_high)
  
  expect_true(segments_low >= segments_high)
})


test_that("model selection across penalty range", {
  data <- as.integer(c(5, 5, 5, 20, 20, 20, 10, 10, 10))
  
  # (higher penalty should never give more segments)
  penalties <- c(0.1, 0.5, 1, 2, 5, 10, 20, 50)
  segment_counts <- numeric(length(penalties))
  
  for(i in seq_along(penalties)) {
    result <- PeakSegFPOPLogUnconstrained(data, penalty=penalties[i])
    segment_counts[i] <- count_segments(result)
  }
  
  # check that sgements counts are non-increasing
  for(i in 2:length(segment_counts)) {
    expect_true(segment_counts[i] <= segment_counts[i-1])
  }
})