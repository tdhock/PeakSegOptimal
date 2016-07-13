library(coseg)
data.vec <- as.integer(c(1, 10, 14, 13))
fit <- PeakSegFPOP(data.vec, rep(1L, 4), 0)
library(testthat)
test_that("no penalty is OK", {
  best.cost <- min(fit$cost.mat[,4])
  exp.cost <- PoissonLoss(data.vec, rev(fit$mean.vec))
  expect_equal(best.cost, exp.cost)
})

data.vec <- as.integer(c(1, 10, 14, 5))
weight.vec <- rep(1L, 4)
pdpa <- PeakSegPDPA(data.vec, weight.vec, 3L)
fpop <- PeakSegFPOP(data.vec, weight.vec, 0)
test_that("FPOP computes same model as PDPA", {
  pdpa.cost <- min(pdpa$cost.mat)
  fpop.cost <- min(fpop$cost.mat)
  expect_equal(fpop.cost, pdpa.cost)
})

fit <- PeakSegFPOP(data.vec, rep(1L, 4), 1e6)
test_that("huge penalty recovers no changes", {
  exp.mean <- mean(data.vec)
  exp.mean.vec <- c(exp.mean, Inf, Inf, Inf)
  expect_equal(fit$mean.vec, exp.mean.vec)
  exp.loss <- PoissonLoss(data.vec, exp.mean)
  expect_equal(min(fit$cost.mat), exp.loss)
})

data(H3K4me3_XJ_immune_chunk1)
H3K4me3_XJ_immune_chunk1$count <- H3K4me3_XJ_immune_chunk1$coverage
by.sample <- split(
  H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)

test_that("FPOP recovers the same models as PDPA", {
  one.name <- "McGill0004"
  one <- by.sample[[one.name]]
  pdpa <- PeakSegPDPAchrom(one, 9L)
  segs.by.peaks <- split(pdpa$segments, pdpa$segments$peaks)
  some.models <- pdpa$modelSelection[-1,]
  for(model.i in 1:nrow(some.models)){
    model.row <- some.models[model.i,]
    lambda <- with(model.row, if(max.lambda==Inf){
      min.lambda+1
    }else{
      (min.lambda+max.lambda)/2
    })
    exp.segs <- segs.by.peaks[[paste(model.row$peaks)]]
    rownames(exp.segs) <- NULL
    fpop <- PeakSegFPOPchrom(one, lambda)
    expect_equal(fpop$segments, exp.segs)
  }
})
