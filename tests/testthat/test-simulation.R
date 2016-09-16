library(testthat)
context("simulation")
##R -d "valgrind --leak-check=full --show-reachable=yes" -f small_test.R > log.txt 2>&1
require(coseg)
##require(Segmentor3IsBack)

### simulate
bckg <- 1
up <- 10
true.signal <- rep(c(bckg, up, bckg, up, bckg), c(1000, 10, 10, 10, 1000))
true.change <- which(diff(true.signal)!=0)

set.seed(1)
Y <- rpois(lambda=true.signal, n=length(true.signal))

i_unique <- which(diff(c(-1, Y))!=0)
w_unique <- diff(c(i_unique, length(Y)+1))
## sum(rep(Y[i_unique], w_unique) != Y) ## sanity check
Y.w <- Y[i_unique]

test_that("no problem with weights", {
  fpop <- PeakSegFPOP(Y.w, w_unique, penalty = 1*log(length(Y)))
  pdpa <- PeakSegPDPA(Y.w, w_unique, max.segments = 5L)
})

test_that("no problem without weights", {
  fpop <- PeakSegFPOP(Y, penalty = 1*log(length(Y)))
  pdpa <- PeakSegPDPA(Y, max.segments = 5L)
})




