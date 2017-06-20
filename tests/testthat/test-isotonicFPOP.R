require(Iso)
require(fpop)
require(coseg)

tests <- function(x) {
# Check that Isotonic FPOP with penalty 0 matches 
# exact PAVA algo. from Iso package
fit <- IsotonicFPOP(x, penalty = 0)
pava.fit <- pava(x)

test_that("matching FPOP lambda = 0 and pava", {expect_equal(fit$fitted.values, pava.fit)})


# Check that Isotonic FPOP with general penalty matches 
# FPOP if FPOP solution is feasible 

# Fpop wrapper
fpopWrap <- function(data.vec, penalty) {
  fit <- Fpop(data.vec, lambda = penalty)
  end.vec <- fit$t.est
  change.vec <- end.vec[-length(end.vec)]
  start.vec <- c(1, change.vec+1)
  segs.list <- list()
  for(seg.i in seq_along(start.vec)){
    start <- start.vec[seg.i]
    end <- end.vec[seg.i]
    seg.data <- data.vec[start:end]
    seg.mean <- mean(seg.data)
    segs.list[[seg.i]] <- data.frame(
      start, end,
      mean=seg.mean,
      seg.cost=sum((seg.data-seg.mean)^2))
  }
  segs <- do.call(rbind, segs.list)
  fitted <- rep(segs$mean, segs$end - segs$start + 1)
  N <- length(data.vec)
  feasible = T
  if (sum(fitted[2:N] < fitted[1:(N-1)]) > 0) {
    feasible = F}
  
  return(list(fitted.values = fitted, feasible = feasible))
}

pens <- 10 ^ seq(0.1, 3, length.out = 10)

for (pen in pens) {
  fit <- IsotonicFPOP(x, penalty = pen)
  fpop.fit <- fpopWrap(x, penalty = pen)
  if (fpop.fit$feasible) {
    test_that(paste0("matching FPOP lambda = , ", pen, 
                      "and feasible fpop"),
              {expect_equal(fit$fitted.values, fpop.fit$fitted.values)})
  }
}
}


set.seed(1)
n <- 1000

## increasing data
x <- rnorm(n) + 0.5 * (1:n) 
tests(x)

## decreasing data
x <- rnorm(n) - 0.5 * (1:n) 
tests(x)

## normal data
x <- rnorm(n) 
tests(x)
