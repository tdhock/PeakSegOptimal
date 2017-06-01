# library(Iso)
# library(fpop)
# 
# n <- 100
# x <- rnorm(n) + 0.5 * (1:n) + 1
# fit <- IsotonicFPOP(x, penalty = 0)
# yFit <- pava(x)
# 
# 
# pens <- c(0, 100, 1000)
# 
# par(mfrow = c(3, 1))
# for (pen in pens) {
#   plot(x, cex = 0.5, pch = 20, main = paste0("Lambda = ", pen))
#   lines(yFit, lwd = 2, col = "red")
#   fit <- IsotonicFPOP(x, penalty = pen)
#   lines(fit$fitted.values, lwd = 2, col = "blue")
#   
#   fitFpop <- fpopWrap(x, penalty = pen)
#   if (fitFpop$feasible)
#     lines(fitFpop$fitted.values, lwd = 2, col = "purple")
#   
# }
# 
# 
# 
# fitFpop <- fpopWrap(x, penalty = 10)
# plot(x, cex = 0.5, pch = 20)
# lines(fitFpop$fitted.values, col = "purple", lwd = 2)
# 
