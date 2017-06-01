# require(fpop)
# fpopWrap <- function(data.vec, penalty) {
#   fit <- Fpop(data.vec, lambda = penalty)
#   end.vec <- fit$t.est
#   change.vec <- end.vec[-length(end.vec)]
#   start.vec <- c(1, change.vec+1)
#   segs.list <- list()
#   for(seg.i in seq_along(start.vec)){
#     start <- start.vec[seg.i]
#     end <- end.vec[seg.i]
#     seg.data <- data.vec[start:end]
#     seg.mean <- mean(seg.data)
#     segs.list[[seg.i]] <- data.frame(
#       start, end,
#       mean=seg.mean,
#       seg.cost=sum((seg.data-seg.mean)^2))
#   }
#   segs <- do.call(rbind, segs.list)
#   fitted <- rep(segs$mean, segs$end - segs$start + 1)
#   N <- length(data.vec)
#   feasible = T
#   if (sum(fitted[2:N] < fitted[1:(N-1)]) > 0) {
#     feasible = F}
#     
#   return(list(fitted.values = fitted, feasible = feasible))
# }
