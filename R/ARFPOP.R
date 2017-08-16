ARFPOP <- structure(function(dat, gam, lambda, constraint = FALSE) {
  stopifnot(gam > 0 && gam <= 1)
  stopifnot(!is.null(gam))
  weight.vec <- rep(1, length(dat))
  n.data <- length(dat)
  stopifnot(3 <= n.data)
  stopifnot(is.numeric(lambda))
  stopifnot(length(lambda) == 1)
  stopifnot(0 <= lambda)
  cost.mat <- double(n.data)
  ends.vec <- integer(n.data)
  mean.vec <- double(n.data)
  intervals.mat <- integer(n.data)
  success <- 1
  result.list <- .C(
    "ARFPOP_interface",
    dat.vec = as.numeric(dat),
    n.data = as.integer(n.data),
    penalty = as.numeric(lambda),
    gam = as.numeric(gam),
    cost.mat = as.double(cost.mat),
    ends.vec = as.integer(ends.vec),
    mean.vec = as.double(mean.vec),
    intervals.mat = as.integer(intervals.mat),
    constraint = constraint,
    success = as.integer(success),
    PACKAGE = "FastLZeroInference"
  )
  
  ## 1-indexed segment ends!
  result.list$ends.vec <-
    result.list$ends.vec + 1L
  result.list$cost.mat <- matrix(result.list$cost.mat, 1, n.data, byrow = TRUE)
  result.list$intervals.mat <- 
    matrix(result.list$intervals.mat, 1, n.data, byrow = TRUE)
  
  changePts <- rev(unique(result.list$ends.vec))
  spikes <- changePts[-1] + 1
  
  constraint_str <- ""
  if (constraint) {
    constraint_str <- "-pos-constrained"  
  }
  
  if (result.list$success) {
    out <-
      list(
        spikes = spikes,
        fittedValues = rev(result.list$mean.vec),
        dat = dat,
        type = paste0("ar1-fpop", constraint_str), 
        changePts = changePts,
        call = match.call(),
        gam = gam,
        lambda = lambda,
        cost = as.numeric(result.list$cost.mat),
        nIntervals = as.numeric(result.list$intervals.mat)
      )
    class(out) <- "estimatedSpikes"
    return(out)
  } else {
    stop("Numerical issues!")
  }

  
})