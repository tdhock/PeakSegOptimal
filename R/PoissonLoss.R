### Compute the weighted Poisson loss function, which is seg.mean -
### count * log(seg.mean). The edge case is when the mean is zero, in
### which case the probability mass function takes a value of 1 when
### the data is 0 (and 0 otherwise). Thus the log-likelihood of a
### maximum likelihood segment with mean zero must be zero.
PoissonLoss <- structure(function(count, seg.mean, weight=1){
  stopifnot(is.numeric(count))
  stopifnot(is.numeric(seg.mean))
  stopifnot(is.numeric(weight))
  n.data <- length(count)
  if(length(seg.mean) == 1){
    seg.mean <- rep(seg.mean, n.data)
  }
  if(length(weight) == 1){
    weight <- rep(weight, n.data)
  }
  stopifnot(n.data == length(seg.mean))
  stopifnot(n.data == length(weight))
  if(any(weight < 0)){
    stop("PoissonLoss undefined for negative weight")
  }
  if(any(seg.mean < 0)){
    stop("PoissonLoss undefined for negative segment mean")
  }
  not.integer <- round(count) != count
  not.positive <- count < 0
  loss <-
    ifelse(not.integer|not.positive, Inf,
           ifelse(seg.mean == 0,
                  ifelse(count == 0, 0, Inf),
                  seg.mean - count * log(seg.mean)))
  sum(loss*weight)
}, ex=function(){
  PoissonLoss(1, 1)
  PoissonLoss(0, 0)
  PoissonLoss(1, 0)
  PoissonLoss(0, 1)
})

