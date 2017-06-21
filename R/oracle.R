### Compute Oracle model complexity from paper of Cleynen et al.
oracleModelComplexity <- function(bases, segments){
  stopifnot(is.numeric(bases))
  stopifnot(is.numeric(segments))
  under.sqrt <- 1.1 + log(bases/segments)
  in.square <- 1 + 4 * sqrt(under.sqrt)
  segments * in.square * in.square
### numeric vector of model complexity values.
}
 
