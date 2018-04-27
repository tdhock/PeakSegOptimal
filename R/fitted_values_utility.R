## quick R wrapper
fit_changepoint_model <- function(y, chgpts, gam) { 
  n <- length(y)
  chgpts <- c(chgpts, n)
  k <- length(chgpts)
  fitted_values <- NULL
  for (i in 1:(k - 1)) {
    y_subset <- y[(chgpts[i] + 1):chgpts[i + 1]]
    nn <- length(y_subset)
    initial_value <- regression_coef(y_subset, gam)
    fitted_values <- c(fitted_values, fit_from_regression(initial_value, nn, gam))
  }
  return(fitted_values)
}

regression_coef <- function(y, gam) { 
  n <- length(y)
  prefactor <- (1 - gam ^ 2) / (1 - gam ^ (2 * n))
  ss <- 0
  for (i in 1:n) {
    ss <- ss + y[i] * gam ^ (i - 1)
  }
  return(ss * prefactor)
}

fit_from_regression <- function(initial_value, n, gam) {
  fitted_values <- initial_value * rep(1, n)
  fitted_values <- fitted_values * gam ^ (0:(n - 1))
  return(fitted_values)  
}
