## quick R wrapper
fit_changepoint_model <- function(y, chgpts, gam, TOL) { 
  n <- length(y)
  chgpts <- c(chgpts, n)
  k <- length(chgpts)
  fitted_values_full <- NULL
  for (i in 1:(k - 1)) {
    y_subset <- y[(chgpts[i] + 1):chgpts[i + 1]]
    nn <- length(y_subset)
    rss_best <- Inf
    for (start_i in 1:nn){
      y_subset_trim <- y_subset[1:start_i]
      m <- length(y_subset_trim)
      initial_value <- regression_coef(y_subset_trim, gam)
      fitted_values <- pmax(fit_from_regression(initial_value, m, gam), TOL)
      fitted_values <- c(fitted_values, rep(TOL, nn - start_i))
      rss_now <- 0.5 * sum((fitted_values - y_subset) ^ 2)
      if (rss_now < rss_best) fitted_values_best <- fitted_values
    }
    fitted_values_full <- c(fitted_values_full, fitted_values_best)
  }
  return(fitted_values_full)
}


regression_coef <- function(y, gam) { 
  n <- length(y)
  prefactor <- (1 - gam ^ 2) / (1 - gam ^ (2 * n))
  ss <- 0
  for (i in 1:n) {
    ss <- ss + y[i] * gam ^ (i - 1)
  }
  if (ss < 0) ss <- 0
  return(ss * prefactor)
}

fit_from_regression <- function(initial_value, n, gam) {
  fitted_values <- initial_value * rep(1, n)
  fitted_values <- fitted_values * gam ^ (0:(n - 1))
  return(fitted_values)  
}
