## implement crops algo 
## https://www.tandfonline.com/doi/abs/10.1080/10618600.2015.1116445

# CROPS algorithm
# input: A dataset y1:n = (y1,y2,...,yn);
# Minimum and maximum values of the penalty,
# βmin and βmax;
# CPD, an algorithm such as PELT, for solving
# the penalized optimization problem.
# output: The details of optimal segmentations 


## get arfpop stats
## input: arfpop object 
## output: dataframe with cols (lambda, m(lambda), Q_m(lambda)(y1_n))
arfpop_stats <- function(fit) {
  df <- data.frame(lambda = fit$lambda,
                   changepoints_n = length(fit$spikes), 
                   cost = fit$cost[length(fit$cost)] - length(fit$spikes) * fit$lambda)
  return(df)
}

update_path_stats <- function(df, new_fit) {
  return(rbind(df, arfpop_stats(new_fit)))
}

get_num_changepts <- function(lambda, df) {
  df$changepoints_n[df$lambda == lambda][1]
}

get_cost <- function(lambda, df) {
  df$cost[df$lambda == lambda][1]
}

## use ARFPOP here
crops <- function(y, gam, lambda_min, lambda_max, constraint, EPS, sensitivity = 1) {
  lambdas_used <- c(lambda_min, lambda_max)
  path_fits <- list()
  path_stats <- NULL
  
  # 1. Run CPD for penalty values βmin and βmax; 
  path_fits[[1]] <- ARFPOP(y, gam, lambda_min, constraint, EPS)
  path_stats <- update_path_stats(path_stats, path_fits[[1]])
  path_fits[[2]] <- ARFPOP(y, gam, lambda_max, constraint, EPS)
  path_stats <- update_path_stats(path_stats, path_fits[[2]])
  
  n_fits <- 2
  
  # 2. Set β∗ = {[βmin, βmax]};
  lambda_star <- list(lambdas_used)
  
  while (length(lambda_star) > 0) { 
    # 3. Choose an element of β∗; denote this element as [β0,β1];
    # here always take the first element of list
    current_interal <- lambda_star[[1]]
    print(paste0("testing interval ", current_interal))
    print(paste0("number of chg pts interval 1 ", get_num_changepts(current_interal[1], path_stats), 
                 "number of chg pts interval 2 ", get_num_changepts(current_interal[2], path_stats)))
    if (get_num_changepts(current_interal[1], path_stats) > 
        get_num_changepts(current_interal[2], path_stats) + sensitivity) {
      # print(paste("cost diff", get_cost(current_interal[2], path_stats) - get_cost(current_interal[1], path_stats)))
      lambda_int <- (get_cost(current_interal[2], path_stats) - 
                       get_cost(current_interal[1], path_stats)) / (
                         get_num_changepts(current_interal[1], path_stats) - 
                           get_num_changepts(current_interal[2], path_stats)
                       )
      
      n_fits <- n_fits + 1
      # print(paste('interval ', current_interal))
      # print(paste('-------- Testing tuning value: ', lambda_int, ' --------'))
      path_fits[[n_fits]] <- ARFPOP(y, gam, lambda_int, constraint, EPS)
      path_stats <- update_path_stats(path_stats, path_fits[[n_fits]])
      
      if (get_num_changepts(lambda_int, path_stats) != 
          get_num_changepts(current_interal[1], path_stats)) {
        # Set β∗ = {β∗,[β0,βint),[βint,β1]}.;  
        n_intervals <- length(lambda_star)
        lambda_star[[n_intervals + 1]] <- c(current_interal[1], lambda_int)
        lambda_star[[n_intervals + 2]] <- c(lambda_int, current_interal[2])
      }
    }  
    lambda_star[[1]] <- NULL
  }
  return(list(path_stats = path_stats[sort.int(path_stats$lambda, index.return = T)$ix, ], 
              path_fits = path_fits))
}
