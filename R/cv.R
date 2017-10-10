#' Cross-validate and optimize model parameters
#' This documentation is out-of-date!!
#' @details
#' We perform cross-validation over a one-dimensional grid of \eqn{\lambda} values.
#'  For each value of \eqn{\lambda} in this grid, we solve the corresponding optimization problem, that is, one of
#'
#' \strong{AR(1)-model:}
#' minimize_{c1,...,cT} 0.5 sum_{t=1}^T ( y_t - c_t )^2 + lambda sum_{t=2}^T 1_{c_t neq gamma c_{t-1} }
#' for the global optimum, where $y_t$ is the observed fluorescence at the tth
#' timepoint.
#'
#' \strong{AR(1) with intercept:}
#' minimize_{c1,...,cT,b1,...,bT} 0.5 sum_{t=1}^T (y_t - c_t - b_t)^2 + lambda sum_{t=2}^T 1_{c_t neq gamma c_{t-1}, b_t neq b_{t-1} }
#' where the indicator variable 1_{(A,B)} equals 1 if the event A cup B holds, and equals zero otherwise.
#'
#' on a training set using a candidate value for \eqn{\gamma}. Given the resulting set of changepoints, we solve a constrained optimization problem for \eqn{\gamma}. We then refit the optimization problem with the optimized value of \eqn{\gamma},
#' and then evaluate the mean squared error (MSE) on a hold-out set. Note that in the final output of the algorithm,
#' we take the square root of the optimal value of \eqn{\gamma} in order to address the fact that the cross-validation
#' scheme makes use of training and test sets consisting of alternately-spaced timesteps.
#'
#' If there is a tuning parameter lambdaT in the path [lambdaMin, lambdaMax] that produces a fit with
#' less than 1 spike per 10,000 timesteps the path is truncated to [lambdaMin, lambdaT] and a warning is produced.
#'
#' See Algorithm 3 of Jewell and Witten (2017) <arXiv:1703.08644>
#' @seealso
#' \strong{Estimate spikes:}
#' \code{\link{estimateSpikes}},
#' \code{\link{print.estimatedSpikes}},
#' \code{\link{plot.estimatedSpikes}}.
#'
#' \strong{Cross validation:}
#' \code{\link{cv.estimateSpikes}},
#' \code{\link{print.cvSpike}},
#' \code{\link{plot.cvSpike}}.
#'
#' \strong{Simulation:}
#' \code{\link{simulateAR1}},
#' \code{\link{plot.simdata}}.
#'
#' @examples
#' # Not run
#' # sim <- simulateAR1(n = 500, gam = 0.998, poisMean = 0.009, sd = 0.05, seed = 1)
#' # plot(sim)
#'
#' # AR(1) model
#' # outAR1 <- cv.estimateSpikes(sim$fl, type = "ar1")
#' # plot(outAR1)
#' # print(outAR1)
#' # fit <- estimateSpikes(sim$fl, gam = outAR1$optimalGam[outAR1$index1SE, 1],
#' # lambda = outAR1$lambda1SE, type = "ar1")
#' # plot(fit)
#' # print(fit)
#'
#' # AR(1) + intercept model
#' # outAR1Intercept <- cv.estimateSpikes(sim$fl, type = "intercept",
#' # lambdas = seq(0.1, 5, length.out = 10))
#' # plot(outAR1Intercept)
#' # print(outAR1Intercept)
#' # fit <- estimateSpikes(sim$fl, gam = outAR1Intercept$optimalGam[outAR1Intercept$index1SE, 1],
#' # lambda = outAR1Intercept$lambda1SE, type = "intercept")
#' # plot(fit)
#' # print(fit)

#' @param dat fluorescence trace (a vector)
#' @param type type of model, must be one of AR(1) 'ar1', or AR(1) with intercept 'intercept'
#' @param gam a scalar value for the AR(1)/AR(1) + intercept decay parameter
#' @param nLambdas number of tuning parameters to estimate the model (grid of values is automatically produced)
#' @param lambdas vector of tuning parameters to use in cross-validation
#' @param hardThreshold boolean specifying whether the calcium concentration must be non-negative (in the AR-1 problem)
#'
#' @return A list of values corresponding to the 2-fold cross-validation:
#' @return \code{cvError} the MSE for each tuning parameter
#' @return \code{cvSE} the SE for the MSE for each tuning parameter
#' @return \code{lambdas} tuning parameters
#' @return \code{optimalGam} matrix of (optimized) parameters, rows correspond to tuning parameters, columns correspond to optimized parameter
#' @return \code{lambdaMin} tuning parameter that gives the smallest MSE
#' @return \code{lambda1SE} 1 SE tuning parameter
#' @return \code{indexMin} the index corresponding to lambdaMin
#' @return \code{index1SE} the index corresponding to lambda1SE
#'
#' @export
#'
cv.ARFPOP <- function(dat, gam = NULL,
                              lambdas = NULL, nLambdas = 10) {
  type = "ar1"
  hardThreshold = TRUE
  k <- 2  ## number of folds
  n <- length(dat)
  
  nDataEven <- (n %% 2) == 0
  
  if (is.null(lambdas)) {
    lambdas <- createLambdaSequence(n, nLambdas)
  } else {
    nLambdas <- length(lambdas)
  }
  
  if (is.null(gam))
  {
    optimizeGams = TRUE
    ## Modified parameters for CV
    params <- 0.998
  } else {
    optimizeGams = FALSE
    params <- gam
  }
  
  paramsTilde <- modifyParams(params, type, "fwd")
  
  foldid <- rep(seq(1, k), n)[1:n]
  cvMSE <- matrix(0, nrow = nLambdas, ncol = k)
  
  ## store all intermediate values for model parameters
  nParams <- length(params)
  if (optimizeGams) {
    paramsOut <- list()
    for (i in 1:nParams) paramsOut[[i]] <-
        matrix(0, nrow = nLambdas, ncol = k)
  } else paramsOut <- gam
  
  for (fold in 1:k) {
    trainInd <- which(foldid != fold)
    testInd <- which(foldid == fold)
    trainDat <- dat[trainInd]
    testDat <- dat[testInd]
    nn <- length(trainInd)
    nTest <- length(testDat)
    
    for (lambdaInd in 1:nLambdas) {
      if (optimizeGams) {
        segments <- ARFPOP(trainDat, paramsTilde,
                                   lambdas[lambdaInd], FALSE)
        paramsTilde <- optimParams(paramsTilde, trainDat,
                                   segments$changePts,
                                   lambdas[lambdaInd], type,
                                   hardThreshold)
      }
      segments <- ARFPOP(trainDat, paramsTilde,
                                 lambdas[lambdaInd], FALSE)
      yhatTrain <- segments$fittedValues
      nnInd <- ceiling(nn)
      yhatTest <- 0.5 * (yhatTrain[1:(nnInd - 1)] + yhatTrain[2:nnInd])
      
      if (nDataEven) {
        if (fold == 1) {
          yTest <- testDat[2 : nTest]
        }  else {
          yTest <- testDat[1 : (nTest - 1)]
        }
      } else {
        if (fold == 1) {
          yTest <- testDat[2 : (nTest - 1)]
        } else {
          yTest <- testDat
        }
      }
      
      stopifnot(length(yTest) == length(yhatTest))
      cvMSE[lambdaInd, fold] <- mean( (yhatTest - yTest) ^ 2)
      
      if (optimizeGams) {
        for (i in 1:nParams) paramsOut[[i]][lambdaInd, fold] <-
            as.numeric(paramsTilde[i])
      }
      if (length(segments$spikes) < nn / 10000 &&
          lambdaInd < nLambdas) {
        warning("Cross validation stopped early as less than 1 spike per 10,000 timesteps estimated. Rerun with smaller lambdas.")
        lambdas <- lambdas[1:lambdaInd]
        nLambdas <- length(lambdas)
        break
      }
      
    }
  }
  
  cvErr <- rowMeans(cvMSE[1: lambdaInd, ])
  cvse <- apply(cvMSE[1: lambdaInd, ], 1, sd) / sqrt(k)
  
  indexMin <- which.min(cvErr)
  lambdaMin <- lambdas[indexMin]
  lambda1SE <- max(lambdas[cvErr <= cvErr[indexMin] + cvse[indexMin]])
  index1SE <- which(lambdas == lambda1SE)
  
  if (optimizeGams) {
    paramsOutTmp <- matrix(0, ncol = nParams, nrow = nLambdas)
    for (i in 1:nParams)
      paramsOutTmp[, i] <- as.matrix(modifyParams(rowMeans(paramsOut[[i]][1: lambdaInd, ]),
                                                  type, "bck"))
    paramsOut <- paramsOutTmp
    colnames(paramsOut) <- colnames(params)
  }
  
  out <- list(cvError = cvErr, cvSE = cvse, lambdas = lambdas,
              optimalGam = paramsOut, lambdaMin = lambdaMin,
              lambda1SE = lambda1SE, indexMin = indexMin,
              index1SE = index1SE,
              call = match.call(),
              optimized = optimizeGams,
              type = type)
  class(out) <- "cvSpike"
  return(out)
  
}

yhatMSE <- function(params, y, changePts, type, hardThreshold) {
  return(mean( (y - computeFittedValues(y, changePts, params, type, hardThreshold)) ^ 2))
}

createLambdaSequence <- function(n, nLambdas = 10) {
  return(10^(seq(-1, 1, length.out = nLambdas)))
}

optimParams <- function(params, dat, changePts, penalty, type, hardThreshold) {
  if (type %in% c("ar1", "intercept")) {
    return(optimize(f = yhatMSE, interval = c(0.9, 1), y = dat, changePts = changePts, type = type, hardThreshold)$minimum)
  }
}

modifyParams <- function(params, type, direction) {
  if (type %in% c("ar1", "intercept")) {
    if (direction == "fwd")
      return(params^2)
    if (direction == "bck")
      return(sqrt(params))
  }
}


computeFittedValues <- function(dat, changePts, params, type, hardThreshold = FALSE) {
  n <- length(dat)
  nSegments <- length(changePts)
  changePts <- c(changePts, n)
  if (type == "ar1") {
    X <- matrix(0, nrow = n, ncol = nSegments)
    for (k in 1:nSegments) {
      X[(changePts[k] + 1):changePts[k + 1], k] <- params^(0:(changePts[k + 1] - (changePts[k] +
                                                                                    1)))
    }
    fit <- lm(dat ~ X - 1)
    
    if (hardThreshold == T) {
      return(pmax(fit$fitted.values, 0))
    } else {
      
      if (sum(fit$coefficients < 0) > 0) {
        warning("Check model fit carefully. In some segments calcium concentration may not 'decay' as expected. Most observed datapoints should be positive.")
      }
      return(fit$fitted.values)
    }
    
    
  }
  
  if (type == "intercept") {
    X <- matrix(0, nrow = n, ncol = 2 * nSegments)
    ind <- c(1, 2)
    for (k in 1:nSegments) {
      X[(changePts[k] + 1):changePts[k + 1], ind] <- cbind(params^(0:(changePts[k + 1] - (changePts[k] +
                                                                                            1))), rep(1, changePts[k + 1] - changePts[k]))
      ind <- c(ind[2] + 1, ind[2] + 2)
    }
    
    fit <- lm(dat ~ X - 1)
    return(fit$fitted.values)
  }
  
}

#### grid by grid


cv.full.grid.ARFPOP <- function(dat, gammas, lambdas, constrained) {
  type = "ar1"
  k <- 2  ## number of folds
  n <- length(dat)
  nGammas <- length(gammas)
  nLambdas <- length(lambdas)
  
  nDataEven <- (n %% 2) == 0
  
  foldid <- rep(seq(1, k), n)[1:n]
  cvMSE <- numeric(k)
  cvMean <- matrix(0, nrow = nLambdas, ncol = nGammas)
  cvSE <- matrix(0, nrow = nLambdas, ncol = nGammas)
  
  for (lambdaInd in 1:nLambdas) {
    for (gammaInd in 1:nGammas) {
      paramsTilde <- modifyParams(gammas[gammaInd], type, "fwd")
      
      for (fold in 1:k) {
        trainInd <- which(foldid != fold)
        testInd <- which(foldid == fold)
        trainDat <- dat[trainInd]
        testDat <- dat[testInd]
        nn <- length(trainInd)
        
        segments <- ARFPOP(trainDat, paramsTilde, lambdas[lambdaInd], constrained)   
        yhatTrain <- segments$fittedValues
        nnInd <- ceiling(nn)
        
        ## predict test points based on training fits
        yhatTest <- 0.5 * (yhatTrain[1:(nnInd - 1)] + yhatTrain[2:nnInd])
        if (nDataEven) {
          if (fold == 1) {
            yhatTest <- c(yhatTrain[1], yhatTest)
          }  else {
            yhatTest <- c(yhatTest, yhatTrain[nn])
          }
        } else {
          if (fold == 1) {
            yhatTest <- c(yhatTrain[1], yhatTest, yhatTrain[nn])
          }  
        }
        
        cvMSE[fold] <- mean( (yhatTest - testDat) ^ 2)
      }
      
      cvMean[lambdaInd, gammaInd] <- mean(cvMSE)
      cvSE[lambdaInd, gammaInd] <- sd(cvMSE) / sqrt(k)
      
    }
  }

  ind <- which(cvMean == min(cvMean), arr.ind = T)
  lambdaMin <- lambdas[ind[1]]
  gamMin <- gammas[ind[2]]
  
  out <- list(cvError = cvMean, cvSE = cvSE, 
              lambdas = lambdas,
              gammas = gammas, 
              gammaMin = gamMin, 
              lambdaMin = lambdaMin,
              call = match.call(),
              constrained = constrained,
              type = type)
  class(out) <- "cvSpike2"
  return(out)
  
}

## cv on calcium estimation

cv.calcium <- function(dat, gam, lambdas, constrained) {
  nLambda <- length(lambdas)
  fold_thin <- 10
  cvMSE <- matrix(0, nrow = nLambda, ncol = fold_thin)
  gam_tilde <- gam ^ fold_thin
  lambda_error <- numeric(nLambda)
  for (lambda_i in 1:nLambda) {
   tryCatch({
     fit_all_dat <- 
      ARFPOP(dat = dat, gam = gam, lambda = lambdas[lambda_i], 
             constraint = constrained)  
    chgpts <- c(fit_all_dat$changePts, length(dat))
    nChgpts <- length(chgpts)
    df <- NULL
    chgPts_error <- numeric(nChgpts - 1)
    for (chg_i in 1:(nChgpts - 1)){
      lower_end <- chgpts[chg_i] + 1
      upper_end <- chgpts[chg_i + 1]
      
      segment_data <- dat[lower_end:upper_end]
      segment_length <- length(segment_data)
      
      if (segment_length > fold_thin) {
      
      fold_ind <- rep(seq(1, fold_thin), segment_length)[1:segment_length]      
      fold_error <- numeric(fold_thin)
      for (k in 1:fold_thin) {
        train_ind <- (fold_ind == k)
        test_ind <- (fold_ind != k)
        
        train_dat <- segment_data[train_ind]
        n_train <- length(train_dat)
        test_dat <- segment_data[test_ind]
        
        est_init_ca <- sum(train_dat * gam_tilde ^ (0:(n_train - 1))) / sum( gam_tilde ^ (2 * (0:(n_train - 1))))
        
        fitted_points <- est_init_ca * gam ^ ((1 - k):(segment_length - k))  
        estimated_test_points <- fitted_points[test_ind]
        
        df <- rbind(df, 
                    data.frame( error = (estimated_test_points - test_dat) ^ 2, 
                                fold = k, 
                                lambda = lambdas[lambda_i]))
      }
      
      
      }
    }
    
    means <- df %>% group_by(fold) %>% summarize(mu = mean(error)) 
    
    cvMSE[lambda_i, ] <- means$mu
   }
   
   , 
   error = function(cond) {
     
   })
  }
  cvErr <- rowMeans(cvMSE)
  cvse <- apply(cvMSE, 1, sd) / sqrt(fold_thin)
  
  indexMin <- which.min(cvErr)
  lambdaMin <- lambdas[indexMin]
  lambda1SE <- max(lambdas[cvErr <= cvErr[indexMin] + cvse[indexMin]])
  index1SE <- which(lambdas == lambda1SE)
  
  out <- list(cvError = cvErr, cvSE = cvse, lambdas = lambdas,
              lambdaMin = lambdaMin,
              lambda1SE = lambda1SE, 
              indexMin = indexMin,
              index1SE = index1SE,
              call = match.call())
  class(out) <- "cvSpike--on-calcium"
  return(out)
}
