#' Plot the solution to an L0 segmentation problem
#' @param x output from running estimatedSpikes
#' @param xlims optional parameter to specify the x-axis limits
#' @param ... arguments to be passed to methods
#'
#' @seealso
#' \code{\link{estimate_spikes}},
#' \code{\link{estimate_calcium}},
#' @export
#' @import graphics
#'
plot.estimated_spikes <- function(x, xlims = NULL, ...) {
  if (sum(is.na(x$estimated_calcium))) {
    stop("Calcium concentration must be estimated before plotting. (Run estimate_calcium(fit).)")
  }
  ind <- 1:length(x$dat)
  rng <- range(c(x$dat, x$estimated_calcium))
  ylims <- rng 
  if (is.null(xlims)){
    plot(ind, x$dat, cex = 0.5, pch = 20, col = "darkgrey", ylab = "", ylim = ylims, xlab = "Index")
  } else {
    plot(ind, x$dat, cex = 0.5, pch = 20, col = "darkgrey", ylab = "", ylim = ylims, xlim = xlims, xlab = "Time")
  }
  lines(ind, x$estimated_calcium, col = "blue", lwd = 2)
  
  hh <- 0.01 * diff(ylims)
  for (spike in x$spikes)
  {
    segments(x0 = ind[spike], x1 = ind[spike], y0 = ylims[1] - hh,
             y1 = hh + ylims[1], col = "blue", lwd = 1)
  }
}

#' Plot simulated data
#' @param x output data from simulateAR1
#' @param xlims optional parameter to specify the x-axis limits
#' @param ... arguments to be passed to methods
#' @seealso
#' \code{\link{estimate_spikes}},
#' \code{\link{estimate_calcium}},
#' @return Plot with simulated fluorescence (dark grey circles), calcium concentration (dark green line) and spikes (dark green tick marks on x-axis)
#'
#' @examples
#'
#' sim <- simulate_ar1(n = 500, gam = 0.998, poisMean = 0.009, sd = 0.05, seed = 1)
#' plot(sim)
#'
#' @export
#'
plot.simdata <- function(x, xlims = NULL, ...) {
  rng <- range(x$fl)
  ylims <- c(floor(rng[1]), ceiling(rng[2]))
  if (is.null(xlims)){
    plot(x$fl, cex = 0.5, pch = 20, col = "darkgrey", ylab = "", ylim = ylims)
  } else {
    plot(x$fl, cex = 0.5, pch = 20, col = "darkgrey", ylab = "", ylim = ylims, xlim = xlims)
  }
  lines(x$conc, col = "darkgreen", lwd = 2)
  
  hh <- 0.01 * diff(ylims)
  for (spike in x$spikes)
  {
    segments(x0 = spike, x1 = spike, y0 = ylims[1] - hh,
             y1 = hh + ylims[1], col = "darkgreen", lwd = 1)
  }
}

#' Print simulated data
#' @param x simulated data
#' @param ... arguments to be passed to methods
print.simdata <- function(x, ...){
  cat("\n Call: \n")
  dput(x$call)
  cat("\n Output: \n")
  cat("Number of spikes \t", length(x$spikes), "\n")
  
  cat("\n Settings: \n")
  cat("Data length \t\t", length(x$fl), "\n")
  cat("Model type \t\t", x$type, "\n")
  cat("Gamma \t\t\t", x$gam, "\n")
  cat("Pois mean \t\t", x$poisMean, "\n")
  cat("SD \t\t\t", x$sd, "\n")
}

#' Print estimated spikes
#'
#' @param x estimated spikes
#' @param ... arguments to be passed to methods
#' @export
print.estimated_spikes <- function(x, ...)
{
  cat("\n Call: \n")
  dput(x$call)
  cat("\n Output: \n")
  cat("Number of spikes \t", length(x$spikes), "\n")
  
  cat("\n Settings: \n")
  cat("Data length \t\t", length(x$dat), "\n")
  cat("Model type \t\t", x$type, "\n")
  cat("Gamma \t\t\t", x$gam, "\n")
  cat("Lambda \t\t\t", x$lambda, "\n")
}