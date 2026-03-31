#' dtnorm
#'
#' Computes the density function for the tensor variate normal.
#'
#' @param x An array of quantiles.
#' @param mu An array of the mean.
#' @param sigmas A list of covariance matrices.
#' @param log Defaults to FALSE. If TRUE, return the log of the density.
#'
#' @return The log likelihood.
#'
#' @examples
#' dtnorm(x = array(1), mu = 0, sigmas = list(matrix(1)))
#' dnorm(1)
#' dtnorm(x = array(1:12, dim = c(2, 3)), mu = array(0, dim = c(2, 3)),
#'        sigmas = lapply(c(2, 3), diag))
#' @export
dtnorm <- function(x, mu = NULL, sigmas = NULL, log = FALSE) {
  d <- dim(x)
  o <- length(d)
  n_star <- prod(d)

  if(is.null(d)) {
    d <- 1
    o <- 1
    n_star <- length(x)
  }

  all_det <- 0

  xm_tmp <- x - mu

  ve_xm <- matrix(c(x - mu), nrow = n_star)

  for (d in length(sigmas):1) {
    sigd <- sigmas[[d]]

    curr_inv <- invert_safe(sigd)

    xm_tmp <- n_prod(xm_tmp, curr_inv, d)

    all_det <- all_det + (n_star / (2 * nrow(sigd))) * (2 * sum(log(det(sigd))))
  }

  xm_tmp <- as.numeric(xm_tmp)

  delta <- sum(xm_tmp * ve_xm)

  loglik <- - (n_star)/2 * log(2 * pi) - all_det + (-1/2 * delta)

  if(log == FALSE) loglik <- exp(loglik)

  loglik
}
