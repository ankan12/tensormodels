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
  dims <- dim(x)
  o <- length(dims)
  n_star <- prod(dims)

  .validate_same_dims(mu, dims, "mu")
  sigmas <- .prepare_sigmas(sigmas, dims)

  x <- array(x, dim = dims)
  mu <- array(mu, dim = dims)

  all_det <- 0

  xm_tmp <- x - mu

  ve_xm <- matrix(c(x - mu), nrow = n_star)

  for (k in length(sigmas):1) {
    sigd <- sigmas[[k]]

    curr_inv <- invert_safe(sigd)

    xm_tmp <- n_prod(xm_tmp, curr_inv, k)

    det_sigd <- det(sigd)

    if(det_sigd == 0) det_sigd <-1e-8

    log_det <- log(det_sigd)

    all_det <- all_det + (n_star / (2 * nrow(sigd))) * (2 * log_det)
  }

  xm_tmp <- as.numeric(xm_tmp)

  delta <- sum(xm_tmp * ve_xm)

  loglik <- - (n_star)/2 * log(2 * pi) - all_det + (-1/2 * delta)

  if(log == FALSE) loglik <- exp(loglik)

  loglik
}
