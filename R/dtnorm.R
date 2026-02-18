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
#' one_draw <- rnorm(n = 1)
#' dtnorm(one_draw, mu = 0, sigmas = list(matrix(1)))
#' dnorm(one_draw)
#' @export
dtnorm <- function(x, mu, sigmas, log = FALSE) {
  d <- dim(x)

  num_dim <- length(d)
  n_star <- prod(d)

  if(is.null(d)) d <- 1; num_dim <- 1; n_star <- length(x)

  centered <- x - mu

  all_det <- 0
  kroneck_sigmas <- 1

  for (d in length(sigmas):1) {
    s_d <- sigmas[[d]]
    U  <- chol(s_d)

    all_det <- all_det + (n_star / (2 * nrow(s_d))) * (2 * sum(log(diag(U))))
    kroneck_sigmas <- kronecker(kroneck_sigmas, chol2inv(U))
  }

  loglik <- 0

  ve_xm <- matrix(c(centered), nrow = n_star)

  delta <- t(ve_xm) %*% kroneck_sigmas %*% ve_xm

  loglik <- loglik - (n_star)/2 * log(2 * pi) - all_det + (-1/2 * delta)

  if(log == FALSE) loglik <- exp(loglik)

  loglik
}
