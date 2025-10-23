#' rtnorm
#'
#' Simulate random draws from the tensor variate normal distribution
#'
#' @param n An integer stating the sample size.
#' @param mu An array containing the mean values of and dims of each draw.
#' @param sigmas A list of covariance matrices. Defaults to the identity.
#'
#' @return An array containing n draws of the tensor variate normal distribution.
#'
#' @examples
#' univar_norm <- rtnorm(n = 10000, mu = -2)
#' mean(univar_norm)
#' sd(univar_norm)
#' matrix_var_norm <- matrix_var_norm <- rtnorm(n = 10000, mu = matrix(1:6, nrow = 2, ncol = 3))
#' @export
rtnorm <- function(n = 50, mu = 0, sigmas = NULL) {
  all_dims <- dim(mu)

  # mu was a scalar
  if(is.null(all_dims)) all_dims <- 1

  # no sigmas provided, use identity
  if (is.null(sigmas)) {
    sigmas <- lapply(seq_along(all_dims), function(k) diag(all_dims[k]))
  }

  # get all Z draws
  X <- array(rnorm(prod(c(n, all_dims))), c(n, all_dims))

  # compute Cholesky
  chol_sigmas <- lapply(sigmas, chol)

  for (k in seq_along(chol_sigmas)) {
    X <- n_mode_prod(X, t(chol_sigmas[[k]]), k+1) # multiply
  }
  X + array(rep(mu, each = n), dim = c(n, all_dims)) # add mean
}
