#' rtnorm
#'
#' Simulate random draws from the tensor variate normal distribution.
#' \deqn{\frac{(x-\mu)^2}{2\sigma^2}}
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
rtnorm <- function(n, mu = 0, sigmas = 1) {
  dims <- dim(mu)

  # mu was a scalar
  if(is.null(dims)) dims <- 1

  # no sigmas provided, use identity
  if (length(sigmas) == 1 && sigmas == 1) {
    sigmas <- lapply(dims, diag)
  }

  # get all Z draws
  X <- array(rnorm(prod(c(n, dims))), c(n, dims))

  # compute Cholesky
  chol_sigmas <- lapply(sigmas, chol)

  for (k in seq_along(chol_sigmas)) {
    X <- n_prod(X, t(chol_sigmas[[k]]), k+1) # multiply
  }
  X + array(rep(mu, each = n), dim = c(n, dims)) # add mean
}
