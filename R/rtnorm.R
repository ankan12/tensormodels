#' rtnorm
#'
#' Simulate random draws from the tensor variate normal distribution
#'
#' @param n An integer stating the sample size.
#' @param mu An array containing the mean values of and dims of each draw.
#' @param listSigmas A list of covariance matrices. Defaults to the identity.
#'
#' @return An array containing n draws of the tensor variate normal.
#'
#' @examples
#' univarNorm <- rtnorm(n = 10000, mu = -2)
#' mean(univarNorm)
#' sd(univarNorm)
#' matrixVarNorm <- matrixVarNorm <- rtnorm(n = 10000, mu = matrix(1:6, nrow = 2, ncol = 3),
#' listSigmas = list(diag(2), diag(3)))
#' @export
rtnorm <- function(n = 50, mu = 0, list_sigmas = NULL) {
  allDims <- dim(mu)

  if(is.null(allDims)) { #mu was a scalar
    allDims <- c(1)
  }

  if (is.null(list_sigmas)) { #no sigmas provided, use identity
    list_sigmas <- lapply(seq_along(allDims), function(k) diag(allDims[k]))
  }

  X <- array(rnorm(prod(c(n, allDims))), c(n, allDims)) #get all Z draws
  cholSigma <- lapply(list_sigmas, chol) #compute Cholesly

  for (k in seq_along(cholSigma)) {
    X <- n_mode_prod(X, t(cholSigma[[k]]), k+1) #multiply
  }
  X + array(rep(mu, each = n), dim = c(n, allDims)) #add mean
}
