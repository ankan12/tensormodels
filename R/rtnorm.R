#' rtnorm
#'
#' Simulates random draws from the tensor variate normal distribution.
#' The density function of the multilinear normal distribution is
#' \deqn{f(x) = (2\pi)^{-p*/2} \left(\prod_{i=1}^{k} |\Sigma_i|^{-p*/(2\pi)}\right) \exp\left\{-\frac{1}{2} (x-\mu)^{T} \Sigma_{1:k}^{-1} (x-\mu)\right\}}
#'
#' @param n An integer stating the sample size.
#' @param mu An array containing the mean values of and dims of each draw.
#' @param sigmas A list of covariance matrices. Defaults to the identity.
#'
#' @return An array containing n draws of the tensor variate normal distribution.
#' @details Let \eqn{\textbf{Z}} be an array of independent standard normal entries.
#' Let \eqn{\mathcal{Y} = \mathcal{U} + \textbf{Z} \times \bold{\Sigma}} with
#' \eqn{\mathcal{U} \in \mathbb{R}^{m_{1} \times \dots \times m_{K}}} and \eqn{\bold{\Sigma}}
#' is a list of covariance matrices. For the covariance matrices,
#' \eqn{\bold{\Sigma}_{1} \in \mathbb{R}^{m_{1} \times m_{1}}, \bold{\Sigma}_{2} \in \mathbb{R}^{m_{2} \times m_{2}}, \dots \bold{\Sigma}_{K}
#' \in \mathbb{R}^{m_{K} \times m_{K}}.} \cr \cr
#' We can also write this model as
#' \eqn{\textbf{X} \sim N(\textbf{U}, \bold{\Sigma}_{1} \otimes \bold{\Sigma}_{2} \dots \otimes \bold{\Sigma}_{K}).} \cr
#' To simulate draws from the multilinear normal, instead of multiplying by the
#' full covariance matrix, we can use the Cholesky decomposition to obtain the
#' covariance in an easier equation. \cr \cr
#' For each \eqn{A_{k}}, we compute the Cholesky to obtain \eqn{\bold{\Sigma}_{k} = R_{k}^{T} RL_{k}.}
#' We can then use those Cholesky versions of the covariance matrices to sample
#' from a multilinear normal distribution. We can rewrite
#' \eqn{\mathcal{Y} = \mathcal{U} + \textbf{Z} \times \bold{\Sigma}} as
#' \eqn{\mathcal{Y} = \mathcal{U} + \textbf{Z} \times \bold{R}_{k}.}
#'
#' @examples
#' univar_norm <- rtnorm(n = 10000, mu = -2)
#' mean(univar_norm)
#' sd(univar_norm)
#' matrix_var_norm <- rtnorm(n = 10000, mu = matrix(1:6, nrow = 2, ncol = 3))
#'
#' @seealso [tucker()] to create the list of cores from the Tucker decomposition.
#' @export
rtnorm <- function(n, mu = 0, sigmas = 1) {
  dims <- dim(mu)

  # mu was a scalar
  if(is.null(dims)) dims <- length(mu)

  # no sigmas provided, use identity
  if (length(sigmas) == 1 && nrow(sigmas[[1]]) == 1 && sigmas[[1]] == 1) {
    sigmas <- lapply(dims, diag)
  }

  d <- length(dims)

  # get all Z draws
  X <- array(rnorm(prod(c(n, dims))), c(n, dims))

  # compute Cholesky
  chol_sigmas <- lapply(sigmas, chol)

  for (k in seq_along(chol_sigmas)) {
    X <- n_prod(X, t(chol_sigmas[[k]]), k+1) # multiply
  }
  norm_array <- X + array(rep(mu, each = n), dim = c(n, dims)) # add mean

  asplit(norm_array, 1) # return as list of arrays
}
