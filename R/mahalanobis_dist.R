#' mahalanobis_dist
#'
#' Computes the Mahalanobis distances for data using the MLE.
#'
#' @param data A list of arrays containing the independent, identically distributed draws.
#'
#' @return A vector of distances.
#'
#' @examples
#' s1 <- matrix(rnorm(4), nrow = 2) |> crossprod()
#' s2 <- matrix(rnorm(16), nrow = 4) |> crossprod()
#' matnorm_draws <- rtnorm(n = 10, mu = matrix(1:8, nrow = 2), sigmas = list(s1, s2))
#' (mahalanobis_dist(matnorm_draws))
#'
#' @export
mahalanobis_dist <- function(data) {
  n <- length(data)
  dims <- dim(data[[1]])
  n_star <- prod(dims)
  o <- length(dims)

  # compute tensor Mahalanobis first
  sample_stats <- tensor_mle(data, model = "normal") # compute MLE estimates

  sample_mean <- sample_stats$mu

  tensor_dist <- rep(0, n)

  sigma_inv <- lapply(sample_stats$sigmas, invert_safe)

  for(i in 1:n) {
    centered_curr <- data[[i]] - sample_mean # center draws

    centered_multiply <- centered_curr

    # standardize with sigma estimates
    for(j in 1:o) {
      centered_multiply <- n_prod(centered_multiply, j, sigma_inv[[j]])
    }

    tensor_dist[i] <- c(centered_curr) %*% c(centered_multiply)
  }

  tensor_dist
}
