#' mahalanobis_dist
#'
#' Computes the vector and tensor Mahalanobis distance for draws using the MLE.
#'
#' @param data A list of arrays containing the independent, identically distributed draws.
#'
#' @return A data frame with one column containing the vector distances and
#'         one column containing the tensor distances.
#'
#' @examples
#' s1 <- matrix(rnorm(4), nrow = 2) |> crossprod()
#' s2 <- matrix(rnorm(16), nrow = 4) |> crossprod()
#' matnorm_draws <- rtnorm(n, mu = matrix(1:8, nrow = 2), sigmas = list(s1, s2))
#' (mahalanobis_dist(matnorm_draws))
#'
#' @seealso [plot_dist()] to create a plot of the Malahanobis distances.
#'
#' @export
mahalanobis_dist <- function(draws) {
  n <- length(draws)
  dim_draws <- dim(draws[[1]])
  num_elem <- prod(dim_draws)
  o <- length(dim_draws)

  # compute tensor Mahalanobis first
  sample_stats <- tensor_mle(draws, model = "normal") # compute MLE estimates

  sample_mean <- sample_stats$mu

  tensor_dist <- rep(0, n)

  sigma_inv <- lapply(sample_stats$sigmas, invert_safe)

  for(i in 1:n) {
    centered_curr <- draws[[i]] - sample_mean # center draws

    centered_multiply <- centered_curr

    # standardize with sigma estimates
    for(j in 1:o) {
      centered_multiply <- n_prod(centered_multiply, sigma_inv[[j]], j)
    }

    tensor_dist[i] <- c(centered_curr) %*% c(centered_multiply)
  }

  # transform to array of n draws
  draws_array <- simplify2array(draws) |>
    aperm(c(o+1, 1:o))

  # compute vector Mahalanobis
  vec_draws <-
    matrix(data = draws_array, nrow = n, ncol = num_elem)

  mu_vec <- apply(vec_draws, MARGIN = 2, FUN = mean)

  vec_dist <- mahalanobis(vec_draws, center = mu_vec,
                          cov = invert_safe(cov(vec_draws)), inverted = TRUE)

  # scaling factor is divide by number of elements in a draw
  tensor_dist <- tensor_dist * mean(vec_dist) / mean(tensor_dist)

  data.frame(vec = vec_dist, tensor = tensor_dist)
}
