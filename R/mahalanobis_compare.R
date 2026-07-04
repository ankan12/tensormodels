#' mahalanobis_compare
#'
#' Computes the vector and tensor Mahalanobis distance for data using the MLE.
#'
#' @param data A list of arrays containing the independent, identically distributed draws.
#'
#' @return A data frame with one column containing the vector distances and
#'         one column containing the tensor distances.
#'
#' @examples
#' s1 <- matrix(rnorm(4), nrow = 2) |> crossprod()
#' s2 <- matrix(rnorm(16), nrow = 4) |> crossprod()
#' matnorm_draws <- rtnorm(n = 10, mu = matrix(1:8, nrow = 2), sigmas = list(s1, s2))
#' (mahalanobis_compare(matnorm_draws))
#' @seealso [plot_dist()] to create a plot of the Malahanobis distances.
#'
#' @export
mahalanobis_compare <- function(data) {
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
      centered_multiply <- n_prod(centered_multiply, sigma_inv[[j]], j)
    }

    tensor_dist[i] <- c(centered_curr) %*% c(centered_multiply)
  }

  # transform to array of n draws
  if(o == 1 && dims == 1) {
    data_array <- simplify2array(data)
  } else {
    data_array <- simplify2array(data) |> aperm(c(o+1, 1:o))
  }

  # compute vector Mahalanobis
  vec_data <- matrix(data = data_array, nrow = n, ncol = n_star)

  mu_vec <- apply(vec_data, MARGIN = 2, FUN = mean)

  vec_dist <- mahalanobis(vec_data, center = mu_vec,
                          cov = invert_safe(cov(vec_data)), inverted = TRUE)

  # scaling factor is divide by number of elements in a draw
  tensor_dist <- tensor_dist * mean(vec_dist) / mean(tensor_dist)

  data.frame(vec = vec_dist, tensor = tensor_dist)
}
