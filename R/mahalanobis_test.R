#' mahalanobis_test
#'
#' Computes the Kolmogorov-Smirnov test on a data frame of tensor and vector
#' Mahalanobis distance.
#'
#' @param distances A data frame with one column containing the vector distances
#'                  and one column containing the tensor distances.
#' @return A list containing the statistic, p-value, method, and data name.
#'
#' @examples
#' s1 <- matrix(rnorm(4), nrow = 2) |> crossprod()
#' s2 <- matrix(rnorm(16), nrow = 4) |> crossprod()
#' matnorm_draws <- rtnorm(n, mu = matrix(1:8, nrow = 2), sigmas = list(s1, s2))
#' distances <- mahalanobis_dist(matnorm_draws)
#' mahalanobis_test(distances)
#'
#' @seealso [mahalanobis_dist()] to compute the distances and
#'          [plot_dist()] to create a plot of the Malahanobis distances.
#'
#' @export
mahalanobis_test <- function(distances) {
  res <- with(distances, ks.test(vec, tensor))

  res$method <- "Asymptotic two-sample Kolmogorov-Smirnov Test
                  for Vector vs Tensor Malahanobis Distances"

  res
}
