#' plot_dist
#'
#' Plots the vector and tensor Mahalanobis distance.
#'
#' @param distances A data frame with one column containing the vector distances and
#'                  one column containing the tensor distances.
#' @return A plot of the vector and tensor Mahalanobis distances.
#'
#' @examples
#' s1 <- matrix(rnorm(4), nrow = 2) |> crossprod()
#' s2 <- matrix(rnorm(16), nrow = 4) |> crossprod()
#' matnorm_draws <- rtnorm(n, mu = matrix(1:8, nrow = 2), sigmas = list(s1, s2))
#' distances <- mahalanobis_dist(matnorm_draws)
#' plot_dist(distances)
#'
#' @seealso [mahalanobis_dist()] to compute the distances.
#'
#' @export
plot_dist <- function(distances) {
  plot(distances, xlab = "Vector Dists", ylab = "Tensor Dists")
  abline(a = 0, b = 1, col = "red", lwd = 2)
}
