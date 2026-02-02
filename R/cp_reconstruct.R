#' cp_reconstruct
#'
#' Computes the approximated tensor from a canonical polyadic decomposition.
#'
#' @param list_cores A list of cores generated from [cp()]
#'
#' @return The reconstructed tensor generated from the canonical polyadic decomposition.
#' @details See T. Kolda, B. Bader, "Tensor decomposition and applications".
#' SIAM Applied Mathematics and Applications 2009.
#'
#' @examples
#' A <- array(1:24, dim = c(2, 3, 4))
#' res <- cp_als(A, R = 2)
#' cp_als_reconstruct(res)
#'
#' @seealso [cp()] to create the list of cores from the higher-order SVD.
#'
#' @export
cp_reconstruct <- function(cp_als_list) {
  lambda <- cp_als_list$lambda

  mats <- cp_als_list$mats

  dims <- lapply(mats, nrow) |> unlist()
  order <- length(dims)

  A_hat <- array(0, dim = dims)

  for(r in 1:2) {
    curr <- outer(mats[[1]][, r], mats[[2]][, r])

    for(k in 3:order) {
      curr <- outer(curr, mats[[k]][, r])
    }

    A_hat <- A_hat + lambda[r] * curr
  }
  A_hat
}
