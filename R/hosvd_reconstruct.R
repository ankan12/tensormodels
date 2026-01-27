#' hosvd_reconstruct
#'
#' Computes the approximated tensor from a higher-order SVD.
#'
#' @param list_cores A list of cores generated from [hosvd()]
#'
#' @return The reconstructed tensor generated from the higher-order SVD.
#' @details See Lathauwer, L. and De Moor, B. (2000),
#' “A Multi-Linear Singular Value Decomposition,”
#'  Society for Industrial and Applied Mathematics, 21, 1253–1278. for more details.
#'
#' @examples
#' A <- array(1:24, dim = c(2, 3, 4))
#' res <- hosvd(A, ranks = c(2, 2, 2))
#' hosvd_reconstruct(res)
#'
#' @seealso [hosvd()] to create the list of cores from the higher-order SVD.
#'
#' @export
hosvd_reconstruct <- function(hosvd_cores) {
  G_core <- hosvd_cores$G

  mats <- hosvd_cores$mats

  dims <- length(mats)

  reconstruct <- n_prod(G_core, mats[[1]], 1)

  for(i in 2:dims) {
    reconstruct <- n_prod(reconstruct, mats[[i]], i)
  }
  reconstruct
}

