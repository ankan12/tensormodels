#' tt_reconstruct
#'
#' Computes the approximated tensor from a tensor-train decomposition.
#'
#' @param list_cores A list of cores generated from [tt()]
#'
#' @return The reconstructed tensor generated from the tensor-train decomposition.
#' @details See I. V. Oseledets, (2011). Tensor-Train Decomposition.
#' SIAM J. SCI. COMPUT. for more details.
#'
#' @examples
#' A <- array(1:24, dim = c(2, 3, 4))
#' res <- tt_svd(A, ranks = c(2, 2, 2))
#' tt_svd_reconstruct(res)
#'
#' @seealso [tt()] to create the list of cores from the tensor-train decomposition.
#'
#' @export
tt_reconstruct <- function(list_cores) {
  order <- length(list_cores)

  reconstruct <- nm_prod(list_cores[[1]], list_cores[[2]], 2, 1)

  for(i in 3:order) {
    reconstruct <- nm_prod(reconstruct, list_cores[[i]], i, 1)
  }
  reconstruct
}
