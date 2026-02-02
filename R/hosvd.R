#' hosvd
#'
#' Computes the higher-order SVD decomposition of a tensor.
#'
#' @param A An array containing the tensor to be decomposed.
#' @param ranks A list of ranks used for the approximation. Must be same length as order of A.
#'
#' @return A list containing a core tensor and a list of matrices.
#' @details See Lathauwer, L. and De Moor, B. (2000),
#' “A Multi-Linear Singular Value Decomposition,”
#'  Society for Industrial and Applied Mathematics, 21, 1253–1278. for more details.
#'
#' @examples
#' A <- array(1:24, dim = c(2, 3, 4))
#' (res <- hosvd(A, ranks = c(2, 2, 2)))
#'
#' @seealso [hosvd_reconstruct()] to reconstruct the tensor from the decomposed parts.
#'
#' @export
hosvd <- function(A, ranks) {
  if(missing(ranks)) stop("Must provide a vector of ranks.")

  dims <- dim(A)
  order <- length(dims)

  if(length(ranks) != order) stop("Ranks should be the same length as the order of A.")

  list_cores <- vector(mode = "list", length = order)

  for(k in 1:order) {
    if(ranks[k] > dims[k]) stop("Ranks must be equal or smaller than
                                corresponding mode of A.")

    A_mat <- matricization(A, k)

    A_mat_svd <- svd(A_mat, nu = ranks[k])

    list_cores[[k]] <- A_mat_svd$u
  }

  G_core <- n_prod(A, t(list_cores[[1]]), 1)

  for(k_mult in 2:order) {
    G_core <- n_prod(G_core, t(list_cores[[k_mult]]), k_mult)
  }

  list(G = G_core, mats = list_cores)
}

