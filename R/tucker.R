#' tucker
#'
#' Computes the canonical polyadic decomposition of a tensor using
#' the alternate least squares algorithm.
#'
#' @param A An array containing the tensor to be decomposed.
#' @param R A list of ranks of the same size as order of A. E
#' Each rank must be smaller or equal to corresponding dim.
#' @param max_iter A maximum number of iterations to try.
#' @param tol Tolerance for convergence criterion.
#' Ends algorithm if relative Frobenius norm between true and reconstructed tensor
#' is smaller than tol.
#'
#' @return A list containing a core tensor and a list of matrices.
#' @details See T. Kolda, B. Bader, "Tensor decomposition and applications".
#' SIAM Applied Mathematics and Applications 2009.
#'
#' @examples
#' A <- array(1:24, dim = c(2, 3, 4))
#' (res <- tucker(A, ranks = c(2, 2, 2)))
#'
#' @seealso [tucker_reconstruct()] to reconstruct the tensor from the decomposed parts.
#'
#' @export
tucker <- function(A, ranks = ranks, max_iter = 2000, tol = 1e-3) {
  dims <- dim(A)
  order <- length(dims)

  if(any(ranks > dims)) stop("Rank must be smaller or equal to corresponding dim")

  list_mats <- hosvd(A, ranks = ranks)$mats

  for(i in 1:max_iter) {
    for(k in 1:order) {
      dims_multiply <- setdiff(1:order, k)

      Y <- n_prod(A, t(list_mats[[dims_multiply[1]]]), dims_multiply[1])

      for(k_inner in dims_multiply[2:length(dims_multiply)]) {
        Y <- n_prod(Y, t(list_mats[[k_inner]]), k_inner)
      }

      Y_svd <- matricization(Y, k) |> svd(nu = ranks[k])

      list_mats[[k]] <- svd(matricization(Y, k), nu = ranks[k])$u
    }

    G_core <- n_prod(A, t(list_mats[[1]]), 1)

    for(k_mult in 2:order) {
      G_core <- n_prod(G_core, t(list_mats[[k_mult]]), k_mult)
    }

    recon_A <- tucker_reconstruct(list(G = G_core, mats = list_mats))

    conv <- frob_norm(A - recon_A) / frob_norm(A)

    if(conv < tol) break
  }

  if(i == max_iter) message("Reached max iter ", max_iter)

  list(G = G_core, mats = list_mats)
}

