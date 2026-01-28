#' cp_als
#'
#' Computes the canonical polyadic decomposition of a tensor using
#' the alternate least squares algorithm.
#'
#' @param A An array containing the tensor to be decomposed.
#' @param R A rank. Must be smaller or equal to the smallest dim of A.
#' @param max_iter A maximum number of iterations to try.
#' @param tol Tolerance for convergence criterion.
#' Ends algorithm if relative Frobenius norm between true and reconstructed tensor
#' is smaller than tol.
#'
#' @return A list containing lambda which is scaling factors and a list of matrices.
#' @details See T. Kolda, B. Bader, "Tensor decomposition and applications".
#' SIAM Applied Mathematics and Applications 2009.
#'
#' @examples
#' A <- array(1:24, dim = c(2, 3, 4))
#' (res <- cp_als(A, ranks = c(2, 2, 2)))
#'
#' @seealso [cp_als_reconstruct()] to reconstruct the tensor from the decomposed parts.
#'
#' @export
cp_als <- function(A, R = 2, max_iter = 2000, tol = 1e-3) {
  dims <- dim(A)
  order <- length(dims)

  if(any(R > dims)) stop("R must be smaller than all dims")

  list_mats <- hosvd(A, ranks = rep(R, order))$mats

  for(i in 1:max_iter) {
    for(k in 1:order) {
      dims_multiply <- setdiff(1:order, k)
      dims_backward <- rev(dims_multiply)

      V <- t(list_mats[[dims_multiply[1]]]) %*% list_mats[[dims_multiply[1]]]
      An <- list_mats[[dims_backward[1]]]

      for(k_inner in dims_multiply[2:length(dims_multiply)]) {
        V <- V * (t(list_mats[[k_inner]]) %*% list_mats[[k_inner]])
      }

      for(curr in dims_backward[-1]) {
        An <- khatri_rao(An, list_mats[[curr]])
      }

      Xn <- matricization(A, k)

      An <- Xn %*% An %*% solve(V)

      lambda <- sqrt(colSums(An^2))

      An <- sweep(An, 2, lambda, "/")

      list_mats[[k]] <- An
    }
    recon_A <- cp_als_reconstruct(list(lambda = lambda, mats = list_mats))

    conv <- frob_norm(A - recon_A) / frob_norm(A)

    if(conv < tol) break
  }

  if(i == max_iter) message("Reached max iter ", max_iter)

  list(lambda = lambda, mats = list_mats)
}

