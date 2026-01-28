#' khatri_rao
#'
#' Computes the higher-order SVD decomposition of a tensor.
#'
#' @param A A matrix of size \eqn{I \times K}
#' @param B A matrix of size \eqn{J \times K}
#'
#' @return A matrix of size \eqn{IJ \times K} containing the Khatri-Rao product of A and B.
#' @details Also known as the columnwise Kronecker product.
#' Computes the Kronecker of the columns of A and B.
#'
#' @examples
#' A <- matrix(1:6, nrow = 2)
#' B <- matrix(1:12, nrow = 4)
#' khatri_rao(A, B)
#'
#' @export
khatri_rao <- function(A, B) {
  I <- nrow(A)
  K <- ncol(A)
  J <- nrow(B)

  if(ncol(A) != ncol(B)) stop("Columns of A and B must be same size.")

  res <- matrix(nrow = I*J, ncol = K)

  for(k in 1:K) {
    res[, k] <- kronecker(A[, k], B[, k])
  }

  res
}
