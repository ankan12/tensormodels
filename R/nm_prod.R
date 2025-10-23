#' nm-mode product
#'
#' Compute the nm-mode product of a tensor with a tensor
#'
#' @param tensor_A An array representing the tensor A.
#' @param tensor_B An array representing tensor B.
#' @param n An integer specifying which mode of tensor A to multiply across.
#' @param m An integer specifying which mode of tensor B to multiply across.
#'
#' @return An array: the result of the nm-mode product between A and B.
#'
#' @examples
#' A <- matrix(c(1, 2, 3, 4), nrow = 2)
#' b <- matrix(c(5, 6), nrow = 2)
#' nm_prod(A, b, 1, 1)
#' @export
nm_prod <- function(tensor_A, tensor_B, n, m) {
  mn_mode_prod_cpp(tensor_A, tensor_B, n, m)
}
