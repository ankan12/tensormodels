#' tensor product
#'
#' Compute the tensor product of a tensor with a tensor
#'
#' @param tensor_A An array representing the tensor A.
#' @param tensor_B An array representing tensor B.
#' @param simplify A bool saying whether to simplify the output. Defaults to TRUE.
#'
#' @return An array: the result of the tensor product between A and B.
#'
#' @examples
#' A <- matrix(c(1, 2, 3, 4), nrow = 2)
#' b <- matrix(c(5, 6), nrow = 2)
#' tensor_prod(A, b)
#' @export
tensor_prod <- function(tensor_A, tensor_B, simplify = TRUE) {
  tensor_prod_cpp(tensor_A, tensor_B, simplify)
}
