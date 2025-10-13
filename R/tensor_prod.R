#' tensor product
#'
#' Compute the tensor product of a tensor with a tensor
#'
#' @param tensorA An array representing the tensor A.
#' @param tensorB An array representing tensor B.
#' @param simplify A bool saying whether to simplify the output. Defaults to TRUE.
#'
#' @return An array: the result of the tensor product between A and B.
#'
#' @examples
#' A <- matrix(c(1, 2, 3, 4), nrow = 2)
#' b <- matrix(c(5, 6), nrow = 2)
#' tensor_prod(A, b)
#' @export
tensor_prod <- function(tensorA, tensorB, simplify = TRUE) {
  tensor_prod_cpp(tensorA, tensorB, simplify)
}
