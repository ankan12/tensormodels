#' tensor
#'
#' Creates a tensor object.
#'
#' @param data A vector containing data for the tensor.
#' @param dims A vector of dimensions for the resulting tensor.
#'
#' @return A tensor with values from data and dimensions from dims.
#'
#' @examples
#' (A <- tensor(1:24, dim = c(2, 3, 4)))
#'
#' @export
tensor <- function(data, dim = dim(data)) {
  structure(
    array(data = data, dim = dim),
    class = "tensor"
  )
}
