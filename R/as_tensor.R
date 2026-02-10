#' Coerce an array to a tensor.
#'
#' @param x An array.
#'
#' @return An object of class `tensor`.
#' @export
as_tensor <- function(x) {
  if (!is.array(x)) {
    stop("as_tensor() requires an array.")
  }

  class(x) <- c("tensor", class(x))
  x
}
