#' Coerce an array to a tensor.
#'
#' @param x An array.
#' @param obs An optional observation dimension containing multiple draws.
#'
#' @return An object of class `tensor`.
#' @export
as_tensor <- function(x, obs = NULL) {
  if (!is.array(x)) {
    stop("as_tensor() requires an array.")
  }

  tensor(x, obs = obs)
}
