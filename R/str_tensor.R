#' Structure of a tensor.
#'
#' @param A A tensor object (S3 wrapper around an array).
#'
#' @return The structure of a tensor.
#' @export
str.tensor <- function(object, ...) {
  d <- dim(object)

  paste0(sprintf("%d order tensor ", length(d)),
         paste(d, collapse = " x "),
         sprintf(" <%s>", typeof(unclass(object))))
}
