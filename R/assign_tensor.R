#' Assign a tensor.
#'
#' @param A A tensor object (S3 wrapper around an array).
#' @param ... Indices.
#' @param value Tensor that will be assigned to A.
#'
#' @return A tensor that was assigned.
#'
#' @examples
#' A <- tensor(array(1:24, dim = c(2, 3, 4)))
#' A[1, ] <- array(1:12, dim = c(3, 4))
#' @export
#' @method [ tensor
`[<-.tensor` <- function(x, ..., value) {
  dims <- dim(unclass(x))

  args <- as.list(substitute(list(...)))[-1]

  # fill missing dimensions with TRUE
  if (length(args) < length(dims)) {
    args <- c(args, rep(list(TRUE), length(dims) - length(args)))
  }

  # replace empty args with TRUE  (e.g., A[1, ])
  for (i in seq_along(args)) {
    if (length(args[[i]]) == 0) args[[i]] <- TRUE
  }

  # capture calling environment (so j in for-loop is found)
  caller <- parent.frame()

  # Evaluate indices safely
  eval_indices <- lapply(args, function(a) {
    if (is.symbol(a) && as.character(a) == "") TRUE
    else eval(a, caller)
  })

  arr <- unclass(x)
  out <- do.call(`[<-`, c(list(arr), eval_indices, list(value = value)))

  class(out) <- class(x)
  out
}
