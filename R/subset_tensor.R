#' Subset a tensor.
#'
#' @param A A tensor object (S3 wrapper around an array).
#' @param ... Indices.
#' @param drop Logical; passed to base `[` (defaults to FALSE).
#'
#' @return A tensor that was subset.
#'
#' @examples
#' A <- tensor(array(1:24, dim = c(2, 3, 4)))
#' A[1, ]
#' @export
#' @method [ tensor
`[.tensor` <- function(x, ..., drop = TRUE) {
  dims <- dim(x) #dimensions of input

  args <- as.list(substitute(list(...)))[-1] #args provided

  # fill missing with TRUE
  if (length(args) < length(dims)) {
    args <- c(args, rep(list(TRUE), length(dims) - length(args)))
  }

  # replace empty arguments (e.g., `A[1, ]`) with TRUE
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

  # subset the underlying array
  result <- do.call(`[`, c(list(unclass(x)), eval_indices, list(drop = drop)))

  # return as tensor class
  if (is.array(result)) class(result) <- "tensor"
  result
}
