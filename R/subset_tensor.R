#' Subset a tensor.
#'
#' @param A A tensor object (S3 wrapper around an array).
#' @param ... Indices.
#' @param drop Logical; passed to base `[` (defaults to FALSE).
#'
#' @return A tensor that was subset.
#' @export
#' @method [ tensor
`[.tensor` <- function(x, ..., drop = TRUE) {
  d <- dim(x)
  if (is.null(d)) stop("tensor must wrap an array with dim()")

  nd  <- length(d)
  idx <- as.list(substitute(list(...)))[-1]
  n   <- length(idx)

  # If fewer indices than dims, fill remaining with "missing" (= take all)
  if (n < nd) {
    idx <- c(idx, rep(list(quote(expr = )), nd - n))
  }

  # Subset underlying array
  out <- do.call(`[`, c(list(unclass(x)), idx, list(drop = drop)))

  # NumPy-like: if result is still 3+ dims, keep tensor class; otherwise return base
  if (is.array(out) && length(dim(out)) >= 3) {
    class(out) <- class(x)
  }

  out
}
