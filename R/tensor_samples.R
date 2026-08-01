#' Print tensor-valued observations
#' @param x A tensor object.
#' @param ... Unused.
#' @param n Number of draws to preview.
#' @param entries Number of entries per draw to preview.
#' @param digits Significant digits for the preview.
#' @return `x`, invisibly.
#' @export
print.tensor <- function(x, ..., n = getOption("tensortools.print_draws", 3L),
                         entries = getOption("tensortools.print_entries", 6L),
                         digits = getOption("digits", 7L)) {
  size <- n_draws(x); shape <- draw_shape(x)
  label <- format(size, big.mark = ",", scientific = FALSE)
  cat(sprintf("<tensor[%s]>\n%s %s of shape %s\n", label, label,
              if (size == 1L) "draw" else "draws", paste(shape, collapse = " × ")))
  if (size == 0L) return(invisible(x))
  n <- .tensor_print_count(n, size, "n")
  entries <- .tensor_print_count(entries, prod(shape), "entries")
  vals <- lapply(seq_len(n), function(i) {
    as.vector(.tensor_single_draw_array(pull_draw(x, i)))[seq_len(entries)]
  })
  preview <- do.call(rbind, vals)
  dn <- dimnames(x)[[1L]]; if (is.null(dn)) dn <- paste("draw", seq_len(n)) else dn <- dn[seq_len(n)]
  rownames(preview) <- dn
  colnames(preview) <- .tensor_index_labels(shape, entries)
  cat(sprintf("\n# Preview: %d %s × %d indexed %s\n", n,
              if (n == 1L) "draw" else "draws", entries,
              if (entries == 1L) "entry" else "entries"))
  print(preview, quote = FALSE, digits = digits)
  if (size > n || prod(shape) > entries) {
    omitted <- character()
    if (size > n) omitted <- c(omitted, sprintf("%d more draws", size - n))
    if (prod(shape) > entries) omitted <- c(omitted, sprintf("%d more entries per draw", prod(shape) - entries))
    cat(sprintf("\n# ... %s\n", paste(omitted, collapse = " and ")))
  }
  invisible(x)
}

.tensor_print_count <- function(n, maximum, arg) {
  n <- vctrs::vec_cast(n, integer())
  if (length(n) != 1L || is.na(n) || n < 1L) stop(sprintf("`%s` must be one positive integer.", arg), call. = FALSE)
  min(n, maximum)
}
.tensor_index_labels <- function(shape, entries) {
  apply(arrayInd(seq_len(entries), .dim = shape), 1L,
        function(i) paste0("[", paste(i, collapse = ","), "]"))
}

#' Number of tensor draws
#' @param x A tensor object.
#' @return Number of draws.
#' @export
n_draws <- function(x) {
  if (!inherits(x, "tensor")) stop("`x` must be a `tensor` object.", call. = FALSE)
  dim(unclass(x))[1L]
}

#' Tensor draw dimensions
#'
#' `dim()` reports the shape of one draw. Use [n_draws()] for the number of
#' IID observations; the leading observation dimension remains in the internal
#' array returned by `unclass()`.
#'
#' @param x A tensor object.
#' @return Dimensions of one tensor draw.
#' @export
dim.tensor <- function(x) draw_shape(x)

#' Shape of each tensor draw
#' @param x A tensor object.
#' @return Dimensions of one draw.
#' @export
draw_shape <- function(x) {
  if (!inherits(x, "tensor")) stop("`x` must be a `tensor` object.", call. = FALSE)
  dim(unclass(x))[-1L]
}

#' Extract one tensor draw
#' @param x A tensor object.
#' @param i Draw location.
#' @return A one-draw tensor.
#' @export
pull_draw <- function(x, i) {
  if (!inherits(x, "tensor")) stop("`x` must be a `tensor` object.", call. = FALSE)
  i <- vctrs::vec_as_location2(i, n = n_draws(x), names = dimnames(x)[[1L]], arg = "i")
  .new_tensor_array(.tensor_slice_array(x, i))
}

#' Select tensor draws
#' @param x A tensor object.
#' @param i Draw locations.
#' @return A tensor containing selected draws.
#' @export
slice_draws <- function(x, i) {
  if (!inherits(x, "tensor")) stop("`x` must be a `tensor` object.", call. = FALSE)
  i <- vctrs::vec_as_location(i, n = n_draws(x), names = dimnames(x)[[1L]], arg = "i")
  .new_tensor_array(.tensor_slice_array(x, i))
}
