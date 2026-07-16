#' Construct a sample of tensor-valued observations
#'
#' `tensor_samples()` stores identically shaped tensor observations in one
#' array. The selected observation axis is moved to the first mode; the
#' remaining modes describe the shape of each observation.
#'
#' @param x An array containing identically shaped tensor observations.
#' @param obs A single numeric or character location identifying the
#'   observation mode Defaults to the first mode
#'
#' @return A `tensor_samples` object backed by an array with
#'   observations stored along its first mode
#' @export
#'
#' @examples
#' input <- array(1:24, dim = c(3, 4, 2))
#' x <- tensor_samples(input, obs = 3)
#' n_draws(x)
#' draw_shape(x)
#' pull_draw(x, 1)
tensor_samples <- function(x, obs = 1L) {
  if (!is.array(x) || length(dim(x)) < 2L) {
    stop("`x` must be an array with at least two dimensions.", call. = FALSE)
  }

  obs <- vctrs::vec_as_location2(
    obs,
    n = length(dim(x)),
    names = names(dim(x)),
    arg = "obs"
  )

  x <- unclass(x)

  if (obs != 1L) {
    permutation <- c(obs, setdiff(seq_along(dim(x)), obs))
    x <- aperm(x, permutation)
  }

  structure(
    x,
    class = c("tensor_samples", "array")
  )
}


#' Print tensor samples
#'
#' @param x Tensor samples.
#' @param ... Additional arguments, currently unused.
#' @param n Number of observations to include in the preview.
#' @param entries Number of indexed entries to show per observation.
#' @param digits Number of significant digits used to print numeric values.
#'
#' @return `x`, invisibly.
#' @export
print.tensor_samples <- function(
    x,
    ...,
    n = getOption("tensormodels.print_draws", 3L),
    entries = getOption("tensormodels.print_entries", 6L),
    digits = getOption("digits", 7L)
) {
  size <- vctrs::vec_size(x)
  shape <- draw_shape(x)
  n_label <- format(size, big.mark = ",", scientific = FALSE)

  cat(
    sprintf(
      "<tensor_samples[%s]>\n%s observations of shape %s\n",
      n_label,
      n_label,
      paste(shape, collapse = " \u00d7 ")
    )
  )

  if (size == 0L) {
    return(invisible(x))
  }

  preview_n <- .tensor_samples_print_count(
    n = n,
    maximum = vctrs::vec_size(x),
    arg = "n"
  )
  preview_entries <- .tensor_samples_print_count(
    n = entries,
    maximum = prod(shape),
    arg = "entries"
  )

  preview <- .tensor_samples_preview(
    x,
    n = preview_n,
    entries = preview_entries
  )

  cat(
    sprintf(
      "\n# Preview: %d %s \u00d7 %d indexed %s\n",
      preview_n,
      if (preview_n == 1L) "observation" else "observations",
      preview_entries,
      if (preview_entries == 1L) "entry" else "entries"
    )
  )
  print(preview, quote = FALSE, digits = digits)

  remaining_observations <- size - preview_n
  remaining_entries <- prod(shape) - preview_entries

  if (remaining_observations > 0L || remaining_entries > 0L) {
    omitted <- character()

    if (remaining_observations > 0L) {
      omitted <- c(
        omitted,
        sprintf(
          "%s more %s",
          format(remaining_observations, big.mark = ",", scientific = FALSE),
          if (remaining_observations == 1L) "observation" else "observations"
        )
      )
    }

    if (remaining_entries > 0L) {
      omitted <- c(
        omitted,
        sprintf(
          "%s more %s per observation",
          format(remaining_entries, big.mark = ",", scientific = FALSE),
          if (remaining_entries == 1L) "entry" else "entries"
        )
      )
    }

    cat(sprintf("\n# ... %s\n", paste(omitted, collapse = " and ")))
  }

  cat("# Use `pull_draw(x, i)` to inspect one observation.\n")

  invisible(x)
}

.tensor_samples_print_count <- function(n, maximum, arg) {
  n <- vctrs::vec_cast(n, integer())

  if (length(n) != 1L || is.na(n) || n < 1L) {
    stop(sprintf("`%s` must be one positive integer.", arg), call. = FALSE)
  }

  min(n, maximum)
}

.tensor_samples_preview <- function(x, n, entries) {
  values <- lapply(seq_len(n), function(i) {
    as.vector(pull_draw(x, i))[seq_len(entries)]
  })
  preview <- do.call(rbind, values)

  observation_names <- dimnames(x)[[1L]]
  if (is.null(observation_names)) {
    observation_names <- paste("obs", seq_len(n))
  } else {
    observation_names <- observation_names[seq_len(n)]
  }

  rownames(preview) <- observation_names
  colnames(preview) <- .tensor_samples_index_labels(draw_shape(x), entries)
  preview
}

.tensor_samples_index_labels <- function(shape, entries) {
  indices <- arrayInd(seq_len(entries), .dim = shape)

  apply(indices, 1L, function(index) {
    paste0("[", paste(index, collapse = ","), "]")
  })
}

#' Number of tensor observations
#'
#' @param x Tensor samples.
#'
#' @return A single integer giving the number of observations.
#' @export
n_draws <- function(x) {
  if (!inherits(x, "tensor_samples")) {
    stop("`x` must be a `tensor_samples` object.", call. = FALSE)
  }

  vctrs::vec_size(x)
}

#' Shape of each tensor observation
#'
#' @param x Tensor samples.
#'
#' @return An integer vector containing the dimensions of one observation.
#' @export
draw_shape <- function(x) {
  if (!inherits(x, "tensor_samples")) {
    stop("`x` must be a `tensor_samples` object.", call. = FALSE)
  }

  dim(x)[-1L]
}

#' Extract one tensor draw
#'
#' @param x Tensor samples.
#' @param i A single numeric or character location identifying a draw.
#'
#' @return A `tensor` containing one tensor observation.
#' @export
pull_draw <- function(x, i) {
  if (!inherits(x, "tensor_samples")) {
    stop("`x` must be a `tensor_samples` object.", call. = FALSE)
  }

  i <- vctrs::vec_as_location2(
    i,
    n = vctrs::vec_size(x),
    names = dimnames(x)[[1L]],
    arg = "i"
  )

  out <- unclass(vctrs::vec_slice(x, i))
  out_dimnames <- dimnames(out)

  if (!is.null(out_dimnames)) {
    out_dimnames <- out_dimnames[-1L]
    dimnames(out) <- NULL
  }

  dim(out) <- dim(out)[-1L]

  if (!is.null(out_dimnames)) {
    dimnames(out) <- out_dimnames
  }

  tensor(out)
}

#' Select tensor draws
#'
#' @param x Tensor samples.
#' @param i Locations identifying draws to retain.
#'
#' @return A `tensor_samples` object containing the selected observations.
#' @export
slice_draws <- function(x, i) {
  if (!inherits(x, "tensor_samples")) {
    stop("`x` must be a `tensor_samples` object.", call. = FALSE)
  }

  vctrs::vec_slice(x, i)
}
