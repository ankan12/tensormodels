#' Prints a tensor.
#'
#' @param A A tensor object (S3 wrapper around an array).
#'
#' @return Prints the tensor.
#' @export
print.tensor <- function(x, ..., max_slices = getOption("tensor.max_slices", 3),
                         edge = getOption("tensor.edge", 6),
                         stats = getOption("tensor.stats", TRUE)) {

  d <- dim(x)
  nd <- length(d)
  cat("<tensor>", paste(d, collapse = " x "),
      "|", typeof(unclass(x)), "\n")

  if (stats && is.numeric(x)) {
    xx <- as.vector(unclass(x))
    cat("  range:", paste0(range(xx, na.rm = TRUE), collapse = " .. "),
        "| mean:", signif(mean(xx, na.rm = TRUE), 2),
        "| NA:", sum(is.na(xx)), "\n")
  } else {
    cat("  NA:", sum(is.na(unclass(x))), "\n")
  }

  # Helper: show matrix with edge truncation
  print_mat <- function(m) {
    nr <- nrow(m); nc <- ncol(m)
    rr <- if (nr <= 2*edge) 1:nr else c(1:edge, (nr-edge+1):nr)
    cc <- if (nc <= 2*edge) 1:nc else c(1:edge, (nc-edge+1):nc)

    mm <- m[rr, cc, drop = FALSE]
    print(mm)

    if (nr > 2*edge || nc > 2*edge) {
      cat(sprintf("  ... (%d x %d shown; full is %d x %d)\n",
                  nrow(mm), ncol(mm), nr, nc))
    }
  }

  # Order 0/1/2: let base handle but avoid giant prints by truncating matrices
  if (is.null(d) || nd == 0L) {
    print(unclass(x))
    return(invisible(x))
  }
  if (nd == 1L) {
    v <- unclass(x)
    n <- length(v)
    if (n <= 2*edge) print(v) else {
      print(c(v[seq_len(edge)], NA, v[(n-edge+1):n]))
      cat(sprintf("  ... (%d shown; full length %d)\n", 2*edge + 1, n))
    }
    return(invisible(x))
  }
  if (nd == 2L) {
    print_mat(unclass(x))
    return(invisible(x))
  }

  # Order >= 3: print a few slices
  # We'll show slices along mode 3 by default, and fix higher modes to 1.
  arr <- unclass(x)

  # choose which k (mode-3) slices to show
  kdim <- d[3]
  ks <- if (kdim <= max_slices) seq_len(kdim) else unique(c(1:min(2, kdim), kdim))
  ks <- ks[seq_along(ks) <= max_slices]  # respect max_slices

  # indices for higher modes fixed at 1
  fixed <- if (nd > 3) rep(1L, nd - 3) else integer(0)

  for (k in ks) {
    cat(sprintf("\n[, , %d", k))
    if (nd > 3) cat(paste0(", ", paste(fixed, collapse = ", ")))
    cat("]\n")

    # build index list: [ , , k, 1, 1, ... ]
    idx <- c(list(quote(expr = ), quote(expr = ), k), as.list(fixed))
    m <- do.call(`[`, c(list(arr), idx, list(drop = TRUE)))
    # ensure it's a matrix view if possible
    if (is.null(dim(m))) {
      print(m)
    } else if (length(dim(m)) == 2L) {
      print_mat(m)
    } else {
      # if still >2D (happens if mode-1 or mode-2 is 1), just print base
      print(m)
    }
  }

  if (kdim > length(ks)) {
    cat(sprintf("\n  ... %d more slices along mode 3 not shown\n", kdim - length(ks)))
  }
  if (nd > 3) {
    cat("  (higher modes fixed at 1 for preview; index explicitly to view others)\n")
  }

  invisible(x)
}
