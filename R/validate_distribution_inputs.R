.tensor_dims <- function(x) {
  if (inherits(x, "tensor")) {
    if (n_draws(x) != 1L) {
      stop("Tensor-valued model parameters must contain exactly one draw.", call. = FALSE)
    }
    return(draw_shape(x))
  }

  dims <- dim(x)
  if(is.null(dims)) {
    return(length(x))
  }

  # Treat one-row/one-column matrices as one-mode multivariate objects.
  if(length(dims) == 2 && any(dims == 1)) {
    return(max(dims))
  }

  dims
}

.normalize_dims <- function(dims) {
  if(is.null(dims)) 1L else dims
}

.is_default_sigmas <- function(sigmas) {
  is.list(sigmas) &&
    length(sigmas) == 1 &&
    is.matrix(sigmas[[1]]) &&
    nrow(sigmas[[1]]) == 1 &&
    ncol(sigmas[[1]]) == 1 &&
    sigmas[[1]][1, 1] == 1
}

.validate_sigmas <- function(sigmas, dims) {
  dims <- .normalize_dims(dims)

  if(!is.list(sigmas) || length(sigmas) != length(dims)) {
    stop("sigmas must be a list with one covariance matrix for each dimension.")
  }

  for(k in seq_along(sigmas)) {
    sigk <- sigmas[[k]]

    if(!is.matrix(sigk) || nrow(sigk) != dims[k] || ncol(sigk) != dims[k]) {
      stop(sprintf(
        "sigmas[[%d]] must be a %d x %d covariance matrix.",
        k, dims[k], dims[k]
      ))
    }
  }

  sigmas
}

.prepare_sigmas <- function(sigmas, dims) {
  dims <- .normalize_dims(dims)

  if(.is_default_sigmas(sigmas)) {
    sigmas <- lapply(dims, diag)
  }

  .validate_sigmas(sigmas, dims)
}

.validate_same_dims <- function(value, dims, name, reference = "x") {
  if(!identical(.tensor_dims(value), dims)) {
    stop(sprintf("%s must have the same dimensions as %s.", name, reference))
  }
}
