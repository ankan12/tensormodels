#' Rosenblatt transform for tensor-variate distributions
#'
#' Applies a sequential Rosenblatt transformation to observations from any
#' distribution supported by [tensor_mle()]. The tensor normal transformation
#' uses exact modewise Cholesky innovations. The four skewed distributions use
#' deterministic adaptive quadrature over their updated generalized inverse
#' Gaussian mixing distributions.
#'
#' @param draws A `tensor` object containing one or more observations.
#' @param model One of `"normal"`, `"skewt"`, `"vargamma"`, `"invgauss"`,
#'   or `"genhyper"`.
#' @param parameters A named list of model parameters. The output of
#'   [tensor_mle()] can be passed directly. Every model requires `mu` and
#'   `sigmas`; skewed models also require `skew`. In addition, `"skewt"`
#'   requires `nu`, `"vargamma"` requires `gamma` (or `scale`),
#'   `"invgauss"` requires `kappa`, and `"genhyper"` requires `lambda` and
#'   `omega`.
#' @param rel.tol Relative tolerance passed to [stats::integrate()] for skewed
#'   models.
#' @param abs.tol Absolute tolerance passed to [stats::integrate()] for skewed
#'   models.
#' @param subdivisions Maximum number of subintervals used by
#'   [stats::integrate()] for skewed models.
#' @param show_progress Whether to report progress after each tensor draw.
#'
#' @return A `tensor` object with the same draw count and draw shape as
#'   `draws`. With known parameters under the selected model, its entries are
#'   independent Uniform(0, 1) variables in R's column-major tensor order.
#'
#' @details
#' For a skewed model, modewise whitening gives
#' \deqn{Y = W a + \sqrt{W} Z,}
#' where the entries of \eqn{Z} are independent standard normals and the
#' initial mixing distribution is \eqn{GIG(\lambda_0,\chi_0,\psi_0)}.
#' After conditioning on the first \eqn{r} vectorized entries,
#' \deqn{\lambda_r = \lambda_0-r/2,\quad
#'       \chi_r = \chi_0+\sum_{j=1}^r y_j^2,\quad
#'       \psi_r = \psi_0+\sum_{j=1}^r a_j^2.}
#' The next conditional CDF is the expectation of
#' \eqn{\Phi\{(y_{r+1}-W a_{r+1})/\sqrt{W}\}} under this updated GIG
#' distribution.
#'
#' The quadrature is evaluated on a log-\eqn{W} scale centered near the mode
#' of the updated GIG distribution. This is more stable than integrating
#' directly over \eqn{(0,\infty)} when the sequential posterior is highly
#' concentrated.
#'
#' If `parameters` were estimated from `draws`, the returned values are not
#' exactly independent uniforms. A formal goodness-of-fit p-value therefore
#' requires a fitted-parameter calibration.
#'
#' @examples
#' set.seed(1)
#' mu <- array(0, dim = c(2, 2))
#' sigmas <- list(diag(2), diag(2))
#' draws <- rtnorm(10, mu = mu, sigmas = sigmas)
#'
#' uniforms <- tensor_rosenblatt(
#'   draws,
#'   model = "normal",
#'   parameters = list(mu = mu, sigmas = sigmas)
#' )
#' uniforms
#' @export
tensor_rosenblatt <- function(
    draws,
    model = c("normal", "skewt", "vargamma", "invgauss", "genhyper"),
    parameters,
    rel.tol = 1e-7,
    abs.tol = 1e-9,
    subdivisions = 100L,
    show_progress = FALSE) {
  model <- match.arg(model)

  if (!inherits(draws, "tensor")) {
    stop("`draws` must be a `tensor` object.", call. = FALSE)
  }
  if (n_draws(draws) < 1L) {
    stop("`draws` must contain at least one observation.", call. = FALSE)
  }
  if (!is.list(parameters)) {
    stop("`parameters` must be a named list.", call. = FALSE)
  }

  mu <- .tensor_rosenblatt_parameter(parameters, "mu", model)
  sigmas <- .tensor_rosenblatt_parameter(parameters, "sigmas", model)
  dims <- draw_shape(draws)

  .validate_same_dims(mu, dims, "parameters$mu", "draws")
  sigmas <- .prepare_sigmas(sigmas, dims)
  .tensor_rosenblatt_validate_controls(
    rel.tol,
    abs.tol,
    subdivisions,
    show_progress
  )

  subdivisions <- as.integer(subdivisions)
  mu <- array(as.numeric(mu), dim = dims)
  whitening <- .tensor_rosenblatt_whitening(sigmas)

  if (model == "normal") {
    return(.tensor_rosenblatt_normal(
      draws = draws,
      mu = mu,
      whitening = whitening,
      show_progress = show_progress
    ))
  }

  skew <- .tensor_rosenblatt_parameter(parameters, "skew", model)
  .validate_same_dims(skew, dims, "parameters$skew", "draws")
  skew <- array(as.numeric(skew), dim = dims)

  initial_gig <- .tensor_rosenblatt_initial_gig(model, parameters)

  .tensor_rosenblatt_mixture(
    draws = draws,
    mu = mu,
    skew = skew,
    whitening = whitening,
    lambda = initial_gig$lambda,
    chi = initial_gig$chi,
    psi = initial_gig$psi,
    rel.tol = rel.tol,
    abs.tol = abs.tol,
    subdivisions = subdivisions,
    show_progress = show_progress
  )
}

.tensor_rosenblatt_parameter <- function(parameters, name, model) {
  value <- parameters[[name]]
  if (is.null(value)) {
    stop(
      "`parameters` must contain `", name, "` for model \"",
      model, "\".",
      call. = FALSE
    )
  }
  value
}

.tensor_rosenblatt_positive_parameter <- function(value, name, model) {
  if (length(value) != 1L || !is.numeric(value) || is.na(value) ||
      !is.finite(value) || value <= 0) {
    stop(
      "`parameters$", name, "` must be one positive finite number for ",
      "model \"", model, "\".",
      call. = FALSE
    )
  }
  as.numeric(value)
}

.tensor_rosenblatt_initial_gig <- function(model, parameters) {
  if (model == "skewt") {
    nu <- .tensor_rosenblatt_positive_parameter(
      .tensor_rosenblatt_parameter(parameters, "nu", model),
      "nu",
      model
    )
    return(list(lambda = -nu / 2, chi = nu, psi = 0))
  }

  if (model == "vargamma") {
    gamma <- parameters$gamma
    scale <- parameters$scale

    if (!is.null(gamma) && !is.null(scale) &&
        !isTRUE(all.equal(gamma, scale))) {
      stop(
        "`parameters$gamma` and `parameters$scale` must agree when both ",
        "are supplied.",
        call. = FALSE
      )
    }
    if (is.null(gamma)) {
      gamma <- scale
    }
    if (is.null(gamma)) {
      stop(
        "`parameters` must contain `gamma` or `scale` for model ",
        "\"vargamma\".",
        call. = FALSE
      )
    }

    gamma <- .tensor_rosenblatt_positive_parameter(
      gamma,
      if (is.null(parameters$gamma)) "scale" else "gamma",
      model
    )
    return(list(lambda = gamma, chi = 0, psi = 2 * gamma))
  }

  if (model == "invgauss") {
    kappa <- .tensor_rosenblatt_positive_parameter(
      .tensor_rosenblatt_parameter(parameters, "kappa", model),
      "kappa",
      model
    )
    return(list(lambda = -0.5, chi = 1, psi = kappa^2))
  }

  lambda <- .tensor_rosenblatt_parameter(parameters, "lambda", model)
  omega <- .tensor_rosenblatt_positive_parameter(
    .tensor_rosenblatt_parameter(parameters, "omega", model),
    "omega",
    model
  )
  if (length(lambda) != 1L || !is.numeric(lambda) || is.na(lambda) ||
      !is.finite(lambda)) {
    stop(
      "`parameters$lambda` must be one finite number for model ",
      "\"genhyper\".",
      call. = FALSE
    )
  }

  list(lambda = as.numeric(lambda), chi = omega, psi = omega)
}

.tensor_rosenblatt_validate_controls <- function(rel.tol,
                                                  abs.tol,
                                                  subdivisions,
                                                  show_progress) {
  if (length(rel.tol) != 1L || !is.numeric(rel.tol) || is.na(rel.tol) ||
      !is.finite(rel.tol) || rel.tol <= 0) {
    stop("`rel.tol` must be one positive finite number.", call. = FALSE)
  }
  if (length(abs.tol) != 1L || !is.numeric(abs.tol) || is.na(abs.tol) ||
      !is.finite(abs.tol) || abs.tol <= 0) {
    stop("`abs.tol` must be one positive finite number.", call. = FALSE)
  }
  if (length(subdivisions) != 1L || !is.numeric(subdivisions) ||
      is.na(subdivisions) || !is.finite(subdivisions) ||
      subdivisions < 1L || subdivisions != floor(subdivisions)) {
    stop("`subdivisions` must be one positive integer.", call. = FALSE)
  }
  if (length(show_progress) != 1L || is.na(show_progress) ||
      !is.logical(show_progress)) {
    stop("`show_progress` must be TRUE or FALSE.", call. = FALSE)
  }
}

.tensor_rosenblatt_whitening <- function(sigmas) {
  lapply(seq_along(sigmas), function(mode) {
    sigma <- (sigmas[[mode]] + t(sigmas[[mode]])) / 2
    lower_chol <- tryCatch(
      t(chol(sigma)),
      error = function(error) {
        stop(
          "`sigmas[[", mode, "]]` is not positive definite: ",
          conditionMessage(error),
          call. = FALSE
        )
      }
    )
    forwardsolve(lower_chol, diag(nrow(lower_chol)))
  })
}

.tensor_rosenblatt_whiten <- function(value, whitening) {
  for (mode in seq_along(whitening)) {
    value <- n_prod(value, mode, whitening[[mode]])
  }
  value
}

.tensor_rosenblatt_output <- function(values, draws) {
  output <- array(
    values,
    dim = c(n_draws(draws), draw_shape(draws)),
    dimnames = dimnames(unclass(draws))
  )
  .new_tensor_array(output)
}

.tensor_rosenblatt_normal <- function(draws,
                                      mu,
                                      whitening,
                                      show_progress) {
  n <- n_draws(draws)
  p <- prod(draw_shape(draws))
  uniforms <- matrix(NA_real_, nrow = n, ncol = p)

  for (draw_index in seq_len(n)) {
    current <- .tensor_single_draw_array(
      pull_draw(draws, draw_index)
    ) - mu
    current <- .tensor_rosenblatt_whiten(current, whitening)
    uniforms[draw_index, ] <- stats::pnorm(as.numeric(current))

    if (show_progress) {
      message("Transformed draw ", draw_index, " of ", n)
    }
  }

  .tensor_rosenblatt_output(uniforms, draws)
}

.tensor_rosenblatt_mixture <- function(draws,
                                       mu,
                                       skew,
                                       whitening,
                                       lambda,
                                       chi,
                                       psi,
                                       rel.tol,
                                       abs.tol,
                                       subdivisions,
                                       show_progress) {
  whitened_skew <- as.numeric(
    .tensor_rosenblatt_whiten(skew, whitening)
  )
  n <- n_draws(draws)
  p <- prod(draw_shape(draws))
  uniforms <- matrix(NA_real_, nrow = n, ncol = p)

  for (draw_index in seq_len(n)) {
    current <- .tensor_single_draw_array(
      pull_draw(draws, draw_index)
    ) - mu
    current <- as.numeric(
      .tensor_rosenblatt_whiten(current, whitening)
    )

    cumulative_squared_draw <- 0
    cumulative_squared_skew <- 0

    for (coordinate in seq_len(p)) {
      previous_coordinates <- coordinate - 1L

      uniforms[draw_index, coordinate] <- tryCatch(
        .tensor_rosenblatt_gig_cdf(
          value = current[coordinate],
          skew = whitened_skew[coordinate],
          lambda = lambda - previous_coordinates / 2,
          chi = chi + cumulative_squared_draw,
          psi = psi + cumulative_squared_skew,
          rel.tol = rel.tol,
          abs.tol = abs.tol,
          subdivisions = subdivisions
        ),
        error = function(error) {
          stop(
            "Quadrature failed for draw ", draw_index,
            ", coordinate ", coordinate, ": ",
            conditionMessage(error),
            call. = FALSE
          )
        }
      )

      cumulative_squared_draw <-
        cumulative_squared_draw + current[coordinate]^2
      cumulative_squared_skew <-
        cumulative_squared_skew + whitened_skew[coordinate]^2
    }

    if (show_progress) {
      message("Transformed draw ", draw_index, " of ", n)
    }
  }

  .tensor_rosenblatt_output(uniforms, draws)
}

.tensor_rosenblatt_gig_log_center <- function(lambda, chi, psi) {
  if (chi == 0) {
    if (lambda <= 0 || psi <= 0) {
      stop("Invalid gamma-boundary GIG parameters.", call. = FALSE)
    }
    return(log(2 * lambda / psi))
  }

  if (psi == 0) {
    if (lambda >= 0 || chi <= 0) {
      stop("Invalid inverse-gamma-boundary GIG parameters.", call. = FALSE)
    }
    return(log(-chi / (2 * lambda)))
  }

  if (chi < 0 || psi < 0) {
    stop("GIG parameters `chi` and `psi` cannot be negative.",
         call. = FALSE)
  }

  root <- sqrt(lambda^2 + chi * psi)

  if (lambda >= 0) {
    log(root + lambda) - log(psi)
  } else {
    log(chi) - log(root - lambda)
  }
}

.tensor_rosenblatt_gig_log_bounds <- function(lambda,
                                              chi,
                                              psi,
                                              log_center,
                                              rel.tol,
                                              abs.tol) {
  center <- exp(log_center)
  scaled_chi <- if (chi == 0) 0 else chi / center
  scaled_psi <- if (psi == 0) 0 else psi * center
  log_drop <- max(
    40,
    -log(min(rel.tol, abs.tol)) + 10
  )

  relative_log_kernel <- function(log_ratio) {
    lower_tail <- if (scaled_chi == 0) {
      0
    } else {
      scaled_chi * expm1(-log_ratio)
    }
    upper_tail <- if (scaled_psi == 0) {
      0
    } else {
      scaled_psi * expm1(log_ratio)
    }
    value <- lambda * log_ratio - (lower_tail + upper_tail) / 2

    if (is.nan(value)) -Inf else value
  }

  expand_bound <- function(direction) {
    width <- 1

    while (width < 1024 &&
           relative_log_kernel(direction * width) > -log_drop) {
      width <- width * 2
    }

    width
  }

  c(
    lower = -expand_bound(-1),
    upper = expand_bound(1)
  )
}

.tensor_rosenblatt_gig_cdf <- function(value,
                                       skew,
                                       lambda,
                                       chi,
                                       psi,
                                       rel.tol,
                                       abs.tol,
                                       subdivisions) {
  log_center <- .tensor_rosenblatt_gig_log_center(
    lambda,
    chi,
    psi
  )
  log_min <- log(.Machine$double.xmin)
  log_max <- log(.Machine$double.xmax)

  integrand <- function(log_ratio) {
    answer <- numeric(length(log_ratio))
    log_w <- log_center + log_ratio
    valid <- is.finite(log_w) & log_w > log_min & log_w < log_max

    if (any(valid)) {
      current_log_w <- log_w[valid]
      current_w <- exp(current_log_w)
      normal_argument <- (
        value - current_w * skew
      ) / sqrt(current_w)

      log_integrand <-
        stats::pnorm(normal_argument, log.p = TRUE) +
        GIGrvg::dgig(
          current_w,
          lambda = lambda,
          chi = chi,
          psi = psi,
          log = TRUE
        ) +
        current_log_w

      answer[valid] <- exp(log_integrand)
      answer[!is.finite(answer)] <- 0
    }

    answer
  }

  result <- tryCatch(
    stats::integrate(
      integrand,
      lower = -Inf,
      upper = Inf,
      rel.tol = rel.tol,
      abs.tol = abs.tol,
      subdivisions = subdivisions,
      stop.on.error = TRUE
    ),
    error = function(infinite_range_error) {
      bounds <- .tensor_rosenblatt_gig_log_bounds(
        lambda = lambda,
        chi = chi,
        psi = psi,
        log_center = log_center,
        rel.tol = rel.tol,
        abs.tol = abs.tol
      )
      fallback_subdivisions <- max(200L, 2L * subdivisions)

      finite_piece <- function(lower, upper) {
        stats::integrate(
          integrand,
          lower = lower,
          upper = upper,
          rel.tol = rel.tol,
          abs.tol = abs.tol / 2,
          subdivisions = fallback_subdivisions,
          stop.on.error = TRUE
        )
      }

      tryCatch(
        {
          left <- finite_piece(bounds[["lower"]], 0)
          right <- finite_piece(0, bounds[["upper"]])
          left$value <- left$value + right$value
          left$abs.error <- left$abs.error + right$abs.error
          left
        },
        error = function(finite_range_error) {
          stop(
            conditionMessage(infinite_range_error),
            "; finite-range fallback also failed: ",
            conditionMessage(finite_range_error),
            call. = FALSE
          )
        }
      )
    }
  )

  probability <- result$value
  if (!is.finite(probability)) {
    stop("The conditional CDF was not finite.", call. = FALSE)
  }
  if (probability < -abs.tol || probability > 1 + abs.tol) {
    stop(
      "The conditional CDF fell outside [0, 1].",
      call. = FALSE
    )
  }

  min(max(probability, 0), 1)
}
