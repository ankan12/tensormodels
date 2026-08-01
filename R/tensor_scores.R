#' Per-observation scores for fitted tensor distributions
#'
#' Computes one likelihood-score row for each tensor observation. For the
#' normal model the scores are analytic. For the normal variance-mean mixture
#' models, the function uses Fisher's identity and one evaluation of the
#' conditional mixing-variable moments at the supplied fit. It does not
#' perturb and re-evaluate every model parameter.
#'
#' @param draws A `tensor` object containing IID observations.
#' @param model The fitted distribution: `"normal"`, `"skewt"`,
#'   `"vargamma"`, `"invgauss"`, or `"genhyper"`.
#' @param fit An optional parameter list returned by [tensor_mle()]. When
#'   `NULL`, the model is fitted once.
#' @param max_iter,tol,quiet,restrict Arguments passed to [tensor_mle()] when
#'   `fit` is `NULL`. Restricted covariance fits are not currently supported.
#'
#' @return A numeric matrix with one row per observation and one column per
#'   identifiable fitted parameter. The fitted parameters and model are
#'   available in the `"fit"` and `"model"` attributes.
#'
#' @details
#' The covariance factors use the same identifiable working parameterization
#' as [tensor_score_gof()]: Cholesky coordinates, with unit mean diagonal for
#' all covariance factors except the last. Shape parameters constrained to be
#' positive are represented on the log scale.
#'
#' The returned scores are evaluated at `fit`. Consequently, column means
#' should be close to zero when `fit` is a sufficiently converged maximum
#' likelihood estimate.
#'
#' @examples
#' draws <- rtnorm(
#'   50,
#'   mu = array(0, dim = c(2, 2)),
#'   sigmas = list(diag(2), diag(2))
#' )
#' fit <- tensor_mle(draws, model = "normal")
#' scores <- tensor_scores(draws, model = "normal", fit = fit)
#' dim(scores)
#'
#' @export
tensor_scores <- function(
    draws,
    model = c("normal", "skewt", "vargamma", "invgauss", "genhyper"),
    fit = NULL,
    max_iter = 1000L,
    tol = 1e-6,
    quiet = TRUE,
    restrict = NULL) {
  model <- match.arg(model)

  if (!inherits(draws, "tensor")) {
    stop("`draws` must be a `tensor` object.", call. = FALSE)
  }
  if (n_draws(draws) < 1L) {
    stop("`draws` must contain at least one observation.", call. = FALSE)
  }
  if (length(restrict) > 0L) {
    stop(
      "Restricted covariance fits are not currently supported by ",
      "`tensor_scores()`.",
      call. = FALSE
    )
  }

  if (is.null(fit)) {
    fit <- tensor_mle(
      draws,
      model = model,
      max_iter = max_iter,
      tol = tol,
      quiet = quiet,
      restrict = restrict
    )
  }
  if (!is.list(fit)) {
    stop("`fit` must be `NULL` or a fitted parameter list.", call. = FALSE)
  }

  dims <- draw_shape(draws)
  parameterization <- .tensor_score_gof_parameterization(
    fit = fit,
    model = model,
    dims = dims
  )
  parameters <- parameterization$decode(parameterization$eta)
  scores <- .tensor_scores_analytic(
    draws = draws,
    model = model,
    parameters = parameters,
    parameter_names = parameterization$names
  )

  attr(scores, "fit") <- parameters
  attr(scores, "model") <- model
  scores
}

.tensor_scores_posterior_moments <- function(delta, rho, model, parameters,
                                             tensor_dimension) {
  n <- length(delta)

  if (model == "normal") {
    return(list(
      Ew = rep(1, n),
      Einvw = rep(1, n),
      Elogw = rep(0, n)
    ))
  }

  specification <- switch(
    model,
    skewt = list(
      lambda = -(parameters$nu + tensor_dimension) / 2,
      chi = delta + parameters$nu,
      psi = rep(rho, n)
    ),
    vargamma = list(
      lambda = parameters$gamma - tensor_dimension / 2,
      chi = delta,
      psi = rep(rho + 2 * parameters$gamma, n)
    ),
    invgauss = list(
      lambda = -(1 + tensor_dimension) / 2,
      chi = delta + 1,
      psi = rep(rho + parameters$kappa^2, n)
    ),
    genhyper = list(
      lambda = parameters$lambda - tensor_dimension / 2,
      chi = delta + parameters$omega,
      psi = rep(rho + parameters$omega, n)
    )
  )

  lambda <- specification$lambda
  chi <- specification$chi
  psi <- specification$psi
  Ew <- Einvw <- Elogw <- rep(NA_real_, n)

  ordinary <- chi > 0 & psi > 0
  if (any(ordinary)) {
    argument <- sqrt(chi[ordinary] * psi[ordinary])
    ratio <- .besselK_asym_ratio(
      argument,
      numerator_nu = lambda + 1,
      denominator_nu = lambda
    )
    Ew[ordinary] <- sqrt(chi[ordinary] / psi[ordinary]) * ratio
    Einvw[ordinary] <-
      sqrt(psi[ordinary] / chi[ordinary]) * ratio -
      2 * lambda / chi[ordinary]
    Elogw[ordinary] <-
      0.5 * log(chi[ordinary] / psi[ordinary]) +
      .dlog_besselK_asym_dnu(argument, lambda, eps = 1e-5)
  }

  inverse_gamma <- chi > 0 & psi == 0 & lambda < 0
  if (any(inverse_gamma)) {
    shape <- -lambda
    rate <- chi[inverse_gamma] / 2
    if (shape <= 1) {
      stop(
        "The fitted conditional mixing distribution has no finite mean.",
        call. = FALSE
      )
    }
    Ew[inverse_gamma] <- rate / (shape - 1)
    Einvw[inverse_gamma] <- shape / rate
    Elogw[inverse_gamma] <- log(rate) - digamma(shape)
  }

  gamma_case <- chi == 0 & psi > 0 & lambda > 0
  if (any(gamma_case)) {
    shape <- lambda
    rate <- psi[gamma_case] / 2
    Ew[gamma_case] <- shape / rate
    Einvw[gamma_case] <- if (shape > 1) {
      rate / (shape - 1)
    } else {
      Inf
    }
    Elogw[gamma_case] <- digamma(shape) - log(rate)
  }

  moments <- cbind(Ew = Ew, Einvw = Einvw, Elogw = Elogw)
  if (any(!is.finite(moments))) {
    stop(
      "Could not calculate finite conditional mixing-variable moments at ",
      "the supplied fit.",
      call. = FALSE
    )
  }

  list(Ew = Ew, Einvw = Einvw, Elogw = Elogw)
}

.tensor_scores_covariance_derivatives <- function(sigma, constrained) {
  dimension <- nrow(sigma)
  lower <- t(chol((sigma + t(sigma)) / 2))
  indices <- which(
    lower.tri(matrix(0, dimension, dimension), diag = TRUE),
    arr.ind = TRUE
  )
  if (constrained) {
    indices <- indices[
      !(indices[, 1L] == 1L & indices[, 2L] == 1L),
      ,
      drop = FALSE
    ]
  }

  derivatives <- vector("list", nrow(indices))
  raw_sigma <- tcrossprod(lower)
  raw_scale <- mean(diag(raw_sigma))

  for (index in seq_len(nrow(indices))) {
    row <- indices[index, 1L]
    column <- indices[index, 2L]
    lower_derivative <- matrix(0, dimension, dimension)
    lower_derivative[row, column] <- if (row == column) {
      lower[row, column]
    } else {
      1
    }
    raw_derivative <-
      lower_derivative %*% t(lower) +
      lower %*% t(lower_derivative)

    derivatives[[index]] <- if (constrained) {
      raw_derivative / raw_scale -
        raw_sigma * mean(diag(raw_derivative)) / raw_scale^2
    } else {
      raw_derivative
    }
  }

  derivatives
}

.tensor_scores_analytic <- function(draws, model, parameters,
                                    parameter_names) {
  n <- n_draws(draws)
  dims <- draw_shape(draws)
  modes <- length(dims)
  tensor_dimension <- prod(dims)
  inverse_sigmas <- lapply(parameters$sigmas, invert_safe)
  residuals <- vector("list", n)
  precision_residuals <- vector("list", n)
  delta <- numeric(n)

  for (observation in seq_len(n)) {
    residual <- .tensor_single_draw_array(
      pull_draw(draws, observation)
    ) - parameters$mu
    precision_residual <- residual
    for (mode in seq_len(modes)) {
      precision_residual <- n_prod(
        precision_residual,
        mode,
        inverse_sigmas[[mode]]
      )
    }
    residuals[[observation]] <- residual
    precision_residuals[[observation]] <- precision_residual
    delta[observation] <- max(0, sum(residual * precision_residual))
  }

  if (model == "normal") {
    skew <- array(0, dim = dims)
    precision_skew <- skew
    rho <- 0
  } else {
    skew <- parameters$skew
    precision_skew <- skew
    for (mode in seq_len(modes)) {
      precision_skew <- n_prod(
        precision_skew,
        mode,
        inverse_sigmas[[mode]]
      )
    }
    rho <- max(0, sum(skew * precision_skew))
  }

  moments <- .tensor_scores_posterior_moments(
    delta = delta,
    rho = rho,
    model = model,
    parameters = parameters,
    tensor_dimension = tensor_dimension
  )

  scores <- matrix(
    0,
    nrow = n,
    ncol = length(parameter_names),
    dimnames = list(NULL, parameter_names)
  )

  mu_columns <- seq_len(tensor_dimension)
  for (observation in seq_len(n)) {
    scores[observation, mu_columns] <-
      moments$Einvw[observation] *
      as.numeric(precision_residuals[[observation]]) -
      as.numeric(precision_skew)
  }

  cursor <- tensor_dimension
  if (model != "normal") {
    skew_columns <- cursor + seq_len(tensor_dimension)
    for (observation in seq_len(n)) {
      scores[observation, skew_columns] <-
        as.numeric(precision_residuals[[observation]]) -
        moments$Ew[observation] * as.numeric(precision_skew)
    }
    cursor <- cursor + tensor_dimension
  }

  for (mode in seq_len(modes)) {
    dimension <- dims[mode]
    other_modes <- setdiff(seq_len(modes), mode)
    skew_other_precision <- skew
    for (other_mode in other_modes) {
      skew_other_precision <- n_prod(
        skew_other_precision,
        other_mode,
        inverse_sigmas[[other_mode]]
      )
    }
    skew_flat <- matricization(skew, mode)
    skew_precision_flat <- matricization(skew_other_precision, mode)
    derivatives <- .tensor_scores_covariance_derivatives(
      parameters$sigmas[[mode]],
      constrained = mode < modes
    )
    covariance_columns <- cursor + seq_along(derivatives)

    for (observation in seq_len(n)) {
      residual <- residuals[[observation]]
      residual_other_precision <- residual
      for (other_mode in other_modes) {
        residual_other_precision <- n_prod(
          residual_other_precision,
          other_mode,
          inverse_sigmas[[other_mode]]
        )
      }
      residual_flat <- matricization(residual, mode)
      residual_precision_flat <-
        matricization(residual_other_precision, mode)

      expected_crossproduct <-
        moments$Einvw[observation] *
        tcrossprod(residual_precision_flat, residual_flat) -
        tcrossprod(residual_precision_flat, skew_flat) -
        tcrossprod(skew_precision_flat, residual_flat) +
        moments$Ew[observation] *
        tcrossprod(skew_precision_flat, skew_flat)

      inverse_sigma <- inverse_sigmas[[mode]]
      covariance_score <- 0.5 * inverse_sigma %*%
        (
          expected_crossproduct -
          (tensor_dimension / dimension) * parameters$sigmas[[mode]]
        ) %*%
        inverse_sigma
      covariance_score <-
        (covariance_score + t(covariance_score)) / 2

      if (length(derivatives) > 0L) {
        scores[observation, covariance_columns] <- vapply(
          derivatives,
          function(derivative) sum(covariance_score * derivative),
          numeric(1)
        )
      }
    }
    cursor <- cursor + length(derivatives)
  }

  if (model == "skewt") {
    nu <- parameters$nu
    scores[, cursor + 1L] <- nu / 2 * (
      log(nu / 2) + 1 - digamma(nu / 2) -
        moments$Elogw - moments$Einvw
    )
  } else if (model == "vargamma") {
    gamma <- parameters$gamma
    scores[, cursor + 1L] <- gamma * (
      log(gamma) + 1 - digamma(gamma) +
        moments$Elogw - moments$Ew
    )
  } else if (model == "invgauss") {
    kappa <- parameters$kappa
    scores[, cursor + 1L] <-
      kappa * (1 - kappa * moments$Ew)
  } else if (model == "genhyper") {
    lambda <- parameters$lambda
    omega <- parameters$omega
    scores[, cursor + 1L] <-
      moments$Elogw -
      .dlog_besselK_asym_dnu(omega, lambda, eps = 1e-5)
    ratio_plus <- .besselK_asym_ratio(
      omega,
      numerator_nu = lambda + 1,
      denominator_nu = lambda
    )
    ratio_minus <- .besselK_asym_ratio(
      omega,
      numerator_nu = lambda - 1,
      denominator_nu = lambda
    )
    scores[, cursor + 2L] <- omega / 2 * (
      ratio_plus + ratio_minus -
        moments$Ew - moments$Einvw
    )
  }

  if (any(!is.finite(scores))) {
    stop("The analytic score matrix contains non-finite values.",
         call. = FALSE)
  }
  scores
}
