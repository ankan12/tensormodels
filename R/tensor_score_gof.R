#' Score-adjusted Rosenblatt smooth goodness-of-fit test
#'
#' Performs an experimental finite-basis goodness-of-fit test for any of the
#' tensor-variate distributions supported by [tensor_mle()]. The test applies
#' [tensor_rosenblatt()], evaluates orthonormal smooth functions of the joint
#' uniform vectors, and projects those functions away from the fitted model's
#' nuisance-score space.
#'
#' @param draws A `tensor` object containing independent tensor observations.
#' @param model One of `"normal"`, `"skewt"`, `"vargamma"`, `"invgauss"`,
#'   or `"genhyper"`.
#' @param fit An optional fitted parameter list returned by [tensor_mle()]. If
#'   `NULL`, the model is fitted internally.
#' @param marginal_degree Positive integer giving the highest shifted-Legendre
#'   degree included separately for each Rosenblatt coordinate.
#' @param include_interactions Whether to include products of first-degree
#'   basis functions for pairs of Rosenblatt coordinates. These terms target
#'   residual dependence.
#' @param max_basis Maximum number of smooth basis functions. When the complete
#'   requested basis is larger, marginal terms are retained first and a
#'   deterministic spread of interaction terms is used.
#' @param fd_step Relative central finite-difference step used to calculate
#'   per-observation likelihood scores.
#' @param eigen_tol Relative eigenvalue tolerance used for generalized
#'   inverses and numerical ranks.
#' @param max_iter,tol,quiet,restrict Arguments passed to [tensor_mle()] when
#'   `fit = NULL`. Restricted covariance fits are not currently supported.
#' @param rel.tol,abs.tol,subdivisions Numerical integration controls passed to
#'   [tensor_rosenblatt()] for skewed models.
#'
#' @return An object inheriting from `"htest"`. In addition to the usual test
#'   fields, the object contains the fitted model, Rosenblatt values, basis
#'   names and means, score diagnostics, and the estimated adjusted covariance.
#'
#' @details
#' Let \eqn{U_i} be the joint Rosenblatt vector for tensor observation \eqn{i}
#' and let \eqn{m(U_i)} contain zero-mean orthonormal shifted-Legendre basis
#' functions. Let \eqn{s_i} be the fitted model's likelihood score,
#' \eqn{I=E(s_i s_i^\mathsf{T})}, and
#' \eqn{C=E\{m(U_i)s_i^\mathsf{T}\}}. Estimating the model parameters changes
#' the limiting covariance of the smooth empirical moments. The function uses
#' score-residualized moments to estimate the adjusted covariance
#' \eqn{\Omega}, and calculates
#' \deqn{T=n\bar m^\mathsf{T}\widehat\Omega^+\bar m,}
#' where \eqn{\bar m} is evaluated at the fitted parameter estimates.
#' Its reference distribution is asymptotically chi-squared with degrees of
#' freedom equal to the numerical rank of \eqn{\widehat\Omega}.
#'
#' The likelihood scores are calculated in an identifiable working
#' parameterization: tensor locations and skews are unconstrained, positive
#' shape parameters and Cholesky diagonals use logarithms, and the first
#' covariance factors have fixed trace with overall scale retained in the last
#' factor.
#'
#' This test is experimental. Its chi-squared calibration requires a regular,
#' root-\eqn{n} consistent fit, a full-rank nuisance information matrix, and a
#' fixed basis dimension that is small relative to the sample size. It should
#' not be treated as formally calibrated when the MLE has not converged, the
#' score information is rank deficient, or a shape parameter is effectively
#' on a boundary.
#'
#' @examples
#' set.seed(1)
#' mu <- array(0, dim = 2)
#' sigmas <- list(diag(2))
#' draws <- rtnorm(200, mu = mu, sigmas = sigmas)
#' tensor_score_gof(draws, model = "normal", marginal_degree = 2)
#'
#' @export
tensor_score_gof <- function(
    draws,
    model = c("normal", "skewt", "vargamma", "invgauss", "genhyper"),
    fit = NULL,
    marginal_degree = 2L,
    include_interactions = TRUE,
    max_basis = 100L,
    fd_step = 1e-5,
    eigen_tol = 1e-8,
    max_iter = 1000L,
    tol = 1e-6,
    quiet = TRUE,
    restrict = NULL,
    rel.tol = 1e-7,
    abs.tol = 1e-9,
    subdivisions = 100L) {
  model <- match.arg(model)
  data_name <- deparse(substitute(draws))

  if (!inherits(draws, "tensor")) {
    stop("`draws` must be a `tensor` object.", call. = FALSE)
  }
  n <- n_draws(draws)
  if (n < 3L) {
    stop("At least three tensor observations are required.", call. = FALSE)
  }
  if (length(restrict) > 0L) {
    stop(
      "Restricted covariance fits are not currently supported by ",
      "`tensor_score_gof()`.",
      call. = FALSE
    )
  }
  .tensor_score_gof_validate_controls(
    marginal_degree = marginal_degree,
    include_interactions = include_interactions,
    max_basis = max_basis,
    fd_step = fd_step,
    eigen_tol = eigen_tol
  )

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

  parameterization <- .tensor_score_gof_parameterization(
    fit = fit,
    model = model,
    dims = draw_shape(draws)
  )
  fitted_parameters <- parameterization$decode(parameterization$eta)

  uniforms <- tensor_rosenblatt(
    draws = draws,
    model = model,
    parameters = fitted_parameters,
    rel.tol = rel.tol,
    abs.tol = abs.tol,
    subdivisions = subdivisions
  )
  uniform_matrix <- matrix(
    as.numeric(unclass(uniforms)),
    nrow = n,
    ncol = prod(draw_shape(draws))
  )

  effective_max_basis <- min(as.integer(max_basis), n - 1L)
  if (effective_max_basis < max_basis) {
    warning(
      "`max_basis` was reduced to ", effective_max_basis,
      " so that the adjusted covariance can be estimated.",
      call. = FALSE
    )
  }
  basis <- .tensor_score_gof_basis(
    uniforms = uniform_matrix,
    marginal_degree = as.integer(marginal_degree),
    include_interactions = include_interactions,
    max_basis = effective_max_basis
  )

  scores <- .tensor_score_gof_scores(
    draws = draws,
    model = model,
    parameterization = parameterization,
    fd_step = fd_step
  )
  q <- ncol(scores)
  if (n <= q) {
    warning(
      "The sample size (", n, ") is not larger than the score dimension (",
      q, "); the chi-squared calibration is not reliable.",
      call. = FALSE
    )
  }

  centered_scores <- sweep(scores, 2L, colMeans(scores), "-")
  centered_basis <- sweep(basis$values, 2L, colMeans(basis$values), "-")

  information <- crossprod(centered_scores) / n
  information_inverse <- .tensor_score_gof_eigen_inverse(
    information,
    eigen_tol
  )
  if (information_inverse$rank < q) {
    warning(
      "The nuisance-score information matrix is rank deficient (rank ",
      information_inverse$rank, " of ", q,
      "); interpret the chi-squared p-value cautiously.",
      call. = FALSE
    )
  }

  cross_covariance <- crossprod(centered_basis, centered_scores) / n
  projection <- cross_covariance %*% information_inverse$inverse
  adjusted_basis <- basis$values - scores %*% t(projection)
  basis_mean <- colMeans(basis$values)
  centered_adjusted <- sweep(
    adjusted_basis,
    2L,
    colMeans(adjusted_basis),
    "-"
  )
  adjusted_covariance <- crossprod(centered_adjusted) / n

  covariance_inverse <- .tensor_score_gof_eigen_inverse(
    adjusted_covariance,
    eigen_tol
  )
  if (covariance_inverse$rank < 1L) {
    stop(
      "The score-adjusted basis covariance has numerical rank zero.",
      call. = FALSE
    )
  }

  statistic <- as.numeric(
    n * crossprod(
      basis_mean,
      covariance_inverse$inverse %*% basis_mean
    )
  )
  degrees_freedom <- covariance_inverse$rank
  p_value <- stats::pchisq(
    statistic,
    df = degrees_freedom,
    lower.tail = FALSE
  )

  score_scales <- sqrt(pmax(diag(information), 0) / n)
  standardized_score_mean <- rep(0, q)
  usable_scales <- score_scales > 0 & is.finite(score_scales)
  standardized_score_mean[usable_scales] <-
    colMeans(scores)[usable_scales] / score_scales[usable_scales]
  stationarity <- max(abs(standardized_score_mean))

  if (is.finite(stationarity) && stationarity > 1) {
    warning(
      "The fitted score is not close to zero (maximum standardized score ",
      sprintf("%.2f", stationarity),
      "); the optimizer may not have converged sufficiently.",
      call. = FALSE
    )
  }

  model_label <- c(
    normal = "tensor normal",
    skewt = "tensor skew-t",
    vargamma = "tensor variance-gamma",
    invgauss = "tensor inverse-Gaussian",
    genhyper = "tensor generalized hyperbolic"
  )[[model]]

  structure(
    list(
      statistic = c(T = statistic),
      parameter = c(df = degrees_freedom),
      p.value = p_value,
      method = paste0(
        "Experimental score-adjusted Rosenblatt smooth GOF test for the ",
        model_label
      ),
      data.name = data_name,
      model = model,
      fit = fitted_parameters,
      uniforms = uniforms,
      basis_names = basis$names,
      basis_means = basis_mean,
      adjusted_basis_means = colMeans(adjusted_basis),
      adjusted_covariance = adjusted_covariance,
      score_dimension = q,
      information_rank = information_inverse$rank,
      adjusted_covariance_rank = covariance_inverse$rank,
      score_stationarity = stationarity,
      scores = scores,
      parameter_names = parameterization$names
    ),
    class = c("tensor_score_gof_test", "htest")
  )
}

.tensor_score_gof_validate_controls <- function(marginal_degree,
                                                include_interactions,
                                                max_basis,
                                                fd_step,
                                                eigen_tol) {
  if (length(marginal_degree) != 1L ||
      !is.numeric(marginal_degree) ||
      !is.finite(marginal_degree) ||
      marginal_degree < 1L ||
      marginal_degree != floor(marginal_degree)) {
    stop("`marginal_degree` must be one positive integer.", call. = FALSE)
  }
  if (length(include_interactions) != 1L ||
      is.na(include_interactions) ||
      !is.logical(include_interactions)) {
    stop("`include_interactions` must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(max_basis) != 1L ||
      !is.numeric(max_basis) ||
      !is.finite(max_basis) ||
      max_basis < 1L ||
      max_basis != floor(max_basis)) {
    stop("`max_basis` must be one positive integer.", call. = FALSE)
  }
  if (length(fd_step) != 1L || !is.numeric(fd_step) ||
      !is.finite(fd_step) || fd_step <= 0) {
    stop("`fd_step` must be one positive finite number.", call. = FALSE)
  }
  if (length(eigen_tol) != 1L || !is.numeric(eigen_tol) ||
      !is.finite(eigen_tol) || eigen_tol <= 0 || eigen_tol >= 1) {
    stop("`eigen_tol` must be strictly between zero and one.",
         call. = FALSE)
  }
}

.tensor_score_gof_canonicalize <- function(parameters, dims, model) {
  required <- c("mu", "sigmas")
  if (model != "normal") required <- c(required, "skew")
  missing <- required[vapply(
    required,
    function(name) is.null(parameters[[name]]),
    logical(1)
  )]
  if (length(missing) > 0L) {
    stop(
      "`fit` is missing required parameter",
      if (length(missing) > 1L) "s: " else ": ",
      paste(missing, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  output <- parameters
  .validate_same_dims(output$mu, dims, "fit$mu", "draws")
  output$mu <- array(as.numeric(output$mu), dim = dims)
  if (model != "normal") {
    .validate_same_dims(output$skew, dims, "fit$skew", "draws")
    output$skew <- array(as.numeric(output$skew), dim = dims)
  }
  output$sigmas <- .prepare_sigmas(output$sigmas, dims)

  modes <- length(dims)
  scale_product <- 1
  if (modes > 1L) {
    for (mode in seq_len(modes - 1L)) {
      scale <- mean(diag(output$sigmas[[mode]]))
      output$sigmas[[mode]] <- output$sigmas[[mode]] / scale
      scale_product <- scale_product * scale
    }
    output$sigmas[[modes]] <-
      output$sigmas[[modes]] * scale_product
  }

  if (model == "skewt") {
    .tensor_rosenblatt_positive_parameter(output$nu, "nu", model)
  } else if (model == "vargamma") {
    if (is.null(output$gamma)) output$gamma <- output$scale
    output$gamma <- .tensor_rosenblatt_positive_parameter(
      output$gamma,
      "gamma",
      model
    )
  } else if (model == "invgauss") {
    output$kappa <- .tensor_rosenblatt_positive_parameter(
      output$kappa,
      "kappa",
      model
    )
  } else if (model == "genhyper") {
    if (length(output$lambda) != 1L ||
        !is.numeric(output$lambda) ||
        !is.finite(output$lambda)) {
      stop("`fit$lambda` must be one finite number.", call. = FALSE)
    }
    output$omega <- .tensor_rosenblatt_positive_parameter(
      output$omega,
      "omega",
      model
    )
  }

  output
}

.tensor_score_gof_parameterization <- function(fit, model, dims) {
  base <- .tensor_score_gof_canonicalize(fit, dims, model)
  eta <- as.numeric(base$mu)
  parameter_names <- paste0("mu[", seq_along(eta), "]")
  p <- length(eta)

  if (model != "normal") {
    eta <- c(eta, as.numeric(base$skew))
    parameter_names <- c(
      parameter_names,
      paste0("skew[", seq_len(p), "]")
    )
  }

  lower_chols <- lapply(base$sigmas, function(sigma) {
    t(chol((sigma + t(sigma)) / 2))
  })
  covariance_specs <- vector("list", length(dims))

  for (mode in seq_along(dims)) {
    dimension <- dims[mode]
    indices <- which(lower.tri(
      matrix(0, dimension, dimension),
      diag = TRUE
    ), arr.ind = TRUE)
    if (mode < length(dims)) {
      indices <- indices[
        !(indices[, 1L] == 1L & indices[, 2L] == 1L),
        ,
        drop = FALSE
      ]
    }

    values <- numeric(nrow(indices))
    names_current <- character(nrow(indices))
    for (index in seq_len(nrow(indices))) {
      row <- indices[index, 1L]
      column <- indices[index, 2L]
      diagonal <- row == column
      values[index] <- if (diagonal) {
        log(lower_chols[[mode]][row, column])
      } else {
        lower_chols[[mode]][row, column]
      }
      names_current[index] <- paste0(
        "sigma", mode, "_chol[", row, ",", column, "]",
        if (diagonal) "_log" else ""
      )
    }

    covariance_specs[[mode]] <- list(
      indices = indices,
      start = length(eta) + 1L,
      end = length(eta) + length(values)
    )
    eta <- c(eta, values)
    parameter_names <- c(parameter_names, names_current)
  }

  shape_start <- length(eta) + 1L
  if (model == "skewt") {
    eta <- c(eta, log(base$nu))
    parameter_names <- c(parameter_names, "log_nu")
  } else if (model == "vargamma") {
    eta <- c(eta, log(base$gamma))
    parameter_names <- c(parameter_names, "log_gamma")
  } else if (model == "invgauss") {
    eta <- c(eta, log(base$kappa))
    parameter_names <- c(parameter_names, "log_kappa")
  } else if (model == "genhyper") {
    eta <- c(eta, base$lambda, log(base$omega))
    parameter_names <- c(parameter_names, "lambda", "log_omega")
  }

  decode <- function(current_eta) {
    output <- base
    cursor <- 0L

    output$mu <- array(current_eta[seq_len(p)], dim = dims)
    cursor <- p

    if (model != "normal") {
      output$skew <- array(
        current_eta[cursor + seq_len(p)],
        dim = dims
      )
      cursor <- cursor + p
    }

    reconstructed_sigmas <- vector("list", length(dims))
    for (mode in seq_along(dims)) {
      lower <- lower_chols[[mode]]
      specification <- covariance_specs[[mode]]
      indices <- specification$indices

      if (nrow(indices) > 0L) {
        current_values <- current_eta[
          specification$start:specification$end
        ]
        for (index in seq_len(nrow(indices))) {
          row <- indices[index, 1L]
          column <- indices[index, 2L]
          lower[row, column] <- if (row == column) {
            exp(current_values[index])
          } else {
            current_values[index]
          }
        }
      }

      sigma <- tcrossprod(lower)
      if (mode < length(dims)) {
        sigma <- sigma / mean(diag(sigma))
      }
      reconstructed_sigmas[[mode]] <- sigma
    }
    output$sigmas <- reconstructed_sigmas

    if (model == "skewt") {
      output$nu <- exp(current_eta[shape_start])
    } else if (model == "vargamma") {
      output$gamma <- exp(current_eta[shape_start])
      output$scale <- NULL
    } else if (model == "invgauss") {
      output$kappa <- exp(current_eta[shape_start])
    } else if (model == "genhyper") {
      output$lambda <- current_eta[shape_start]
      output$omega <- exp(current_eta[shape_start + 1L])
    }

    output
  }

  list(
    eta = eta,
    names = parameter_names,
    decode = decode,
    base = base
  )
}

.tensor_score_gof_log_density <- function(draws, model, parameters) {
  values <- switch(
    model,
    normal = dtnorm(
      draws,
      mu = parameters$mu,
      sigmas = parameters$sigmas,
      log = TRUE
    ),
    skewt = dtskewt(
      draws,
      mu = parameters$mu,
      skew = parameters$skew,
      sigmas = parameters$sigmas,
      nu = parameters$nu,
      log = TRUE
    ),
    vargamma = dtvargamma(
      draws,
      mu = parameters$mu,
      skew = parameters$skew,
      sigmas = parameters$sigmas,
      scale = parameters$gamma,
      log = TRUE
    ),
    invgauss = dtinvgauss(
      draws,
      mu = parameters$mu,
      skew = parameters$skew,
      sigmas = parameters$sigmas,
      kappa = parameters$kappa,
      log = TRUE
    ),
    genhyper = dtgenhyper(
      draws,
      mu = parameters$mu,
      skew = parameters$skew,
      sigmas = parameters$sigmas,
      lambda = parameters$lambda,
      omega = parameters$omega,
      log = TRUE
    )
  )

  if (length(values) != n_draws(draws) || any(!is.finite(values))) {
    stop("A perturbed model produced non-finite log densities.",
         call. = FALSE)
  }
  as.numeric(values)
}

.tensor_score_gof_scores <- function(draws,
                                     model,
                                     parameterization,
                                     fd_step) {
  eta <- parameterization$eta
  scores <- matrix(
    NA_real_,
    nrow = n_draws(draws),
    ncol = length(eta),
    dimnames = list(NULL, parameterization$names)
  )

  for (parameter in seq_along(eta)) {
    step <- fd_step * max(1, abs(eta[parameter]))
    succeeded <- FALSE

    for (attempt in seq_len(6L)) {
      eta_plus <- eta
      eta_minus <- eta
      eta_plus[parameter] <- eta_plus[parameter] + step
      eta_minus[parameter] <- eta_minus[parameter] - step

      evaluated <- tryCatch(
        list(
          plus = .tensor_score_gof_log_density(
            draws,
            model,
            parameterization$decode(eta_plus)
          ),
          minus = .tensor_score_gof_log_density(
            draws,
            model,
            parameterization$decode(eta_minus)
          )
        ),
        error = function(error) error
      )

      if (!inherits(evaluated, "error")) {
        scores[, parameter] <-
          (evaluated$plus - evaluated$minus) / (2 * step)
        succeeded <- all(is.finite(scores[, parameter]))
      }
      if (succeeded) break
      step <- step / 2
    }

    if (!succeeded) {
      stop(
        "Could not calculate a finite score for parameter `",
        parameterization$names[parameter],
        "`.",
        call. = FALSE
      )
    }
  }

  scores
}

.tensor_score_gof_legendre <- function(values, degree) {
  transformed <- 2 * values - 1
  output <- vector("list", degree)
  output[[1L]] <- sqrt(3) * transformed

  if (degree == 1L) {
    return(output)
  }

  previous_previous <- matrix(1, nrow(values), ncol(values))
  previous <- transformed

  for (current_degree in 2:degree) {
    current <- (
      (2 * current_degree - 1) * transformed * previous -
        (current_degree - 1) * previous_previous
    ) / current_degree
    output[[current_degree]] <-
      sqrt(2 * current_degree + 1) * current
    previous_previous <- previous
    previous <- current
  }

  output
}

.tensor_score_gof_even_indices <- function(total, count) {
  if (count >= total) return(seq_len(total))
  unique(as.integer(round(seq(1, total, length.out = count))))
}

.tensor_score_gof_basis <- function(uniforms,
                                    marginal_degree,
                                    include_interactions,
                                    max_basis) {
  n <- nrow(uniforms)
  p <- ncol(uniforms)
  legendre <- .tensor_score_gof_legendre(uniforms, marginal_degree)

  marginal_total <- p * marginal_degree
  marginal_count <- min(marginal_total, max_basis)
  marginal_indices <- .tensor_score_gof_even_indices(
    marginal_total,
    marginal_count
  )

  marginal_values <- matrix(NA_real_, nrow = n, ncol = marginal_count)
  marginal_names <- character(marginal_count)
  for (position in seq_along(marginal_indices)) {
    index <- marginal_indices[position]
    degree <- (index - 1L) %/% p + 1L
    coordinate <- (index - 1L) %% p + 1L
    marginal_values[, position] <- legendre[[degree]][, coordinate]
    marginal_names[position] <- paste0(
      "L", degree, "(U", coordinate, ")"
    )
  }

  remaining <- max_basis - marginal_count
  if (!include_interactions || p < 2L || remaining < 1L) {
    return(list(values = marginal_values, names = marginal_names))
  }

  pairs <- utils::combn(p, 2L)
  pair_indices <- .tensor_score_gof_even_indices(
    ncol(pairs),
    min(remaining, ncol(pairs))
  )
  pairs <- pairs[, pair_indices, drop = FALSE]

  interaction_values <- matrix(
    NA_real_,
    nrow = n,
    ncol = ncol(pairs)
  )
  interaction_names <- character(ncol(pairs))
  first_degree <- legendre[[1L]]

  for (pair in seq_len(ncol(pairs))) {
    left <- pairs[1L, pair]
    right <- pairs[2L, pair]
    interaction_values[, pair] <-
      first_degree[, left] * first_degree[, right]
    interaction_names[pair] <- paste0(
      "L1(U", left, ")*L1(U", right, ")"
    )
  }

  list(
    values = cbind(marginal_values, interaction_values),
    names = c(marginal_names, interaction_names)
  )
}

.tensor_score_gof_eigen_inverse <- function(x, eigen_tol) {
  x <- (x + t(x)) / 2
  decomposition <- eigen(x, symmetric = TRUE)
  largest <- max(abs(decomposition$values))

  if (!is.finite(largest) || largest == 0) {
    return(list(
      inverse = matrix(0, nrow(x), ncol(x)),
      rank = 0L,
      values = decomposition$values
    ))
  }

  retained <- decomposition$values > eigen_tol * largest
  inverse <- if (any(retained)) {
    vectors <- decomposition$vectors[, retained, drop = FALSE]
    vectors %*%
      (t(vectors) / decomposition$values[retained])
  } else {
    matrix(0, nrow(x), ncol(x))
  }

  list(
    inverse = inverse,
    rank = sum(retained),
    values = decomposition$values
  )
}
