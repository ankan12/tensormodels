#' Fast score-adjusted Rosenblatt goodness-of-fit test
#'
#' Performs a smooth goodness-of-fit test on tensor Rosenblatt values while
#' accounting for parameters estimated from the same observations. This
#' version uses [tensor_scores()] rather than parameter-by-parameter finite
#' differences.
#'
#' @param draws A `tensor` object containing IID observations.
#' @param model The fitted distribution: `"normal"`, `"skewt"`,
#'   `"vargamma"`, `"invgauss"`, or `"genhyper"`.
#' @param fit An optional parameter list returned by [tensor_mle()]. When
#'   `NULL`, the model is fitted once.
#' @param marginal_degree Positive integer giving the largest marginal
#'   shifted-Legendre degree.
#' @param include_interactions Whether to include pairwise first-degree
#'   interaction terms.
#' @param max_basis Maximum number of basis functions.
#' @param calibration P-value calibration. `"multiplier"` uses a fast
#'   score-adjusted weighted bootstrap and is the default. `"chisq"` uses the
#'   experimental asymptotic chi-squared calibration.
#' @param B Number of weighted-bootstrap replications when
#'   `calibration = "multiplier"`.
#' @param seed Optional random seed for the weighted bootstrap.
#' @param multiplier_chunk_size Number of multiplier replications processed
#'   in one matrix multiplication. This controls temporary memory use without
#'   changing the result.
#' @param eigen_tol Relative eigenvalue threshold used for pseudoinverses and
#'   numerical ranks.
#' @param max_iter,tol,quiet,restrict Arguments passed to [tensor_mle()] when
#'   `fit` is `NULL`. Restricted covariance fits are not currently supported.
#' @param rel.tol,abs.tol,subdivisions Quadrature controls passed to
#'   [tensor_rosenblatt()].
#'
#' @return An object inheriting from `"htest"`. Its `statistic` and `p.value`
#'   entries contain the smooth statistic and its selected calibrated p-value.
#'   The former covariance-whitened statistic and chi-squared p-value are
#'   retained as `asymptotic.statistic` and `asymptotic.p.value`.
#'
#' @details
#' If \eqn{h_i} denotes the selected Rosenblatt basis vector and \eqn{s_i}
#' the fitted-model score, this function constructs the adjusted influence
#' rows
#' \deqn{g_i=h_i-CI^+s_i,}
#' where \eqn{C=Cov(h_i,s_i)} and \eqn{I=Var(s_i)}. For the default weighted
#' bootstrap, the observed statistic is \eqn{n\|\bar h\|^2}. Each null
#' realization is obtained from
#' \deqn{\left\|n^{-1/2}\sum_i(Z_i-\bar Z)g_i\right\|^2,}
#' where the \eqn{Z_i} are independent standard-normal multipliers. This
#' requires no new tensor draws, model fits, or Rosenblatt transformations.
#'
#' The weighted bootstrap is a large-sample calibration and should still be
#' validated for a proposed model and sample-size regime. The chi-squared
#' calibration is retained mainly for comparison.
#'
#' @examples
#' draws <- rtnorm(
#'   100,
#'   mu = array(0, dim = c(2, 2)),
#'   sigmas = list(diag(2), diag(2))
#' )
#' result <- tensor_gof_scores(
#'   draws,
#'   model = "normal",
#'   marginal_degree = 2,
#'   max_basis = 10,
#'   B = 999,
#'   seed = 1
#' )
#' result$statistic
#' result$p.value
#'
#' @export
tensor_gof_scores <- function(
    draws,
    model = c("normal", "skewt", "vargamma", "invgauss", "genhyper"),
    fit = NULL,
    marginal_degree = 2L,
    include_interactions = TRUE,
    max_basis = 100L,
    calibration = c("multiplier", "chisq"),
    B = 1999L,
    seed = NULL,
    multiplier_chunk_size = 250L,
    eigen_tol = 1e-8,
    max_iter = 1000L,
    tol = 1e-6,
    quiet = TRUE,
    restrict = NULL,
    rel.tol = 1e-7,
    abs.tol = 1e-9,
    subdivisions = 100L) {
  model <- match.arg(model)
  calibration <- match.arg(calibration)
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
      "`tensor_gof_scores()`.",
      call. = FALSE
    )
  }
  .tensor_score_gof_validate_controls(
    marginal_degree = marginal_degree,
    include_interactions = include_interactions,
    max_basis = max_basis,
    fd_step = 1e-5,
    eigen_tol = eigen_tol
  )
  if (length(B) != 1L || !is.numeric(B) || is.na(B) ||
      !is.finite(B) || B < 1L || B != floor(B)) {
    stop("`B` must be one positive integer.", call. = FALSE)
  }
  B <- as.integer(B)
  if (length(multiplier_chunk_size) != 1L ||
      !is.numeric(multiplier_chunk_size) ||
      is.na(multiplier_chunk_size) ||
      !is.finite(multiplier_chunk_size) ||
      multiplier_chunk_size < 1L ||
      multiplier_chunk_size != floor(multiplier_chunk_size)) {
    stop(
      "`multiplier_chunk_size` must be one positive integer.",
      call. = FALSE
    )
  }
  multiplier_chunk_size <- as.integer(multiplier_chunk_size)
  if (!is.null(seed) &&
      (length(seed) != 1L || !is.numeric(seed) || is.na(seed) ||
       !is.finite(seed))) {
    stop("`seed` must be `NULL` or one finite number.", call. = FALSE)
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
  scores <- tensor_scores(
    draws = draws,
    model = model,
    fit = fitted_parameters
  )

  score_dimension <- ncol(scores)
  if (n <= score_dimension) {
    warning(
      "The sample size (", n, ") is not larger than the score dimension (",
      score_dimension, "); the chi-squared calibration is not reliable.",
      call. = FALSE
    )
  }

  centered_scores <- sweep(scores, 2L, colMeans(scores), "-")
  centered_basis <- sweep(
    basis$values,
    2L,
    colMeans(basis$values),
    "-"
  )
  information <- crossprod(centered_scores) / n
  information_inverse <- .tensor_score_gof_eigen_inverse(
    information,
    eigen_tol
  )
  if (information_inverse$rank < score_dimension) {
    warning(
      "The nuisance-score information matrix is rank deficient (rank ",
      information_inverse$rank, " of ", score_dimension,
      "); interpret the chi-squared p-value cautiously.",
      call. = FALSE
    )
  }

  cross_covariance <- crossprod(centered_basis, centered_scores) / n
  projection <- cross_covariance %*% information_inverse$inverse
  adjusted_influence <-
    centered_basis - centered_scores %*% t(projection)
  adjusted_covariance <- crossprod(adjusted_influence) / n
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

  basis_mean <- colMeans(basis$values)
  raw_statistic <- n * sum(basis_mean^2)
  asymptotic_statistic <- as.numeric(
    n * crossprod(
      basis_mean,
      covariance_inverse$inverse %*% basis_mean
    )
  )
  degrees_freedom <- covariance_inverse$rank
  asymptotic_p_value <- stats::pchisq(
    asymptotic_statistic,
    df = degrees_freedom,
    lower.tail = FALSE
  )

  multiplier_statistics <- NULL
  if (calibration == "multiplier") {
    if (!is.null(seed)) set.seed(seed)
    multiplier_statistics <- numeric(B)
    starts <- seq.int(1L, B, by = multiplier_chunk_size)

    for (start in starts) {
      end <- min(B, start + multiplier_chunk_size - 1L)
      current_size <- end - start + 1L
      multipliers <- t(matrix(
        stats::rnorm(current_size * n),
        nrow = n,
        ncol = current_size
      ))
      multiplier_process <-
        multipliers %*% adjusted_influence / sqrt(n)
      multiplier_statistics[start:end] <-
        rowSums(multiplier_process^2)
    }

    statistic <- raw_statistic
    p_value <-
      (1 + sum(multiplier_statistics >= statistic)) / (B + 1)
  } else {
    statistic <- asymptotic_statistic
    p_value <- asymptotic_p_value
  }

  score_scales <- sqrt(pmax(diag(information), 0) / n)
  standardized_score_mean <- rep(0, score_dimension)
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
        "Analytic-score-adjusted Rosenblatt smooth GOF test for the ",
        model_label,
        " with ",
        if (calibration == "multiplier") {
          "weighted-multiplier calibration"
        } else {
          "asymptotic chi-squared calibration"
        }
      ),
      data.name = data_name,
      model = model,
      fit = fitted_parameters,
      uniforms = uniforms,
      basis_names = basis$names,
      basis_means = basis_mean,
      calibration = calibration,
      B = if (calibration == "multiplier") B else NULL,
      multiplier_statistics = multiplier_statistics,
      adjusted_covariance = adjusted_covariance,
      adjusted_influence = adjusted_influence,
      asymptotic.statistic = c(T = asymptotic_statistic),
      asymptotic.p.value = asymptotic_p_value,
      score_dimension = score_dimension,
      information_rank = information_inverse$rank,
      adjusted_covariance_rank = covariance_inverse$rank,
      score_stationarity = stationarity,
      scores = scores,
      parameter_names = colnames(scores)
    ),
    class = c("tensor_gof_scores_test", "htest")
  )
}
