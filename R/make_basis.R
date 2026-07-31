#' Construct smooth basis functions from Rosenblatt values
#'
#' Converts uniform Rosenblatt values into an orthonormal shifted-Legendre
#' basis for smooth goodness-of-fit tests. Marginal terms detect departures
#' from univariate uniformity, while optional pairwise products detect
#' dependence between Rosenblatt coordinates.
#'
#' @param uniforms Uniform values returned by [tensor_rosenblatt()], or a
#'   numeric vector, matrix, or array. Matrix rows are observations and columns
#'   are Rosenblatt coordinates. For an array, the first dimension indexes
#'   observations and all remaining dimensions are flattened into coordinates.
#' @param marginal_degree Positive integer giving the highest shifted-Legendre
#'   degree included separately for each coordinate.
#' @param include_interactions Whether to include pairwise products of the
#'   first-degree marginal terms.
#' @param max_basis Maximum number of basis functions. Marginal terms are
#'   retained first. If space remains, a deterministic spread of pairwise
#'   interaction terms is retained.
#'
#' @return A numeric matrix with one row per observation and one column per
#'   basis function. Column names identify the terms.
#'
#' @details
#' If \eqn{U} is uniform on \eqn{(0,1)}, every marginal basis term has
#' theoretical mean zero and variance one. Distinct terms are orthogonal under
#' independent uniforms. The first two terms for one coordinate are
#' \deqn{L_1(U)=\sqrt{3}(2U-1)}
#' and
#' \deqn{L_2(U)=\sqrt{5}\{6U^2-6U+1\}.}
#'
#' The interaction terms currently use products
#' \eqn{L_1(U_j)L_1(U_k)}. Consequently, they target pairwise dependence but
#' are not an exhaustive basis for every possible multivariate alternative.
#'
#' @examples
#' set.seed(1)
#' uniforms <- matrix(stats::runif(200), nrow = 100, ncol = 2)
#' basis <- make_basis(
#'   uniforms,
#'   marginal_degree = 2,
#'   include_interactions = TRUE
#' )
#' colMeans(basis)
#'
#' @export
make_basis <- function(uniforms,
                       marginal_degree = 2L,
                       include_interactions = TRUE,
                       max_basis = 100L) {
  .tensor_score_gof_validate_controls(
    marginal_degree = marginal_degree,
    include_interactions = include_interactions,
    max_basis = max_basis,
    fd_step = 1e-5,
    eigen_tol = 1e-8
  )

  if (inherits(uniforms, "tensor")) {
    number_observations <- n_draws(uniforms)
    coordinate_matrix <- matrix(
      as.numeric(unclass(uniforms)),
      nrow = number_observations
    )
    rownames(coordinate_matrix) <- dimnames(unclass(uniforms))[[1L]]
  } else if (is.numeric(uniforms)) {
    input_dimensions <- dim(uniforms)

    if (is.null(input_dimensions)) {
      coordinate_matrix <- matrix(
        as.numeric(uniforms),
        ncol = 1L,
        dimnames = list(names(uniforms), "U1")
      )
    } else {
      coordinate_matrix <- matrix(
        as.numeric(uniforms),
        nrow = input_dimensions[1L]
      )
      rownames(coordinate_matrix) <- dimnames(uniforms)[[1L]]
    }
  } else {
    stop(
      "`uniforms` must be a numeric vector, matrix, array, or `tensor`.",
      call. = FALSE
    )
  }

  if (nrow(coordinate_matrix) < 1L || ncol(coordinate_matrix) < 1L) {
    stop("`uniforms` must contain at least one value.", call. = FALSE)
  }
  if (any(!is.finite(coordinate_matrix)) ||
      any(coordinate_matrix < 0 | coordinate_matrix > 1)) {
    stop(
      "Every value in `uniforms` must be finite and between zero and one.",
      call. = FALSE
    )
  }

  basis <- .tensor_score_gof_basis(
    uniforms = coordinate_matrix,
    marginal_degree = as.integer(marginal_degree),
    include_interactions = include_interactions,
    max_basis = as.integer(max_basis)
  )
  colnames(basis$values) <- basis$names
  rownames(basis$values) <- rownames(coordinate_matrix)

  basis$values
}
