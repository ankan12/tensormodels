#' plot_density
#'
#' Plots two-dimensional density slices for a 2 by 2 matrix-variate model.
#'
#' @param x12 A numeric vector of fixed values for `X[1, 2]`.
#' @param x22 A numeric vector of fixed values for `X[2, 2]`.
#' @param x11 A numeric vector of values to evaluate for `X[1, 1]`.
#' @param x21 A numeric vector of values to evaluate for `X[2, 1]`.
#' @param mu A 2 by 2 location matrix.
#' @param sigmas A list of covariance matrices.
#' @param model The model name. One of `"normal"`, `"skewt"`, `"vargamma"`,
#'   `"invgauss"`, or `"genhyper"`.
#' @param skew A 2 by 2 skewness matrix. Required for the skewed models.
#' @param nu Degrees of freedom. Required when `model = "skewt"`.
#' @param scale Scale parameter. Required when `model = "vargamma"`.
#' @param kappa Kappa parameter. Required when `model = "invgauss"`.
#' @param lambda Lambda parameter. Required when `model = "genhyper"`.
#' @param omega Omega parameter. Required when `model = "genhyper"`.
#' @param contour_bins Number of contour bins to draw.
#' @param interpolate Logical. If TRUE, interpolate the raster display.
#'
#' @return A ggplot object.
#'
#' @examples
#' mu <- matrix(c(0, 0, 0, 0), nrow = 2)
#' sigmas <- list(diag(2), diag(2))
#' plot_density(
#'   x12 = seq(-2, 2, by = 1),
#'   x22 = seq(-2, 2, by = 1),
#'   x11 = seq(-4, 4, length.out = 25),
#'   x21 = seq(-4, 4, length.out = 25),
#'   mu = mu,
#'   sigmas = sigmas,
#'   model = "normal"
#' )
#'
#' @export
plot_density <- function(x12, x22, x11, x21, mu, sigmas,
                         model = c("normal", "skewt", "vargamma",
                                   "invgauss", "genhyper"),
                         skew = NULL, nu = NULL, scale = NULL,
                         kappa = NULL, lambda = NULL, omega = NULL,
                         contour_bins = NULL, interpolate = TRUE) {
  model <- match.arg(model)

  if (!is.matrix(mu) || !identical(dim(mu), c(2L, 2L))) {
    stop("mu must be a 2 by 2 matrix.", call. = FALSE)
  }

  x12 <- .validate_density_axis(x12, "x12")
  x22 <- .validate_density_axis(x22, "x22")
  x11 <- .validate_density_axis(x11, "x11")
  x21 <- .validate_density_axis(x21, "x21")

  eval_density <- .make_plot_density_evaluator(
    model = model,
    mu = mu,
    sigmas = sigmas,
    skew = skew,
    nu = nu,
    scale = scale,
    kappa = kappa,
    lambda = lambda,
    omega = omega
  )

  density_grid <- expand.grid(
    x11 = x11,
    x21 = x21,
    x12 = x12,
    x22 = x22,
    KEEP.OUT.ATTRS = FALSE
  )

  x_list <- lapply(seq_len(nrow(density_grid)), function(i) {
    matrix(
      c(
        density_grid$x11[i],
        density_grid$x21[i],
        density_grid$x12[i],
        density_grid$x22[i]
      ),
      nrow = 2
    )
  })

  density_grid$log_density <- eval_density(x_list, log = TRUE)
  density_grid$x12 <- factor(density_grid$x12, levels = x12)
  density_grid$x22 <- factor(density_grid$x22, levels = x22)

  ggplot2::ggplot(density_grid, ggplot2::aes(x = x11, y = x21)) +
    ggplot2::geom_raster(
      ggplot2::aes(fill = log_density),
      interpolate = interpolate
    ) +
    ggplot2::geom_contour(
      ggplot2::aes(z = log_density),
      color = "black",
      linewidth = 0.4,
      bins = contour_bins
    ) +
    ggplot2::facet_grid(rows = ggplot2::vars(x22), cols = ggplot2::vars(x12)) +
    ggplot2::scale_fill_viridis_c(option = "turbo") +
    ggplot2::coord_fixed(expand = FALSE) +
    ggplot2::labs(
      x = expression(X[11]),
      y = expression(X[21]),
      fill = "log density",
      title = paste0(.density_model_title(model), " matrix-variate density slices"),
      subtitle = expression(
        "Columns fix " * X[12] * "; rows fix " * X[22] *
          ". Each panel varies " * X[11] * " and " * X[21] * "."
      )
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(size = 7),
      strip.text.y = ggplot2::element_text(size = 7, angle = 0),
      panel.spacing = grid::unit(0.05, "lines"),
      plot.title = ggplot2::element_text(size = 18, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 12),
      legend.position = "right"
    )
}

.validate_density_axis <- function(x, name) {
  if (!is.numeric(x) || length(x) == 0L || any(!is.finite(x))) {
    stop(name, " must be a non-empty numeric vector with finite values.",
         call. = FALSE)
  }

  x
}

.make_plot_density_evaluator <- function(model, mu, sigmas, skew = NULL,
                                         nu = NULL, scale = NULL,
                                         kappa = NULL, lambda = NULL,
                                         omega = NULL) {
  if (model == "normal") {
    return(.make_dtnorm_evaluator(mu = mu, sigmas = sigmas, log = TRUE))
  }

  if (is.null(skew)) {
    stop("skew must be supplied for skewed models.", call. = FALSE)
  }

  if (!is.matrix(skew) || !identical(dim(skew), c(2L, 2L))) {
    stop("skew must be a 2 by 2 matrix.", call. = FALSE)
  }

  switch(
    model,
    skewt = {
      if (is.null(nu)) stop("nu must be supplied for model = 'skewt'.",
                            call. = FALSE)
      .make_dtskewt_evaluator(
        mu = mu, skew = skew, sigmas = sigmas, nu = nu, log = TRUE
      )
    },
    vargamma = {
      if (is.null(scale)) stop("scale must be supplied for model = 'vargamma'.",
                               call. = FALSE)
      .make_dtvargamma_evaluator(
        mu = mu, skew = skew, sigmas = sigmas, scale = scale, log = TRUE
      )
    },
    invgauss = {
      if (is.null(kappa)) stop("kappa must be supplied for model = 'invgauss'.",
                               call. = FALSE)
      .make_dtinvgauss_evaluator(
        mu = mu, skew = skew, sigmas = sigmas, kappa = kappa, log = TRUE
      )
    },
    genhyper = {
      if (is.null(lambda) || is.null(omega)) {
        stop("lambda and omega must be supplied for model = 'genhyper'.",
             call. = FALSE)
      }
      .make_dtgenhyper_evaluator(
        mu = mu, skew = skew, sigmas = sigmas, lambda = lambda,
        omega = omega, log = TRUE
      )
    }
  )
}

.density_model_title <- function(model) {
  switch(
    model,
    normal = "Normal",
    skewt = "Skew-t",
    vargamma = "Variance-gamma",
    invgauss = "Inverse Gaussian",
    genhyper = "Generalized hyperbolic"
  )
}
