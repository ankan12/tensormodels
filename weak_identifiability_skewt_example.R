load_tensormodels <- function() {
  if (requireNamespace("tensormodels", quietly = TRUE)) {
    suppressPackageStartupMessages(library(tensormodels))
    return(invisible(TRUE))
  }

  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
    return(invisible(TRUE))
  }

  stop(
    "Could not load tensormodels. Install the package or install pkgload ",
    "to load the local source tree."
  )
}

scalar_loglik <- function(draws, mu, sigma, nu, skew) {
  sum(vapply(
    draws,
    function(x) {
      dtskewt(
        x = x,
        mu = array(mu, dim = 1),
        skew = array(skew, dim = 1),
        sigmas = list(matrix(sigma, nrow = 1, ncol = 1)),
        nu = nu,
        log = TRUE
      )
    },
    numeric(1)
  ))
}

profile_skew <- function(draws, mu, sigma, nu, interval = c(0.2, 2.0)) {
  opt <- optimize(
    f = function(skew) -scalar_loglik(draws, mu = mu, sigma = sigma,
                                      nu = nu, skew = skew),
    interval = interval
  )

  data.frame(
    nu = nu,
    skew = opt$minimum,
    loglik = -opt$objective
  )
}

load_tensormodels()

set.seed(42)

# A 1 x 1 tensor is enough to expose the weak identification.
# We hold mu and sigma fixed at their true values so the ridge is only in
# the (nu, skew) directions.
n_draws <- 50
mu_true <- 0
sigma_true <- 1
nu_true <- 10
skew_true <- 1

draws <- rtskewt(
  n = n_draws,
  mu = array(mu_true, dim = 1),
  sigmas = list(matrix(sigma_true, nrow = 1, ncol = 1)),
  skew = array(skew_true, dim = 1),
  nu = nu_true
)

draw_values <- vapply(draws, as.numeric, numeric(1))

nu_grid <- seq(4.1, 40, length.out = 150)
skew_grid <- seq(0.3, 2.2, length.out = 160)

loglik_surface <- outer(
  nu_grid,
  skew_grid,
  Vectorize(function(nu, skew) {
    scalar_loglik(draws, mu = mu_true, sigma = sigma_true, nu = nu, skew = skew)
  })
)

max_idx <- which(loglik_surface == max(loglik_surface), arr.ind = TRUE)[1, ]
max_nu <- nu_grid[max_idx[1]]
max_skew <- skew_grid[max_idx[2]]
max_loglik <- loglik_surface[max_idx[1], max_idx[2]]

near_half <- which(loglik_surface >= max_loglik - 0.5, arr.ind = TRUE)
near_two <- which(loglik_surface >= max_loglik - 2, arr.ind = TRUE)

profile_nus <- c(10, 15, 20, 30, 40)
profile_table <- do.call(
  rbind,
  lapply(profile_nus, function(nu) {
    profile_skew(draws, mu = mu_true, sigma = sigma_true, nu = nu)
  })
)
profile_table$diff_from_best <- max(profile_table$loglik) - profile_table$loglik

candidate_a <- profile_table[profile_table$nu == 10, ]
candidate_b <- profile_table[profile_table$nu == 40, ]

x_grid <- seq(min(draw_values) - 1, max(draw_values) + 1, length.out = 500)
d_true <- vapply(
  x_grid,
  function(x) {
    dtskewt(x, mu = array(mu_true, dim = 1), skew = array(skew_true, dim = 1),
            sigmas = list(matrix(sigma_true, 1, 1)), nu = nu_true, log = FALSE)
  },
  numeric(1)
)
d_a <- vapply(
  x_grid,
  function(x) {
    dtskewt(x, mu = array(mu_true, dim = 1),
            skew = array(candidate_a$skew, dim = 1),
            sigmas = list(matrix(sigma_true, 1, 1)), nu = candidate_a$nu,
            log = FALSE)
  },
  numeric(1)
)
d_b <- vapply(
  x_grid,
  function(x) {
    dtskewt(x, mu = array(mu_true, dim = 1),
            skew = array(candidate_b$skew, dim = 1),
            sigmas = list(matrix(sigma_true, 1, 1)), nu = candidate_b$nu,
            log = FALSE)
  },
  numeric(1)
)

png(
  filename = "weak_identifiability_skewt_example.png",
  width = 1400,
  height = 500,
  res = 130
)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

contour(
  x = nu_grid,
  y = skew_grid,
  z = loglik_surface,
  nlevels = 12,
  xlab = expression(nu),
  ylab = "skew",
  main = "Log-likelihood ridge"
)
points(max_nu, max_skew, pch = 19, col = "firebrick")
points(profile_table$nu, profile_table$skew, pch = 1, col = "navy")

plot(
  profile_table$nu,
  profile_table$loglik,
  type = "b",
  pch = 19,
  col = "navy",
  xlab = expression(nu),
  ylab = "Profile log-likelihood",
  main = "Nearly flat profile"
)
abline(h = max(profile_table$loglik) - 0.5, lty = 2, col = "firebrick")
text(profile_table$nu, profile_table$loglik, labels = round(profile_table$skew, 3),
     pos = 3, cex = 0.8)

hist(
  draw_values,
  breaks = 12,
  freq = FALSE,
  border = "white",
  col = "grey85",
  xlab = "x",
  main = "Different pairs fit almost identically"
)
lines(x_grid, d_true, col = "black", lwd = 2)
lines(x_grid, d_a, col = "steelblue", lwd = 2, lty = 2)
lines(x_grid, d_b, col = "firebrick", lwd = 2, lty = 3)
legend(
  "topright",
  bty = "n",
  lwd = 2,
  lty = c(1, 2, 3),
  col = c("black", "steelblue", "firebrick"),
  legend = c(
    sprintf("truth: nu = %.0f, skew = %.2f", nu_true, skew_true),
    sprintf("candidate A: nu = %.0f, skew = %.3f", candidate_a$nu, candidate_a$skew),
    sprintf("candidate B: nu = %.0f, skew = %.3f", candidate_b$nu, candidate_b$skew)
  )
)

dev.off()

cat("\nWeak identification example for rtskewt()\n")
cat("----------------------------------------\n")
cat(sprintf("True parameters used to simulate data: nu = %.1f, skew = %.1f\n",
            nu_true, skew_true))
cat(sprintf("Grid maximum: nu = %.2f, skew = %.3f, logLik = %.3f\n",
            max_nu, max_skew, max_loglik))
cat("\nProfile likelihood by fixing nu and optimizing skew:\n")
print(profile_table, row.names = FALSE)

cat(sprintf(
  "\nWithin 0.5 log-likelihood units of the maximum, nu ranges from %.2f to %.2f\n",
  min(nu_grid[near_half[, 1]]), max(nu_grid[near_half[, 1]])
))
cat(sprintf(
  "and skew ranges from %.3f to %.3f.\n",
  min(skew_grid[near_half[, 2]]), max(skew_grid[near_half[, 2]])
))

cat(sprintf(
  "\nTwo very different parameter pairs:\n  A: nu = %.0f, skew = %.3f\n  B: nu = %.0f, skew = %.3f\n",
  candidate_a$nu, candidate_a$skew, candidate_b$nu, candidate_b$skew
))
cat(sprintf(
  "Their profile log-likelihoods differ by only %.3f.\n",
  abs(candidate_a$loglik - candidate_b$loglik)
))

cat(
  "\nInterpretation: once mu and sigma are fixed, the likelihood still has a broad\n",
  "ridge in the (nu, skew) directions. Larger nu can be offset by larger skew,\n",
  "so the data do not pin those two parameters down separately in this sample.\n",
  sep = ""
)
cat("\nSaved figure to weak_identifiability_skewt_example.png\n")
