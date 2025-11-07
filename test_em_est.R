library(tensormodels)

s1 <- matrix(rnorm(4), nrow = 2); s2 <- matrix(rnorm(9), nrow = 3); s3 <- matrix(rnorm(16), nrow = 4)
s1 <- crossprod(s1, s1); s2 <- crossprod(s2, s2); s3 <- crossprod(s3, s3)
rand_mu <- array(rnorm(24), dim = c(2, 3, 4))
nu <- 10
skew <- 4
rtskewt_draws <- rtskewt(n = 1e3, mu = rand_mu, skew = skew, sigmas = list(s1, s2, s3), nu = nu)
exp_mu <- rand_mu + (nu/(nu-2)) * skew
(res <- em_est(rtskewt_draws, quiet = FALSE))

# test correctness
frob_norm_diff(s1/sum(diag(s1)), res$sigmas[[1]])
frob_norm_diff(s2/sum(diag(s2)), res$sigmas[[2]])
frob_norm_diff(s3/sum(diag(s3)), res$sigmas[[3]])
frob_norm_diff(exp_mu, res$mu)
frob_norm_diff(2, res$skew)
frob_norm_diff(12, res$nu)
compare_draws <- rtskewt(n = 1e3, mu = res$mu, skew = res$skew, sigmas = res$sigmas, nu = res$nu)

library(tidyverse)

tibble("orig" = rtskewt_draws |> c(),
       "sim" = compare_draws |> c()) |>
pivot_longer(cols = everything()) |>
filter(value < 1e2 & value > -1e2) |>
ggplot() +
  geom_histogram(aes(x = value), color = "black", fill = "pink") +
  facet_wrap(~name) +
  theme_minimal()

wggplot() +
  geom_function(fun = update_nu, args = list(b = b, c = c, n = n)) +
  xlim(0.2, 1)
