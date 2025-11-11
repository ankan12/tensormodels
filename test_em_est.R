s1 <- matrix(rnorm(4), nrow = 2); s2 <- matrix(rnorm(9), nrow = 3); s3 <- matrix(rnorm(16), nrow = 4)
s1 <- crossprod(s1, s1); s2 <- crossprod(s2, s2); s3 <- crossprod(s3, s3)
rand_mu <- array(rnorm(24), dim = c(2, 3, 4))
nu <- 20
rand_skew <- array(rnorm(24), dim = c(2, 3, 4))
rtskewt_draws <- rtskewt(n = 1e3, mu = rand_mu, sigmas = list(s1, s2, s3),
                         skew = rand_skew, nu = nu)
exp_mu <- rand_mu + (nu/(nu-2)) * rand_skew
(res <- em_est(rtskewt_draws, quiet = FALSE, model = "skewt"))

# test correctness
tensormodels:::frob_norm_diff(s1/(sum(diag(s1))) * 2, res$sigmas[[1]])
tensormodels:::frob_norm_diff(s2/(sum(diag(s2))) * 3, res$sigmas[[2]])
tensormodels:::frob_norm_diff(s3, res$sigmas[[3]])
tensormodels:::frob_norm_diff(exp_mu, res$mu)
tensormodels:::frob_norm_diff(rand_skew, res$skew)
tensormodels:::frob_norm_diff(nu, res$nu)
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

ggplot() +
  geom_function(fun = update_nu, args = list(b = b, c = c, n = n)) +
  geom_hline(aes(yintercept = 0), lty = 2, color = "red") +
  xlim(1, 100)

theme_minimal()

sim_draws <- rinvgamma(n = 2*3*4, shape = res$nu/2, rate = res$nu/2)

res$mu + res$skew * sim_draws


rand_mu + skew * rinvgamma(n = 2*3*4, shape = nu/2, rate = nu/2)


ggplot() +
  geom_function(fun = dinvgamma, args = list(b = b, c = c, n = n)) +
  xlim(1, 100)
