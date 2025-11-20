
orig_norm <- rnorm(1e5, mean = 10, sd = 4)

sim_norm <- rnorm(1e5, mean = mean(orig_norm), sd = sd(orig_norm))

diff <- (orig_norm - sim_norm) |> scale()

diff^2 |>
  c() |>
  as_tibble() |>
  ggplot() +
  geom_histogram(aes(x = value, y = after_stat(density)), bins = 100) +
  geom_function(fun = dchisq, args = list(df = 1), color = 'red') +
  xlim(0.08, 12)


s1 <- matrix(rnorm(4), nrow = 2); s2 <- matrix(rnorm(9), nrow = 3); s3 <- matrix(rnorm(16), nrow = 4)
s1 <- crossprod(s1, s1); s2 <- crossprod(s2, s2); s3 <- crossprod(s3, s3)
rand_mu <- array(rnorm(24), dim = c(2, 3, 4))

norm_draws <- rtnorm(n = 1e3, mu = rand_mu, sigmas = list(s1, s2, s3))

est <- mle_est(norm_draws)

est$mu

orig_kroneck <- kronecker(s1, s2) |> kronecker(s3)
sim_kroneck <- kronecker(est$sigmas[[1]], est$sigmas[[2]]) |> kronecker(est$sigmas[[3]])

(orig_kroneck - sim_kroneck) |>
  as_tibble() |>
  mutate(row = row_number()) |>
  pivot_longer(cols = -row) |>
  mutate(name = str_remove(name, "V"), name = as.numeric(name)) |>
ggplot() +
  geom_tile(aes(x = row, y = name, fill = value)) +
  scale_fill_gradient(low = "blue", high = "red")

sum((orig_kroneck - sim_kroneck)^2)

sim_data <- rtnorm(n = 1e3, mu = est$mu, sigmas = est$sigmas)


ggplot() +
  geom_histogram(aes(x = norm_draws - sim_data), color = "black", fill = "pink") +
  theme_minimal()

diff_draws <- (norm_draws - sim_data) |>
  apply(MARGIN = 2:4, FUN = function(x) scale(x, center = FALSE))

diff_draws^2 |>
  as_tibble() |>
  ggplot() +
  geom_histogram(aes(x = V1, y = after_stat(density)), bins = 100) +
  geom_function(fun = dchisq, args = list(df = 1), color = 'red') +
  xlim(0.2, 12)

diff_draws^2 |> sum()

diff_t^2 |> sum()

(diff_draws - diff_t) |>
  as_tibble() |>
  mutate(row = row_number()) |>
  pivot_longer(cols = -row) |>
  mutate(name = str_remove(name, "V"), name = as.numeric(name)) |>
  ggplot() +
  geom_tile(aes(x = row, y = name, fill = value)) +
  scale_fill_gradient()
# sim_data <- rnorm(n = 1e3, mean = 10, sd = sqrt(5)) |> scale()
# more_data <- rnorm(n = 1e3)
#
# sum(sim_data - more_data)
#
# ggplot() +
#   geom_histogram(aes(x = (sim_data - more_data)))

sim_t <- rtskewt(n = 1e3, mu = rand_mu, sigmas = list(s1, s2, s3),
                 skew = 10, nu = 12)

test_res <- mle_est(sim_t)

test_t <- rtnorm(n = 1e3, mu = test_res$mu, sigmas = test_res$sigmas)

diff_t <- (sim_t - test_t) #|> apply(MARGIN = 2:4, FUN = scale)

diff_t <- (sim_t - test_t) |>
  apply(MARGIN = 2:4, FUN = function(x) scale(x, center = FALSE))

diff_t^2 |>
  as_tibble() |>
  ggplot() +
  geom_histogram(aes(x = V1, y = after_stat(density)), bins = 100) +
  geom_function(fun = dchisq, args = list(df = 1), color = 'red') +
  xlim(0.2, 12)

sum(diff_draws^2)

sum(diff_t^2)

(diff_t)^2 |> sum()


#Var =5
#alph
