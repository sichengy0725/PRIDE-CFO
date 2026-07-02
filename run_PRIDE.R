


library(Iso)
library(coda)

K <- 5


# choose increasing prior means c_k, then mu=logit(c_k)
c_k <- c(0.1, 0.2, 0.3, 0.4, 0.5)
mu  <- qlogis(c_k)
sce1 <- c(0.33, 0.45, 0.58, 0.70, 0.80)
sce2 <- c(0.18, 0.33, 0.52, 0.60, 0.70)
sce3 <- c(0.11, 0.22, 0.33, 0.40, 0.50)
sce4 <- c(0.01, 0.02, 0.03, 0.33, 0.50)
sce5 <- c(0,0,0.05,0.10,0.33)
sce6 <- c(0.45, 0.55, 0.65, 0.75, 0.85)

out1 <- simulate_PRIDE_design(
  N_pat = 30,
  Nmax_eff = 24,
  C = 3,
  T_assess = 3,
  cycle_max = 3,
  arrival_rate = 1/2,
  p_true = sce1,
  sigma_true_w = 0.4,
  K = K,
  mu = mu,
  TARGET = 0.33,
  cutoff = 0.95,
  model_file = "PRIDE.bug",
  seed = 1,
  verbose = TRUE,
  CFO = FALSE
)

out2 <- simulate_PRIDE_design(
  N_pat = 30,
  Nmax_eff = 30,
  C = 3,
  T_assess = 3,
  cycle_max = 3,
  arrival_rate = 1/2,
  p_true = sce2,
  sigma_true_w = 0.4,
  K = K,
  mu = mu,
  TARGET = 0.30,
  cutoff = 0.95,
  model_file = "PRIDE.bug",
  seed = job_i,
  verbose = TRUE
)
out3 <- simulate_PRIDE_design(
  N_pat = 30,
  Nmax_eff = 30,
  C = 3,
  T_assess = 3,
  cycle_max = 3,
  arrival_rate = 1/2,
  p_true = sce3,
  sigma_true_w = 0.4,
  K = K,
  mu = mu,
  TARGET = 0.30,
  cutoff = 0.95,
  model_file = "PRIDE.bug",
  seed = job_i,
  verbose = TRUE
)
out4 <- simulate_PRIDE_design(
  N_pat = 30,
  Nmax_eff = 30,
  C = 3,
  T_assess = 3,
  cycle_max = 3,
  arrival_rate = 1/2,
  p_true = sce4,
  sigma_true_w = 0.4,
  K = K,
  mu = mu,
  TARGET = 0.30,
  cutoff = 0.95,
  model_file = "PRIDE.bug",
  seed = job_i,
  verbose = TRUE
)
out5 <- simulate_PRIDE_design(
  N_pat = 30,
  Nmax_eff = 30,
  C = 3,
  T_assess = 3,
  cycle_max = 3,
  arrival_rate = 1/2,
  p_true = sce5,
  sigma_true_w = 0.4,
  K = K,
  mu = mu,
  TARGET = 0.30,
  cutoff = 0.95,
  model_file = "PRIDE.bug",
  seed = 1,
  verbose = TRUE
)
out6 <- simulate_PRIDE_design(
  N_pat = 30,
  Nmax_eff = 30,
  C = 3,
  T_assess = 3,
  cycle_max = 3,
  arrival_rate = 1/2,
  p_true = sce6,
  sigma_true_w = 0.4,
  K = K,
  mu = mu,
  TARGET = 0.30,
  cutoff = 0.95,
  model_file = "PRIDE.bug",
  seed = 1,
  verbose = TRUE
)
