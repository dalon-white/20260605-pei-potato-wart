source("R/config.R")
source("R/sampling_rules.R")

n_bin <- binomial_n(0.001, 0.05)
stopifnot(n_bin == 2995)

n_h_500 <- hypergeom_n(N = 500, p0 = 0.001, alpha = 0.05)
stopifnot(n_h_500 <= 500)

n_h_10 <- hypergeom_n(N = 10, p0 = 0.001, alpha = 0.05)
stopifnot(n_h_10 == 10)

n_inf <- inflate_n_for_sensitivity(p0 = 0.001, alpha = 0.05, Se_eff = 0.5)
stopifnot(n_inf > 2995)

message("test_sampling.R passed")
