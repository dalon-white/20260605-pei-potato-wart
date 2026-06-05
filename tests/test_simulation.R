source("R/config.R")
source("R/distributions.R")
source("R/sampling_rules.R")
source("R/upstream_screening.R")
source("R/border_sampling.R")
source("R/metrics.R")
source("R/simulation.R")

cfg <- load_cfg("config/default.yaml")
cfg$shipments_per_year$mean <- 200
cfg$shipments_per_year$dispersion <- 0.1
cfg$prevalence$mu <- 0.1
cfg$prevalence$kappa <- 200
cfg$upstream_screen$f_asym_fixed <- 1
cfg$tiers$N_small_threshold <- 1e9

sim_all <- simulate_year(seed = 123, cfg = cfg, scenario = list(mu = 0.1, kappa = 200, label = "high_p"))
inf_rows <- sim_all$shipment_level$infected_present
stopifnot(all(sim_all$shipment_level$detected[inf_rows]))

cfg_high <- cfg
cfg_high$upstream_screen$f_asym_fixed <- 0
cfg_high$upstream_screen$Se_vis <- 1.0
cfg_low <- cfg
cfg_low$upstream_screen$f_asym_fixed <- 0
cfg_low$upstream_screen$Se_vis <- 0.0

res_high <- simulate_year(seed = 456, cfg = cfg_high, scenario = list(mu = 0.1, kappa = 200, label = "high_vis"))
res_low <- simulate_year(seed = 456, cfg = cfg_low, scenario = list(mu = 0.1, kappa = 200, label = "low_vis"))

risk_high <- res_high$annual_summary$prob_any_infected_pass_undetected[[1]]
risk_low <- res_low$annual_summary$prob_any_infected_pass_undetected[[1]]
stopifnot(risk_high <= risk_low)

message("test_simulation.R passed")
