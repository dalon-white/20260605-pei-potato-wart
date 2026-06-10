source("R/config.R")
source("R/distributions.R")
source("R/sampling_rules.R")
source("R/upstream_screening.R")
source("R/border_sampling.R")
source("R/metrics.R")
source("R/simulation.R")

cfg <- load_cfg("config/default.yaml")
expected_units <- units_per_shipment(shipment_kg = 1000, unit_oz = cfg$unit_size_oz$mean)
stopifnot(identical(shipment_unit_count(shipment_kg = 1000, unit_cfg = cfg$unit_size_oz), expected_units))

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

res_named_vec <- simulate_year(seed = 789, cfg = cfg, scenario = c(mu = 0.1, kappa = 200, label = "vec"))
stopifnot(is.data.frame(res_named_vec$annual_summary), nrow(res_named_vec$annual_summary) == 1)

res_scalar <- simulate_year(seed = 790, cfg = cfg, scenario = 0.1)
stopifnot(is.data.frame(res_scalar$annual_summary), nrow(res_scalar$annual_summary) == 1)
stopifnot("mean_unit_size_oz" %in% names(res_scalar$shipment_level))
stopifnot(!("unit_size_oz" %in% names(res_scalar$shipment_level)))
stopifnot("year" %in% names(res_scalar$shipment_level))
stopifnot("annual_by_year" %in% names(res_scalar))
stopifnot(is.list(res_scalar$years), length(res_scalar$years) == 1)
stopifnot("year" %in% names(res_scalar$annual_summary))
stopifnot("shipments_mean" %in% names(res_scalar$annual_by_year))

cfg_up <- load_cfg("config/default.yaml")
cfg_up$upstream_screen$f_asym_fixed <- 0
cfg_up$upstream_screen$Se_vis <- 1
up_true <- condition_to_border(N_units = 100, infected_idx = c(1L, 2L, 3L), cfg = cfg_up)
stopifnot(isTRUE(up_true$shipment_filtered_upstream))
stopifnot(!isTRUE(up_true$shipment_reaches_border))
stopifnot(up_true$N_asym == 0L)

cfg_fp <- load_cfg("config/default.yaml")
cfg_fp$upstream_screen$f_asym_fixed <- 0
cfg_fp$upstream_screen$Se_vis <- 0
cfg_fp$upstream_screen$Sp_vis <- 0
up_fp <- condition_to_border(N_units = 100, infected_idx = integer(0), cfg = cfg_fp)
stopifnot(!isTRUE(up_fp$shipment_filtered_upstream))
stopifnot(isTRUE(up_fp$shipment_reaches_border))
stopifnot(up_fp$N_asym == 100L)

sim_cols <- names(res_scalar$shipment_level)
stopifnot("upstream_detected" %in% sim_cols)
stopifnot("border_detected" %in% sim_cols)
stopifnot(all(res_scalar$shipment_level$detected == res_scalar$shipment_level$border_detected))

res_multi <- simulate_year(seed = 901, cfg = cfg, scenario = list(mu = 0.1, kappa = 200, label = "multi"), n_years = 3)
stopifnot(length(res_multi$years) == 3)
stopifnot(is.data.frame(res_multi$annual_summary), nrow(res_multi$annual_summary) == 3)
stopifnot(is.data.frame(res_multi$annual_by_year), nrow(res_multi$annual_by_year) == 1)
stopifnot("years_simulated" %in% names(res_multi$annual_by_year))
stopifnot(res_multi$annual_by_year$years_simulated[[1]] == 3)
stopifnot("shipments_sd" %in% names(res_multi$annual_by_year))
stopifnot(length(unique(res_multi$shipment_level$year)) == 3)

message("test_simulation.R passed")
