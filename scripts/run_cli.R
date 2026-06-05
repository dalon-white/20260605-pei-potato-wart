#!/usr/bin/env Rscript

source("R/config.R")
source("R/distributions.R")
source("R/sampling_rules.R")
source("R/upstream_screening.R")
source("R/border_sampling.R")
source("R/metrics.R")
source("R/plotting.R")
source("R/simulation.R")
source("R/report.R")

required_pkgs <- c("yaml", "dplyr", "tidyr", "purrr", "ggplot2", "stringr", "tibble", "readr", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

parse_args <- function(args) {
  parsed <- list(
    config = "config/default.yaml",
    seed = NULL,
    outputs_dir = "outputs",
    scenarios = NULL
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    val <- if (i < length(args)) args[[i + 1]] else NULL
    if (key == "--config") parsed$config <- val
    if (key == "--seed") parsed$seed <- as.integer(val)
    if (key == "--outputs_dir") parsed$outputs_dir <- val
    if (key == "--scenarios") parsed$scenarios <- as.numeric(strsplit(val, ",")[[1]])
    i <- i + 2
  }
  parsed
}

pick_kappa <- function(cfg, mu) {
  if (is.null(cfg$prevalence$scenarios) || length(cfg$prevalence$scenarios) == 0) {
    return(cfg$prevalence$kappa)
  }
  mus <- vapply(cfg$prevalence$scenarios, function(x) x$mu, numeric(1))
  idx <- which.min(abs(mus - mu))
  cfg$prevalence$scenarios[[idx]]$kappa
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
cfg <- load_cfg(args$config)
if (!is.null(args$seed)) cfg$random_seed <- args$seed

out_dir <- args$outputs_dir
csv_dir <- file.path(out_dir, "csv")
fig_dir <- file.path(out_dir, "figures")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

scenario_ps <- args$scenarios %||% cfg$scenario_p

all_shipments <- list()
all_annual <- list()

for (i in seq_along(scenario_ps)) {
  p <- scenario_ps[[i]]
  scenario <- list(mu = p, kappa = pick_kappa(cfg, p), label = sprintf("p=%.4f", p))
  sim <- simulate_year(seed = cfg$random_seed + i - 1, cfg = cfg, scenario = scenario)
  all_shipments[[i]] <- sim$shipment_level
  all_annual[[i]] <- sim$annual_summary
}

shipment_level <- dplyr::bind_rows(all_shipments)
annual_summary <- dplyr::bind_rows(all_annual)
required_tbl <- required_n_table(shipment_level)
oc_tbl <- oc_curve_table(p_values = scenario_ps, n_values = seq(0, max(shipment_level$n_required), by = 25), Se_eff = Se_eff(cfg))

sensitivity_runs <- list(
  list(parameter = "Se_vis", values = c(0.6, 0.8, 0.95)),
  list(parameter = "Se_mol", values = c(0.7, 0.85, 1.0)),
  list(parameter = "Se_subsample", values = c(0.7, 0.85, 1.0)),
  list(parameter = "f_asym_fixed", values = c(0.05, 0.15, 0.25)),
  list(parameter = "kappa", values = c(200, 500, 1000))
)

sensitivity_tbl <- purrr::map_dfr(sensitivity_runs, function(run) {
  purrr::map_dfr(run$values, function(v) {
    cfg_i <- cfg
    if (run$parameter == "Se_vis") cfg_i$upstream_screen$Se_vis <- v
    if (run$parameter == "Se_mol") cfg_i$border_testing$Se_mol <- v
    if (run$parameter == "Se_subsample") cfg_i$border_testing$Se_subsample <- v
    if (run$parameter == "f_asym_fixed") cfg_i$upstream_screen$f_asym_fixed <- v
    if (run$parameter == "kappa") cfg_i$prevalence$kappa <- v

    sim_i <- simulate_year(
      seed = cfg_i$random_seed + as.integer(abs(v) * 100),
      cfg = cfg_i,
      scenario = list(mu = cfg_i$prevalence$mu, kappa = cfg_i$prevalence$kappa, label = "baseline")
    )

    ann <- sim_i$annual_summary
    tibble::tibble(
      parameter = run$parameter,
      value = v,
      prob_any_infected_pass_undetected = ann$prob_any_infected_pass_undetected[[1]]
    )
  })
})

if (isTRUE(cfg$outputs$save_csv)) {
  readr::write_csv(shipment_level, file.path(csv_dir, "shipment_level.csv"))
  readr::write_csv(annual_summary, file.path(csv_dir, "annual_summary.csv"))
}

if (isTRUE(cfg$outputs$save_plots)) {
  plot_required_n_vs_N(required_tbl, file.path(fig_dir, "required_n_vs_N.png"))
  plot_oc_curve(oc_tbl, file.path(fig_dir, "oc_curve.png"))
  plot_annual_risk_sensitivity(sensitivity_tbl, file.path(fig_dir, "annual_risk_sensitivity.png"))
}

if (isTRUE(cfg$outputs$save_report)) {
  write_report(
    cfg = cfg,
    annual_summary = annual_summary,
    required_tbl = required_tbl,
    oc_tbl = oc_tbl,
    sensitivity_tbl = sensitivity_tbl,
    output_path = file.path(out_dir, "report.md")
  )
}

message("Simulation complete. Outputs in: ", normalizePath(out_dir))
