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

# What it does: Chooses dispersion parameter kappa for each mu scenario, with fallback to a default.
# Why it matters: Kappa controls heterogeneity/overdispersion in infection prevalence, which strongly affects undetected-pass risk.
# Limitation: It picks the nearest mu scenario rather than interpolating; this can create step changes and potentially misstate risk between defined points.
pick_kappa <- function(cfg, mu) {
  if (is.null(cfg$prevalence$scenarios) || length(cfg$prevalence$scenarios) == 0) {
    return(cfg$prevalence$kappa)
  }

  get_field_num <- function(x, field) {
    if (is.list(x) && !is.null(x[[field]])) {
      return(as.numeric(x[[field]]))
    }
    if (is.atomic(x)) {
      if (!is.null(names(x)) && field %in% names(x)) {
        return(as.numeric(x[[field]]))
      }
      if (field == "mu" && length(x) == 1) {
        return(as.numeric(x[[1]]))
      }
    }
    NA_real_
  }

  mus <- vapply(cfg$prevalence$scenarios, function(x) get_field_num(x, "mu"), numeric(1))
  if (all(is.na(mus))) {
    return(cfg$prevalence$kappa)
  }
  idx <- which.min(abs(mus - mu))

  kappa_i <- get_field_num(cfg$prevalence$scenarios[[idx]], "kappa")
  if (is.na(kappa_i)) cfg$prevalence$kappa else kappa_i
}

build_runtime_scenarios <- function(cfg, override_mus = NULL) {
  get_field_num <- function(x, field) {
    if (is.list(x) && !is.null(x[[field]])) {
      return(as.numeric(x[[field]]))
    }
    if (is.atomic(x) && !is.null(names(x)) && field %in% names(x)) {
      return(as.numeric(x[[field]]))
    }
    NA_real_
  }

  get_field_chr <- function(x, field) {
    if (is.list(x) && !is.null(x[[field]])) {
      return(as.character(x[[field]]))
    }
    if (is.atomic(x) && !is.null(names(x)) && field %in% names(x)) {
      return(as.character(x[[field]]))
    }
    NA_character_
  }

  if (!is.null(override_mus)) {
    override_mus <- as.numeric(override_mus)
    if (any(!is.finite(override_mus))) {
      stop("`--scenarios` values must be numeric and finite.")
    }

    return(lapply(seq_along(override_mus), function(i) {
      mu_i <- override_mus[[i]]
      list(mu = mu_i, kappa = pick_kappa(cfg, mu_i), label = sprintf("p=%.6g", mu_i))
    }))
  }

  scn_cfg <- cfg$prevalence$scenarios
  if (is.null(scn_cfg) || length(scn_cfg) == 0) {
    return(list(list(
      mu = as.numeric(cfg$prevalence$mu),
      kappa = as.numeric(cfg$prevalence$kappa),
      label = sprintf("mu=%.6g_kappa=%.6g", cfg$prevalence$mu, cfg$prevalence$kappa)
    )))
  }

  out <- lapply(seq_along(scn_cfg), function(i) {
    scn_i <- scn_cfg[[i]]
    mu_i <- get_field_num(scn_i, "mu")
    if (!is.finite(mu_i)) {
      stop(sprintf("prevalence.scenarios[[%s]] is missing a valid `mu`.", i))
    }

    kappa_i <- get_field_num(scn_i, "kappa")
    if (!is.finite(kappa_i)) {
      kappa_i <- as.numeric(cfg$prevalence$kappa)
    }

    label_i <- get_field_chr(scn_i, "label")
    if (!is.finite(nchar(label_i)) || nchar(label_i) == 0) {
      label_i <- sprintf("p=%.6g", mu_i)
    }

    list(mu = as.numeric(mu_i), kappa = as.numeric(kappa_i), label = label_i)
  })

  out
}


# Load and finalize runtime configuration:
args <- parse_args(commandArgs(trailingOnly = TRUE))
cfg <- load_cfg(args$config)
if (!is.null(args$seed)) cfg$random_seed <- args$seed

# Define output directories and ensure they exist:
out_dir <- args$outputs_dir
csv_dir <- file.path(out_dir, "csv")
fig_dir <- file.path(out_dir, "figures")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Build runtime scenarios directly from prevalence.scenarios unless overridden by --scenarios:
runtime_scenarios <- build_runtime_scenarios(cfg, override_mus = args$scenarios)
runtime_mus <- vapply(runtime_scenarios, function(s) as.numeric(s$mu), numeric(1))

all_shipments <- list()
all_annual <- list()


# What it does: For each prevalence p, builds scenario with chosen kappa, runs simulate_year, and stores shipment-level and annual summaries.
# Why it matters: This is the key estimator: it generates simulated outcomes under the protocol and quantifies how often infected shipments pass undetected.
# Limitation: Apparent single run per scenario seed sequence; unless simulate_year internally performs enough Monte Carlo replication, uncertainty may be under-characterized
for (i in seq_along(runtime_scenarios)) {
  sim <- simulate_year(seed = cfg$random_seed + i - 1, cfg = cfg, scenario = runtime_scenarios[[i]])
  all_shipments[[i]] <- sim$shipment_level
  all_annual[[i]] <- sim$annual_summary
}

# What it does: Combines scenario outputs; computes required sampling table and OC curve table based on effective sensitivity.
# Why it matters: Translates raw simulations into operational decision support (how much to sample, and expected detection performance).
# Limitation: OC curve n-grid increments by 25 may be coarse; may miss fine threshold behavior near decision boundaries.
shipment_level <- dplyr::bind_rows(all_shipments)
annual_summary <- dplyr::bind_rows(all_annual)
required_tbl <- required_n_table(shipment_level)
upstream_tbl <- upstream_filtering_table(shipment_level, shipments_possible_basis = cfg$shipments_per_year$mean)
oc_tbl <- oc_curve_table(p_values = runtime_mus, n_values = seq(0, max(shipment_level$n_required), by = 25), Se_eff = Se_eff(cfg))

# What it does: Perturbs key parameters (visual sensitivity, molecular sensitivity, subsample sensitivity, asymptomatic fraction, kappa), reruns simulation, and records annual undetected risk metric.
# Why it matters: Shows robustness of protection likelihood to uncertain assumptions; essential for policy confidence.
# Limitation: Univariate sensitivity only; no joint uncertainty, interactions, or probabilistic calibration across parameters.
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

# Conditionally write tabular results and markdown report
if (isTRUE(cfg$outputs$save_csv)) {
  readr::write_csv(shipment_level, file.path(csv_dir, "shipment_level.csv"))
  readr::write_csv(annual_summary, file.path(csv_dir, "annual_summary.csv"))
  readr::write_csv(upstream_tbl, file.path(csv_dir, "upstream_filtering_summary.csv"))
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
    upstream_tbl = upstream_tbl,
    output_path = file.path(out_dir, "report.md")
  )
}

message("Simulation complete. Outputs in: ", normalizePath(out_dir))


# The script estimates likelihood of protection indirectly through detection failure risk. In practice, you interpret protection as:
# Lower prob_any_infected_pass_undetected implies higher protection likelihood for Commodity B.
#    - This is the key metric for decision-makers: it quantifies the risk that an infected shipment slips through the protocol undetected, which directly relates to biosecurity risk.
# Required sample size tables and OC curves inform operational decisions on how much sampling is needed to achieve certain risk thresholds, given the effective sensitivity of the testing protocol.
# Sensitivity analyses show how robust the protection likelihood is to key assumptions, which is crucial for confidence in the protocol under real-world uncertainty. 

# So the script gets closer to the question by:
# 1. Simulating protocol performance under plausible prevalence conditions.
# 2. Quantifying undetected infection risk annually and by sensitivity assumptions.
# 3. Delivering tables/plots/report to compare scenarios and stress-test assumptions.

# Main methodological limits for the exact question:
# 1. It quantifies biosecurity risk, not direct Commodity B economic outcomes (market access, losses avoided, etc.).
# 2. It does not visibly output confidence intervals in this script layer.
# 3. Sensitivity is deterministic and one-factor-at-a-time, not full probabilistic uncertainty propagation.
# 4. Validity depends heavily on assumptions in config and sourced model functions.