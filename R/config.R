#' Configuration loading and validation
#'
#' @description Load YAML configuration, apply defaults, validate key fields,
#' and expose deterministic seed initialization.

#' Null-coalescing helper.
#' @param x A value.
#' @param y Fallback value.
#' @return \\code{x} when not NULL, otherwise \\code{y}.
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Recursively merge two lists.
#' @param base Base list.
#' @param override Override list.
#' @return Merged list.
merge_lists <- function(base, override) {
  if (is.null(override)) {
    return(base)
  }
  out <- base
  for (nm in names(override)) {
    if (is.list(out[[nm]]) && is.list(override[[nm]])) {
      out[[nm]] <- merge_lists(out[[nm]], override[[nm]])
    } else {
      out[[nm]] <- override[[nm]]
    }
  }
  out
}

#' Validate configuration values.
#' @param cfg Configuration list.
#' @return Validated configuration.
validate_cfg <- function(cfg) {
  stopifnot(is.numeric(cfg$random_seed), length(cfg$random_seed) == 1)
  stopifnot(is.numeric(cfg$simulation_years), length(cfg$simulation_years) == 1, cfg$simulation_years >= 1)
  stopifnot(cfg$alpha > 0, cfg$alpha < 1)
  stopifnot(cfg$p0 > 0, cfg$p0 < 1)
  if (!is.null(cfg$scenario_p)) {
    stopifnot(is.numeric(cfg$scenario_p), all(cfg$scenario_p > 0), all(cfg$scenario_p < 1))
  }

  stopifnot(cfg$shipments_per_year$mean > 0)
  stopifnot(cfg$shipments_per_year$dispersion > 0)

  stopifnot(cfg$shipment_mass_kg$min > 0, cfg$shipment_mass_kg$max >= cfg$shipment_mass_kg$min)
  stopifnot(cfg$unit_size_oz$min > 0, cfg$unit_size_oz$max >= cfg$unit_size_oz$min)

  stopifnot(cfg$prevalence$mu > 0, cfg$prevalence$mu < 1)
  stopifnot(cfg$prevalence$kappa > 0)
  if (!is.null(cfg$prevalence$scenarios) && length(cfg$prevalence$scenarios) > 0) {
    for (i in seq_along(cfg$prevalence$scenarios)) {
      sc <- cfg$prevalence$scenarios[[i]]
      mu_i <- as.numeric(sc$mu %||% NA_real_)
      kappa_i <- as.numeric(sc$kappa %||% cfg$prevalence$kappa)
      stopifnot(is.finite(mu_i), mu_i > 0, mu_i < 1)
      stopifnot(is.finite(kappa_i), kappa_i > 0)
    }
  }

  stopifnot(cfg$upstream_screen$Se_vis >= 0, cfg$upstream_screen$Se_vis <= 1)
  stopifnot(cfg$upstream_screen$Sp_vis >= 0, cfg$upstream_screen$Sp_vis <= 1)

  if (!is.null(cfg$upstream_screen$f_asym_fixed)) {
    stopifnot(cfg$upstream_screen$f_asym_fixed >= 0, cfg$upstream_screen$f_asym_fixed <= 1)
  }

  stopifnot(cfg$border_testing$Se_mol >= 0, cfg$border_testing$Se_mol <= 1)
  stopifnot(cfg$border_testing$Se_subsample >= 0, cfg$border_testing$Se_subsample <= 1)
  stopifnot(cfg$tiers$N_small_threshold >= 0)
  cfg
}

#' Load YAML configuration.
#' @param path Path to YAML config file.
#' @return Validated config list.
load_cfg <- function(path = "config/default.yaml") {
  default_path <- "config/default.yaml"
  defaults <- yaml::read_yaml(default_path)
  cfg <- if (normalizePath(path, mustWork = FALSE) == normalizePath(default_path, mustWork = FALSE)) {
    defaults
  } else {
    user_cfg <- yaml::read_yaml(path)
    merge_lists(defaults, user_cfg)
  }
  validate_cfg(cfg)
}

#' Set global random seed.
#' @param seed Integer-like seed.
#' @return Invisibly returns seed.
set_seed <- function(seed) {
  set.seed(as.integer(seed))
  invisible(seed)
}
