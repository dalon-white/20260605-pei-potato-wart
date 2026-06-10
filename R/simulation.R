#' Year-level simulation workflow

#' Simulate one shipment.
#' @param cfg Configuration list.
#' @param scenario Scenario list with mu, kappa, and label.
#' @return One-row tibble.
simulate_shipment <- function(cfg, scenario) {
  # scenario is from the cfg; various scenarios are selected. They contain mu (prevalence mean) and kappa (overdispersion) that define the underlying infection distribution in the shipment, as well as a label for reporting. The scenario can be a list or a named atomic vector (e.g., from yaml), and it can override the default mu and kappa from the cfg. If scenario is NULL, then defaults from cfg are used.
  if (is.atomic(scenario)) {
    scenario <- as.list(scenario)
  }
  if (!is.list(scenario)) {
    stop("`scenario` must be a list or named atomic vector.")
  }

  mu <- as.numeric(scenario$mu %||% cfg$prevalence$mu) # if the value is not supplied for the scenario, it will fallback to the default mu
  kappa <- as.numeric(scenario$kappa %||% cfg$prevalence$kappa) # if the value is not supplied for the scenario, it will fallback to the default kappa
  label <- as.character(scenario$label %||% sprintf("mu=%.4f_kappa=%s", mu, kappa)) # if it is not named, it will create one

  # randomly sample for the shipment size using the appropriate distribution and parameters from config
  mass <- r_trunc_lognorm(
    n = 1,
    mean = cfg$shipment_mass_kg$mean,
    sd = cfg$shipment_mass_kg$sd,
    min = cfg$shipment_mass_kg$min,
    max = cfg$shipment_mass_kg$max
  )

  # gather mean unit size from config, then calculate the numebr of units using the total mass and mean unit size
  mean_unit_size_oz <- as.numeric(cfg$unit_size_oz$mean)
  N_units <- shipment_unit_count(shipment_kg = mass, unit_cfg = cfg$unit_size_oz)

  # calculate the distribution shape of the Beta distribution for shipment-level prevalence
    # large Kappa: beta is tight around mu (shipment prevalence is likely close to mu); small kappa: beta is wide around mu (shipment prevalence is more variable, with more chance of being much higher or lower than mu)
      # mu = 0.01, kappa = 1000 -> beta is very tight around 0.01, so most shipments will have prevalence close to 1%
      # mu = 0.01, kappa = 1 -> beta is very wide, so shipments could have prevalence anywhere from near 0% to much higher than 1%
  ab <- alpha_beta_from_mean_conc(mu = mu, kappa = kappa)

  # randomly sample beta distributino to get the actual prevalence for this shipment
  p_s <- stats::rbeta(1, shape1 = ab$alpha, shape2 = ab$beta)
  # sample from binomial to get the number of infected units
  infected_n <- stats::rbinom(1, size = N_units, prob = p_s)
  # then randomly assign which units are infected
  infected_idx <- if (infected_n > 0) sample.int(N_units, infected_n) else integer(0)

  up <- condition_to_border(N_units = N_units, infected_idx = infected_idx, cfg = cfg)
  se <- Se_eff(cfg)

  border <- list(samples_taken = 0L, detected = FALSE)
  n_required <- 0L
  infected_present <- FALSE
  pass_undetected <- FALSE

  if (isTRUE(up$shipment_reaches_border)) {
    n_required <- border_required_n(
      N_asym = up$N_asym,
      p0 = cfg$p0,
      alpha = cfg$alpha,
      Se_eff = se,
      tiers_cfg = cfg$tiers
    )

    border <- apply_border_sampling(
      N_asym = up$N_asym,
      n = n_required,
      infected_indices = up$infected_border,
      Se_eff = se
    )

    infected_present <- length(up$infected_border) > 0
    pass_undetected <- infected_present && !border$detected
  }

  tibble::tibble(
    scenario = label,
    mu = mu,
    kappa = kappa,
    shipment_mass_kg = mass,
    mean_unit_size_oz = mean_unit_size_oz,
    N_units = as.integer(N_units),
    N_asym = as.integer(up$N_asym),
    p_shipment = p_s,
    infected_count_pre_upstream = as.integer(infected_n),
    shipment_filtered_upstream = as.logical(up$shipment_filtered_upstream),
    shipment_reaches_border = as.logical(up$shipment_reaches_border),
    upstream_true_detected_count = as.integer(up$upstream_true_detected_count),
    upstream_false_positive_count = as.integer(up$upstream_false_positive_count),
    upstream_detected = as.logical(up$upstream_true_detected_count > 0),
    infected_count_border = as.integer(length(up$infected_border)),
    f_asym = up$f_asym,
    n_required = as.integer(n_required),
    samples_taken = as.integer(border$samples_taken),
    border_detected = as.logical(border$detected),
    detected = as.logical(border$detected),
    infected_present = as.logical(infected_present),
    pass_undetected = as.logical(pass_undetected)
  )
}

#' Simulate one or more years of shipments.
#' @param seed Random seed.
#' @param cfg Configuration list.
#' @param scenario Optional list with mu and kappa overrides.
#' @param n_years Number of years to simulate.
#' @return List with per-year results and across-year summary statistics.
simulate_year <- function(seed, cfg, scenario = NULL, n_years = NULL) {
  set_seed(seed)

  normalize_scenario <- function(scn, cfg) {
    if (is.null(scn)) {
      return(list(
        mu = cfg$prevalence$mu,
        kappa = cfg$prevalence$kappa,
        label = sprintf("mu=%.4f_kappa=%s", cfg$prevalence$mu, cfg$prevalence$kappa)
      ))
    }

    if (is.atomic(scn)) {
      scn <- as.list(scn)
    }

    if (!is.list(scn)) {
      stop("`scenario` must be NULL, a list, or a named atomic vector.")
    }

    mu <- scn$mu %||% cfg$prevalence$mu
    kappa <- scn$kappa %||% cfg$prevalence$kappa
    label <- scn$label %||% sprintf("mu=%.4f_kappa=%s", mu, kappa)

    list(mu = as.numeric(mu), kappa = as.numeric(kappa), label = as.character(label))
  }

  sc <- normalize_scenario(scenario, cfg)

  years_raw <- n_years %||% cfg$simulation_years %||% 1L
  if (!is.numeric(years_raw) || length(years_raw) != 1 || !is.finite(years_raw) || years_raw < 1 || years_raw != as.integer(years_raw)) {
    stop("`n_years` must be a positive integer.")
  }
  years_to_simulate <- as.integer(years_raw)

  simulate_one_year <- function(year_index) {
    n_shipments <- max(1L, r_overdispersed_counts(
      mean = cfg$shipments_per_year$mean,
      dispersion = cfg$shipments_per_year$dispersion
    ))

    shipment_df <- purrr::map_dfr(seq_len(n_shipments), ~simulate_shipment(cfg = cfg, scenario = sc)) |>
      dplyr::mutate(year = as.integer(year_index))

    annual_df <- annual_metrics(
      shipment_df,
      shipments_possible_basis = cfg$shipments_per_year$mean
    ) |>
      dplyr::mutate(year = as.integer(year_index), .before = 1)

    list(shipment_level = shipment_df, annual_summary = annual_df)
  }

  year_results <- purrr::map(seq_len(years_to_simulate), simulate_one_year)
  shipment_df <- purrr::map_dfr(year_results, "shipment_level")
  annual_summary <- purrr::map_dfr(year_results, "annual_summary")
  annual_by_year <- annual_by_year_summary(annual_summary)

  list(
    years = year_results,
    shipment_level = shipment_df,
    annual_by_year = annual_by_year,
    annual_summary = annual_summary
  )
}
