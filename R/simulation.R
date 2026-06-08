#' Year-level simulation workflow

#' Simulate one shipment.
#' @param cfg Configuration list.
#' @param scenario Scenario list with mu, kappa, and label.
#' @return One-row tibble.
simulate_shipment <- function(cfg, scenario) {
  if (is.atomic(scenario)) {
    scenario <- as.list(scenario)
  }
  if (!is.list(scenario)) {
    stop("`scenario` must be a list or named atomic vector.")
  }

  mu <- as.numeric(scenario$mu %||% cfg$prevalence$mu)
  kappa <- as.numeric(scenario$kappa %||% cfg$prevalence$kappa)
  label <- as.character(scenario$label %||% sprintf("mu=%.4f_kappa=%s", mu, kappa))

  mass <- r_trunc_lognorm(
    n = 1,
    mean = cfg$shipment_mass_kg$mean,
    sd = cfg$shipment_mass_kg$sd,
    min = cfg$shipment_mass_kg$min,
    max = cfg$shipment_mass_kg$max
  )

  mean_unit_size_oz <- as.numeric(cfg$unit_size_oz$mean)
  N_units <- shipment_unit_count(shipment_kg = mass, unit_cfg = cfg$unit_size_oz)
  ab <- alpha_beta_from_mean_conc(mu = mu, kappa = kappa)
  p_s <- stats::rbeta(1, shape1 = ab$alpha, shape2 = ab$beta)
  infected_n <- stats::rbinom(1, size = N_units, prob = p_s)
  infected_idx <- if (infected_n > 0) sample.int(N_units, infected_n) else integer(0)

  up <- condition_to_border(N_units = N_units, infected_idx = infected_idx, cfg = cfg)
  se <- Se_eff(cfg)

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
    infected_count_border = as.integer(length(up$infected_border)),
    f_asym = up$f_asym,
    n_required = as.integer(n_required),
    samples_taken = as.integer(border$samples_taken),
    detected = as.logical(border$detected),
    infected_present = as.logical(infected_present),
    pass_undetected = as.logical(pass_undetected)
  )
}

#' Simulate one year of shipments.
#' @param seed Random seed.
#' @param cfg Configuration list.
#' @param scenario Optional list with mu and kappa overrides.
#' @return List with shipment and annual summary data frames.
simulate_year <- function(seed, cfg, scenario = NULL) {
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

  n_shipments <- max(1L, r_overdispersed_counts(
    mean = cfg$shipments_per_year$mean,
    dispersion = cfg$shipments_per_year$dispersion
  ))

  shipment_df <- purrr::map_dfr(seq_len(n_shipments), ~simulate_shipment(cfg = cfg, scenario = sc))
  list(
    shipment_level = shipment_df,
    annual_summary = annual_metrics(shipment_df)
  )
}
