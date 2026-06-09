#' Metrics and operating characteristics

#' Build required sample table by border population size.
#' @param shipment_df Shipment-level simulation data.
#' @return Tibble mapping N_asym to required n summary.
required_n_table <- function(shipment_df) {
  shipment_df |>
    dplyr::group_by(scenario, N_asym) |>
    dplyr::summarise(
      required_n = as.integer(stats::median(n_required)),
      .groups = "drop"
    ) |>
    dplyr::arrange(scenario, N_asym)
}

#' Summarize annual-level metrics.
#' @param shipment_df Shipment-level simulation data.
#' @param shipments_possible_basis Optional baseline shipment count for scaling.
#' @return Tibble with annual summary fields.
annual_metrics <- function(shipment_df, shipments_possible_basis = NA_integer_) {
  shipment_df |>
    dplyr::group_by(scenario) |>
    dplyr::summarise(
      shipments = dplyr::n(),
      shipments_filtered_upstream = sum(shipment_filtered_upstream),
      shipments_reaching_border = sum(shipment_reaches_border),
      upstream_filter_rate = mean(shipment_filtered_upstream),
      shipments_possible_basis = as.integer(shipments_possible_basis),
      projected_filtered_upstream = ifelse(
        is.na(shipments_possible_basis),
        NA_real_,
        upstream_filter_rate * shipments_possible_basis
      ),
      infected_shipments = sum(infected_present),
      detected_infected_shipments = sum(detected & infected_present),
      realized_any_undetected = as.numeric(any(pass_undetected)),
      false_negative_rate_given_infected = ifelse(sum(infected_present) == 0, 0, mean(pass_undetected[infected_present])),
      prob_any_infected_pass_undetected = ifelse(
        infected_shipments == 0,
        0,
        1 - (1 - false_negative_rate_given_infected)^infected_shipments
      ),
      expected_total_tests = sum(samples_taken),
      n_median = stats::median(n_required),
      n_q25 = stats::quantile(n_required, 0.25),
      n_q75 = stats::quantile(n_required, 0.75),
      n_max = max(n_required),
      detection_prob_given_infected = ifelse(sum(infected_present) == 0, NA_real_, mean(detected[infected_present])),
      false_negative_prob_per_shipment = ifelse(sum(infected_present) == 0, NA_real_, mean(pass_undetected[infected_present])),
      .groups = "drop"
    )
}

#' Summarize upstream shipment filtering outcomes.
#' @param shipment_df Shipment-level simulation data.
#' @param shipments_possible_basis Optional baseline shipment count for scaling.
#' @return Tibble with shipment filtering outcomes by scenario.
upstream_filtering_table <- function(shipment_df, shipments_possible_basis = NA_integer_) {
  shipment_df |>
    dplyr::group_by(scenario) |>
    dplyr::summarise(
      shipments_simulated = dplyr::n(),
      shipments_filtered_upstream = sum(shipment_filtered_upstream),
      shipments_reaching_border = sum(shipment_reaches_border),
      shipments_with_true_upstream_detection = sum(upstream_true_detected_count > 0),
      shipments_with_false_positive_visuals = sum(upstream_false_positive_count > 0),
      upstream_filter_rate = mean(shipment_filtered_upstream),
      shipments_possible_basis = as.integer(shipments_possible_basis),
      projected_filtered_upstream = ifelse(
        is.na(shipments_possible_basis),
        NA_real_,
        upstream_filter_rate * shipments_possible_basis
      ),
      projected_reaching_border = ifelse(
        is.na(shipments_possible_basis),
        NA_real_,
        shipments_possible_basis - projected_filtered_upstream
      ),
      .groups = "drop"
    )
}

#' Generate operating characteristic table.
#' @param p_values Prevalence values.
#' @param n_values Sample size values.
#' @param Se_eff Effective sensitivity.
#' @return Tibble with confidence by n and p.
oc_curve_table <- function(p_values, n_values, Se_eff) {
  tidyr::crossing(p = p_values, n = n_values) |>
    dplyr::mutate(confidence = 1 - (1 - p * Se_eff)^n)
}
