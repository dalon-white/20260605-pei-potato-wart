#' Distribution helpers for shipments and units

#' Draw from a clamped lognormal distribution.
#' @param n Number of draws.
#' @param mean Target arithmetic mean of the unclamped distribution.
#' @param sd Log-scale SD (sdlog).
#' @param min Minimum bound.
#' @param max Maximum bound.
#' @return Numeric vector.
r_trunc_lognorm <- function(n, mean, sd, min, max) {
  stopifnot(n >= 1, mean > 0, sd > 0, min > 0, max >= min)
  meanlog <- log(mean) - 0.5 * sd^2
  draws <- stats::rlnorm(n = n, meanlog = meanlog, sdlog = sd)
  pmin(pmax(draws, min), max)
}

#' Convert ounces to grams.
#' @param x Ounces.
#' @return Grams.
to_grams_oz <- function(x) {
  x * 28.3495
}

#' Derive units per shipment from mass and unit weight.
#' @param shipment_kg Shipment mass in kilograms.
#' @param unit_oz Unit mass in ounces.
#' @return Integer unit counts.
units_per_shipment <- function(shipment_kg, unit_oz) {
  units <- floor((shipment_kg * 1000) / to_grams_oz(unit_oz))
  as.integer(pmax(1, units))
}

#' Derive units per shipment from shipment mass and mean unit weight.
#' @param shipment_kg Shipment mass in kilograms.
#' @param unit_cfg Unit-size configuration with mean, sd, min, and max.
#' @return Integer unit counts.
shipment_unit_count <- function(shipment_kg, unit_cfg) {
  stopifnot(is.list(unit_cfg), !is.null(unit_cfg$mean), unit_cfg$mean > 0)
  units_per_shipment(shipment_kg = shipment_kg, unit_oz = unit_cfg$mean)
}

#' Draw overdispersed shipment counts using negative binomial.
#' @param mean Mean count.
#' @param dispersion Dispersion parameter where size = 1/dispersion.
#' @return Integer shipment count.
r_overdispersed_counts <- function(mean, dispersion) {
  stopifnot(mean > 0, dispersion > 0)
  size <- 1 / dispersion
  prob <- size / (size + mean)
  as.integer(stats::rnbinom(n = 1, size = size, prob = prob))
}

#' Convert mean and concentration to Beta shape parameters.
#' @param mu Mean prevalence in (0, 1).
#' @param kappa Concentration (>0).
#' @return List with alpha and beta.
alpha_beta_from_mean_conc <- function(mu, kappa) {
  stopifnot(mu > 0, mu < 1, kappa > 0)
  list(alpha = mu * kappa, beta = (1 - mu) * kappa)
}
