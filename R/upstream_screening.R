#' Upstream visual screening and border conditioning

#' Partition infected indices into symptomatic and asymptomatic.
#' @param indices Infected unit indices.
#' @param f_asym Proportion asymptomatic.
#' @return List with \\code{symptomatic} and \\code{asymptomatic} integer vectors.
partition_infected <- function(indices, f_asym) {
  stopifnot(f_asym >= 0, f_asym <= 1)
  if (length(indices) == 0) {
    return(list(symptomatic = integer(0), asymptomatic = integer(0)))
  }
  asym_flag <- stats::rbinom(length(indices), size = 1, prob = f_asym) == 1
  list(
    symptomatic = indices[!asym_flag],
    asymptomatic = indices[asym_flag]
  )
}

#' Apply visual screening to symptomatic units.
#' @param symptomatic_idx Indices of symptomatic infected units.
#' @param Se_vis Visual sensitivity.
#' @return Indices of symptomatic infected units that leak through.
apply_visual_screen <- function(symptomatic_idx, Se_vis) {
  stopifnot(Se_vis >= 0, Se_vis <= 1)
  if (length(symptomatic_idx) == 0) {
    return(integer(0))
  }
  detected <- stats::rbinom(length(symptomatic_idx), size = 1, prob = Se_vis) == 1
  symptomatic_idx[!detected]
}

#' Condition shipment units to border stage after upstream screening.
#' @param N_units Total shipment units.
#' @param infected_idx Infected unit indices (within 1..N_units).
#' @param cfg Configuration list.
#' @return List with border population size and surviving infected indices in border indexing.
condition_to_border <- function(N_units, infected_idx, cfg) {
  stopifnot(N_units >= 0)
  if (N_units == 0) {
    return(list(N_asym = 0L, infected_border = integer(0), f_asym = NA_real_))
  }

  f_asym <- cfg$upstream_screen$f_asym_fixed
  if (is.null(f_asym)) {
    rng <- cfg$upstream_screen$f_asym_range
    f_asym <- stats::runif(1, min = rng[[1]], max = rng[[2]])
  }

  split <- partition_infected(infected_idx, f_asym = f_asym)
  leaked_symptomatic <- apply_visual_screen(split$symptomatic, Se_vis = cfg$upstream_screen$Se_vis)
  surviving_infected <- sort(c(split$asymptomatic, leaked_symptomatic))

  pass <- rep(TRUE, N_units)
  detected_symptomatic <- setdiff(split$symptomatic, leaked_symptomatic)
  if (length(detected_symptomatic) > 0) {
    pass[detected_symptomatic] <- FALSE
  }

  non_infected <- setdiff(seq_len(N_units), infected_idx)
  if (length(non_infected) > 0 && cfg$upstream_screen$Sp_vis < 1) {
    false_positive <- stats::rbinom(length(non_infected), size = 1, prob = 1 - cfg$upstream_screen$Sp_vis) == 1
    pass[non_infected[false_positive]] <- FALSE
  }

  border_indices <- which(pass)
  infected_border <- match(surviving_infected, border_indices)
  infected_border <- as.integer(infected_border[!is.na(infected_border)])

  list(
    N_asym = as.integer(length(border_indices)),
    infected_border = infected_border,
    f_asym = f_asym
  )
}
