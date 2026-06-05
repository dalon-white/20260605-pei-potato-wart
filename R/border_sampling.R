#' Border sampling operations

#' Sample indices without replacement.
#' @param N_asym Border-stage asymptomatic population size.
#' @param n Number to sample.
#' @return Integer sampled indices.
sample_without_replacement <- function(N_asym, n) {
  stopifnot(N_asym >= 0, n >= 0)
  if (N_asym == 0 || n == 0) {
    return(integer(0))
  }
  sample.int(N_asym, size = min(n, N_asym), replace = FALSE)
}

#' Apply border sampling and molecular detection.
#' @param N_asym Border-stage asymptomatic population size.
#' @param n Number to sample.
#' @param infected_indices Infected positions within the border population.
#' @param Se_eff Effective per-unit sensitivity.
#' @return List with decision and detection details.
apply_border_sampling <- function(N_asym, n, infected_indices, Se_eff) {
  sampled <- sample_without_replacement(N_asym = N_asym, n = n)
  infected_sampled <- sampled[sampled %in% infected_indices]
  positives <- 0L

  if (length(infected_sampled) > 0) {
    positives <- sum(stats::rbinom(length(infected_sampled), size = 1, prob = Se_eff))
  }

  detected <- positives > 0
  list(
    sampled_indices = sampled,
    infected_sampled = infected_sampled,
    positives = as.integer(positives),
    detected = detected,
    accept = !detected,
    samples_taken = as.integer(length(sampled))
  )
}
