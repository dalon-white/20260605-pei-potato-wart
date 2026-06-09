#' Sampling rule calculations

#' Binomial zero-acceptance sample size.
#' @param p0 Design prevalence threshold.
#' @param alpha False-negative tail target.
#' @return Integer sample size.
binomial_n <- function(p0, alpha) {
  stopifnot(p0 > 0, p0 < 1, alpha > 0, alpha < 1)
  as.integer(ceiling(log(alpha) / log(1 - p0)))
}

#' Hypergeometric zero-detection probability.
#' @param N Population size.
#' @param D Number infected.
#' @param n Sample size.
#' @return Probability of zero infected in sample.
hypergeom_zero_prob <- function(N, D, n) {
  # If there are no infected units (or a nonpositive count), then every sample has zero infected. So the probability is exactly 1
  if (D <= 0) {
    return(1)
  }
  # n cannot be larger than N, and cannot be larger than the number of non-infected units (N - D), otherwise zero infected in sample is impossible, so probability is exactly 0
  if (n > N) {
    return(0)
  }
  if (n > (N - D)) {
    return(0)
  }
  exp(lchoose(N - D, n) - lchoose(N, n))
# Quick intuition example
# If N = 100, D = 10, n = 5:

# favorable samples = choose 5 from 90 uninfected
# total samples = choose 5 from 100
# probability is that ratio, about 0.584
# So there is about a 58.4% chance you miss infection entirely in that sample size.
}

#' Exact finite-pop sample size for zero acceptance.
#' @param N Population size.
#' @param p0 Design prevalence threshold.
#' @param alpha False-negative tail target.
#' @return Minimum integer n in [1, N].
hypergeom_n <- function(N, p0, alpha) {
# Hypergeometric sampling generation function - determines how many samples to take from a finite population of size N to have a probability of at most alpha of missing infection when the true prevalence is p0. It finds the smallest n such that the probability of zero infected in the sample is less than or equal to alpha, given D = ceiling(p0 * N) infected units in the population.
#  Used for building test_sampling.R
  stopifnot(N >= 1, p0 > 0, p0 < 1, alpha > 0, alpha < 1)
  N <- as.integer(N)
  D <- max(1L, as.integer(ceiling(p0 * N)))

  # This computes the probability of missing all infected units for every possible sample size from 1 to N.
  probs <- vapply(seq_len(N), function(n) hypergeom_zero_prob(N, D, n), numeric(1))

  # This is an exact finite-population sample-size calculator for a zero-acceptance rule.
  #   “Zero-acceptance” means if any infected unit is detected, the lot fails.
  #   “Finite-population” means the shipment is treated as a fixed set of units, sampled without replacement.
  #   “Exact” means it uses the hypergeometric model rather than an approximation.
  idx <- which(probs <= alpha)
  if (length(idx) == 0) {
    return(N)
  }
  as.integer(idx[[1]])
}

#' Effective per-unit sensitivity.
#' @param cfg Configuration list.
#' @return Effective sensitivity in [0, 1].
Se_eff <- function(cfg) {
# This calculates the effective sensitivity per sampled unit. This is adjusted to account for testing errors
#     Se_mol: sensitivity of the molecular test itself
#     Se_subsample: sensitivity loss from the subsampling process
  cfg$border_testing$Se_mol * cfg$border_testing$Se_subsample
}

#' Inflate sample size for imperfect sensitivity under binomial approximation.
#' @param p0 Design prevalence threshold.
#' @param alpha False-negative tail target.
#' @param Se_eff Effective sensitivity.
#' @return Integer sample size, or Inf when no sensitivity.
inflate_n_for_sensitivity <- function(p0, alpha, Se_eff) {
  # Without imperfect sensitivity, the miss probability is approximately:
  #   (1 - p0)^n
  #
  # With imperfect sensitivity, each sampled unit has probability p0 * Se_eff of being both:
  #   - infected
  #   - successfully detected
  stopifnot(p0 > 0, p0 < 1, alpha > 0, alpha < 1, Se_eff >= 0, Se_eff <= 1)
  if (Se_eff == 0) {
    return(Inf)
  }
  as.integer(ceiling(log(alpha) / log(1 - p0 * Se_eff)))
  # This is a binomial approximation. It does not account for finite-population sampling without replacement. It is usually reasonable, especially for larger populations, but it is not exact in the same way hypergeom_n() is
}


#' Compute required border sample size under tiered policy.
#' @param N_asym Asymptomatic population size at border.
#' @param p0 Design prevalence threshold.
#' @param alpha False-negative tail target.
#' @param Se_eff Effective sensitivity.
#' @param tiers_cfg Tier configuration.
#' @return Integer required n.
border_required_n <- function(N_asym, p0, alpha, Se_eff, tiers_cfg) {
  # If there are no asymptomatic units, then no sampling is needed, so return 0. This also guards against negative values which would not make sense in this context.
  stopifnot(N_asym >= 0)
  if (N_asym == 0) {
    return(0L)
  }

  # If the shipment is small enough, sample all asymptomatic units (and skip calculations).
  if (N_asym <= tiers_cfg$N_small_threshold) {
    return(as.integer(N_asym))
  }

  n_target <- 0L
  # For medium-sized lots, compute two candidate sample sizes:
  #   n_h: exact finite-population sample size from hypergeometric logic
  #   n_s: sensitivity-adjusted sample size from the binomial approximation
  # Then take the larger of the two.
  if (isTRUE(tiers_cfg$use_hypergeom_for_medium)) {
    n_h <- hypergeom_n(N = N_asym, p0 = p0, alpha = alpha)
    n_s <- inflate_n_for_sensitivity(p0 = p0, alpha = alpha, Se_eff = Se_eff)
    n_target <- max(n_h, n_s)
  } else if (isTRUE(tiers_cfg$use_binomial_floor_for_large)) {
    # For large lots, the policy may switch to the sensitivity-adjusted binomial approximation only.
    # Why: for large populations, hypergeometric and binomial results are often close, and the approximation is simpler.
    n_target <- inflate_n_for_sensitivity(p0 = p0, alpha = alpha, Se_eff = Se_eff)
  } else {
    # Fallback case: use the simplest binomial sample-size rule with perfect detection assumption.
    n_target <- binomial_n(p0 = p0, alpha = alpha)
  }

  # Final safeguard: required sample size cannot exceed the number of available asymptomatic units.
  as.integer(min(N_asym, n_target))
    #So if the formula says sample 500, but there are only 320 asymptomatic units, the rule returns 320, b/c you cannot sample more than exist.
}



# Explanation:

# hypergeom_n() asks: how many samples are needed if sampling is random and detection is perfect?
# Se_eff() asks: how good is the detection process once a unit is sampled?
# inflate_n_for_sensitivity() asks: how much more sampling is needed because detection is imperfect?
# border_required_n() combines those ideas into the actual tiered border-sampling rule.