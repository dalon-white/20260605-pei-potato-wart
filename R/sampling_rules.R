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
  if (D <= 0) {
    return(1)
  }
  if (n > N) {
    return(0)
  }
  if (n > (N - D)) {
    return(0)
  }
  exp(lchoose(N - D, n) - lchoose(N, n))
}

#' Exact finite-pop sample size for zero acceptance.
#' @param N Population size.
#' @param p0 Design prevalence threshold.
#' @param alpha False-negative tail target.
#' @return Minimum integer n in [1, N].
hypergeom_n <- function(N, p0, alpha) {
  stopifnot(N >= 1, p0 > 0, p0 < 1, alpha > 0, alpha < 1)
  N <- as.integer(N)
  D <- max(1L, as.integer(ceiling(p0 * N)))
  probs <- vapply(seq_len(N), function(n) hypergeom_zero_prob(N, D, n), numeric(1))
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
  cfg$border_testing$Se_mol * cfg$border_testing$Se_subsample
}

#' Inflate sample size for imperfect sensitivity under binomial approximation.
#' @param p0 Design prevalence threshold.
#' @param alpha False-negative tail target.
#' @param Se_eff Effective sensitivity.
#' @return Integer sample size, or Inf when no sensitivity.
inflate_n_for_sensitivity <- function(p0, alpha, Se_eff) {
  stopifnot(p0 > 0, p0 < 1, alpha > 0, alpha < 1, Se_eff >= 0, Se_eff <= 1)
  if (Se_eff == 0) {
    return(Inf)
  }
  as.integer(ceiling(log(alpha) / log(1 - p0 * Se_eff)))
}

#' Compute required border sample size under tiered policy.
#' @param N_asym Asymptomatic population size at border.
#' @param p0 Design prevalence threshold.
#' @param alpha False-negative tail target.
#' @param Se_eff Effective sensitivity.
#' @param tiers_cfg Tier configuration.
#' @return Integer required n.
border_required_n <- function(N_asym, p0, alpha, Se_eff, tiers_cfg) {
  stopifnot(N_asym >= 0)
  if (N_asym == 0) {
    return(0L)
  }

  if (N_asym <= tiers_cfg$N_small_threshold) {
    return(as.integer(N_asym))
  }

  n_target <- 0L
  if (isTRUE(tiers_cfg$use_hypergeom_for_medium)) {
    n_h <- hypergeom_n(N = N_asym, p0 = p0, alpha = alpha)
    n_s <- inflate_n_for_sensitivity(p0 = p0, alpha = alpha, Se_eff = Se_eff)
    n_target <- max(n_h, n_s)
  } else if (isTRUE(tiers_cfg$use_binomial_floor_for_large)) {
    n_target <- inflate_n_for_sensitivity(p0 = p0, alpha = alpha, Se_eff = Se_eff)
  } else {
    n_target <- binomial_n(p0 = p0, alpha = alpha)
  }

  as.integer(min(N_asym, n_target))
}
