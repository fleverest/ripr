#' Sampling family: the observation model
#'
#' A `family` owns the sampling model p_theta(x) and nothing else: no knowledge
#' of null hypotheses, alternatives, or optimisation. Concrete families provide
#' log-densities (vectorised over outcomes and over parameter columns), the
#' score, simulation, and -- for finite sample spaces -- full support
#' enumeration.
#' @export
family <- new_class("family", abstract = TRUE)

#' Log density log p_theta(x)
#'
#' Contract: `theta` is a parameter vector of length [param_dim()]. With
#' `x = NULL` (the hot path used by engines) the density is evaluated over the
#' family's enumerated support and returned as a length-M numeric vector. With
#' an explicit `x` (an `(N, d)` matrix of outcomes) a length-N vector is
#' returned.
#' @param family A `family`.
#' @param theta Parameter vector.
#' @param x Outcomes, or `NULL` for the enumerated support.
#' @return Numeric vector of log densities.
#' @export
log_density <- new_generic(
  "log_density",
  "family",
  function(family, theta, x = NULL) {
    S7::S7_dispatch()
  }
)

#' Batched log density over parameter columns
#'
#' Contract: `theta_mat` is `(d, C)` -- one parameter vector per column.
#' Returns an `(M, C)` (or `(N, C)` for explicit `x`) matrix with one
#' log-density column per parameter.
#' @param family A `family`.
#' @param theta_mat `(d, C)` matrix of parameter columns.
#' @param x Outcomes, or `NULL` for the enumerated support.
#' @return `(M, C)` matrix of log densities.
#' @keywords internal
log_density_batch <- new_generic(
  "log_density_batch",
  "family",
  function(family, theta_mat, x = NULL) {
    S7::S7_dispatch()
  }
)

#' Score function `d log p_theta(x) / d theta`
#'
#' Contract: returns an `(M, d)` matrix over the family's support (`x = NULL`)
#' or an `(N, d)` matrix for explicit outcomes. Used by the Frank-Wolfe and EM
#' face oracles for chain-ruled gradients.
#' @param family A `family`.
#' @param theta Parameter vector.
#' @param x Outcomes, or `NULL` for the enumerated support.
#' @return `(M, d)` matrix of per-outcome score contributions.
#' @export
score <- new_generic("score", "family", function(family, theta, x = NULL) {
  S7::S7_dispatch()
})

#' Dimension of the parameter vector
#' @param family A `family`.
#' @return Integer parameter dimension.
#' @export
param_dim <- new_generic("param_dim", "family", function(family) {
  S7::S7_dispatch()
})

#' Score used inside the EM M-step gradient
#'
#' Contract: the gradient of the responsibility-weighted log-likelihood
#' `sum_x w(x) log p_theta(x)` with respect to `theta`, as an `(M, d)` matrix of
#' per-outcome contributions. Defaults to [score()]; the multinomial family
#' adds back the `+n` constant because it does not vanish under weighted sums
#' whose total weight varies.
#' @param family A `family`.
#' @param theta Parameter vector.
#' @return `(M, d)` matrix of per-outcome M-step score contributions.
#' @keywords internal
em_score <- new_generic("em_score", "family", function(family, theta) {
  S7::S7_dispatch()
})

method(em_score, family) <- function(family, theta) {
  score(family, theta)
}

#' Draw observations from p_theta
#'
#' Contract: returns an `(n_obs, d)`-like numeric matrix of draws.
#' @param family A `family`.
#' @param theta Parameter vector.
#' @param n_obs Number of draws.
#' @return `(n_obs, d)` numeric matrix of draws.
#' @export
simulate <- new_generic(
  "simulate",
  "family",
  function(family, theta, n_obs) {
    S7::S7_dispatch()
  }
)

#' Enumerated support for finite families
#'
#' Contract: an `(M, d)` matrix listing every outcome with positive probability
#' under some parameter. Families with infinite sample spaces must signal an
#' informative error.
#' @param family A `family`.
#' @return `(M, d)` matrix of outcomes.
#' @export
support <- new_generic("support", "family", function(family) S7::S7_dispatch())

#' Is `p_theta(x)` a Bernstein basis function of `theta`?
#'
#' The capability behind the deterministic certified bound (see
#' [oracle_bound()]). Contract: return `NULL` when the family's density is not a
#' simplicial Bernstein basis function of its parameter, and otherwise
#' `list(n = , K = , tally = )` where `n` is the Bernstein degree, `K` the
#' number of barycentric coordinates (so `theta` ranges over the `(K-1)`-simplex
#' and `param_dim()` is `K`), and `tally` the `(M, K)` matrix of exponent
#' multi-indices **in the same row order as the family's [support()]** -- and
#' hence as `engine@outcomes` for an [exact_engine()] over that family.
#'
#' The capability belongs to the family, not the engine: being exact is
#' necessary but not sufficient. It is the *multinomial* pmf that happens to be
#' the degree-`n` Bernstein basis function at count vector `x`, so an exact
#' engine over some other finite-support family does not qualify. Families
#' without a method inherit the `NULL` default and degrade to the heuristic
#' oracle rather than being silently certified.
#' @param family A `family`.
#' @return `NULL`, or `list(n = , K = , tally = )`.
#' @keywords internal
bernstein_form <- new_generic("bernstein_form", "family", function(family) {
  S7::S7_dispatch()
})

method(bernstein_form, family) <- function(family) {
  NULL
}

#' Number of Bernstein coefficients implied by a [bernstein_form()]
#' @param bf A [bernstein_form()] result, or `NULL`.
#' @return Integer count, or `NA_integer_` for `NULL`.
#' @keywords internal
#' @noRd
n_coefficients <- function(bf) {
  if (is.null(bf)) {
    return(NA_integer_)
  }
  as.numeric(choose(bf$n + bf$K - 1L, bf$K - 1L))
}

#' Multinomial sampling family
#'
#' The n-trial, K-category multinomial. The support enumeration and log
#' multinomial coefficients are built once at construction; the support-path
#' methods serve every expectation as a base R weighted sum in log space.
#'
#' @param n_trials Integer. Total number of draws per observation.
#' @param k Integer. Number of categories.
#' @return A `multinomial_family`.
#' @export
multinomial_family <- new_class(
  "multinomial_family",
  parent = family,
  properties = list(
    n_trials = class_numeric,
    k = class_numeric,
    likelihood = class_list
  ),
  constructor = function(n_trials, k) {
    likelihood <- make_multinomial_likelihood(n_trials, k)
    new_object(
      S7_object(),
      n_trials = as.numeric(n_trials),
      k = as.numeric(k),
      likelihood = likelihood
    )
  }
)

method(log_density, multinomial_family) <- function(family, theta, x = NULL) {
  if (is.null(x)) {
    return(family@likelihood$log_pmf(theta))
  }
  as.vector(mnom_logpmf(
    as_outcome_matrix(x),
    matrix(log(theta), ncol = 1L),
    family@n_trials
  ))
}

method(log_density_batch, multinomial_family) <- function(
  family,
  theta_mat,
  x = NULL
) {
  if (is.null(x)) {
    return(family@likelihood$log_pmf_batch(theta_mat))
  }
  mnom_logpmf(as_outcome_matrix(x), log(theta_mat), family@n_trials)
}

method(score, multinomial_family) <- function(family, theta, x = NULL) {
  if (is.null(x)) {
    return(family@likelihood$score(theta))
  }
  sweep(as_outcome_matrix(x), 2L, theta, "/") - family@n_trials
}

method(param_dim, multinomial_family) <- function(family) {
  as.integer(family@k)
}

method(em_score, multinomial_family) <- function(family, theta) {
  family@likelihood$score(theta) + family@likelihood$n
}

method(simulate, multinomial_family) <- function(family, theta, n_obs) {
  t(rmultinom(n_obs, size = family@n_trials, prob = theta))
}

method(support, multinomial_family) <- function(family) {
  family@likelihood$support
}

method(bernstein_form, multinomial_family) <- function(family) {
  # P_theta(x) = choose(n; x) prod theta_j^{x_j} is exactly the degree-n
  # simplicial Bernstein basis function indexed by the count vector x.
  list(
    n = as.integer(family@n_trials),
    K = as.integer(family@k),
    tally = family@likelihood$support
  )
}

#' Normalise a sampling-model argument into a `family`
#'
#' Accepts a `family` (returned unchanged) and errors otherwise. Kept as the
#' single normalisation seam so problem constructors need not special-case
#' their `family` argument.
#' @param x A `family`.
#' @return A `family`.
#' @keywords internal
as_family <- function(x) {
  if (S7_inherits(x, family)) {
    return(x)
  }
  stop("cannot interpret this object as a sampling family: expected a `family`")
}
