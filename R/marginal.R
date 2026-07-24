#' Distribution over a family's sample space (abstract)
#'
#' The sample-space branch of [distribution]: a law over outcomes `X`. This is
#' the evaluable interface -- [dist_log_density()] and [dist_sample()] dispatch
#' here -- and the type the RIPr problem consumes as the alternative Q. A
#' [marginal] (the pushforward of a [mixing] through a [family]) is one
#' implementation; direct closed-form outcome laws can subclass this too.
#' @export
outcome_distribution <- new_class(
  "outcome_distribution",
  parent = distribution,
  abstract = TRUE
)

#' Outcome distribution induced by mixing a family
#'
#' A `marginal` is a self-contained [outcome_distribution]: the law of `X` when
#' `theta` is drawn from a [mixing] measure and `X ~ p_theta`, i.e.
#' `m(x) = integral p_theta(x) d(mixing)(theta)`. Unlike a bare `mixing`, it
#' carries its family, so [dist_log_density()] and [dist_sample()] need no
#' separate family argument. This is the object the alternative Q and the
#' projection P* become when evaluated on data (see [e_variable()]).
#'
#' Build one with [as_marginal()] rather than the raw constructor.
#' @param mixing A [mixing] measure over the parameter space.
#' @param family The sampling [family] whose kernel `p_theta` is mixed.
#' @return A `marginal`.
#' @export
marginal <- new_class(
  "marginal",
  parent = outcome_distribution,
  properties = list(
    mixing = mixing,
    family = family
  )
)

#' Pair a mixing measure with a family to get its induced outcome distribution
#'
#' @param mixing A [mixing] measure (or anything [as_mixing()] accepts).
#' @param family A sampling [family] (or anything [as_family()] accepts).
#' @return A [marginal] over the family's sample space.
#' @export
as_marginal <- function(mixing, family) {
  marginal(mixing = as_mixing(mixing), family = as_family(family))
}

#' Validate that the alternative is an outcome distribution
#'
#' The RIPr problem's Q lives over the *sample* space, so it must be an
#' [outcome_distribution] (typically a [marginal]). A bare [mixing] measure -- a
#' law over *parameters* -- is rejected with a hint to marginalise it first, so
#' the type at the boundary stays honest.
#' @param x An [outcome_distribution].
#' @return `x`, unchanged.
#' @keywords internal
as_outcome_distribution <- function(x) {
  if (S7_inherits(x, outcome_distribution)) {
    return(x)
  }
  if (S7_inherits(x, mixing)) {
    stop(
      "the alternative Q must be an `outcome_distribution` over the sample ",
      "space, not a bare `mixing` measure over parameters; marginalise it ",
      "first with `as_marginal(mixing, family)`"
    )
  }
  stop("cannot interpret this object as an outcome distribution")
}

#' Log density of a distribution over outcomes
#'
#' Contract: `x` is an `(N, d)` matrix of outcomes, or `NULL` for the family's
#' enumerated support. Returns a length-N numeric vector of `log p(x)`.
#' @param dist An [outcome_distribution].
#' @param x Outcomes, or `NULL` for the enumerated support.
#' @return Numeric vector of log densities.
#' @export
dist_log_density <- new_generic(
  "dist_log_density",
  "dist",
  function(dist, x = NULL) {
    S7::S7_dispatch()
  }
)

#' Draw outcomes from a distribution
#'
#' Contract: returns an `(n_obs, d)` numeric matrix of draws. This is the sample
#' source for Monte Carlo engines.
#' @param dist An [outcome_distribution].
#' @param n_obs Number of draws.
#' @return `(n_obs, d)` numeric matrix of draws.
#' @export
dist_sample <- new_generic(
  "dist_sample",
  "dist",
  function(dist, n_obs) {
    S7::S7_dispatch()
  }
)

method(dist_log_density, marginal) <- function(dist, x = NULL) {
  induced_log_density(dist@mixing, dist@family, x)
}

method(dist_sample, marginal) <- function(dist, n_obs) {
  induced_draw(dist@mixing, dist@family, n_obs)
}

# --- Induced-marginal formulas: double dispatch on (mixing, family) ----------
# These carry the actual integral `m(x) = int p_theta(x) d(mixing)(theta)` and
# its sampler, specialised by mixing kind (and, where a closed form exists, by
# family too -- see gaussian.R). They are the implementation layer behind the
# self-contained `marginal` interface above.

#' Induced log density `log int p_theta(x) d(mixing)(theta)`
#' @param mixing A [mixing] measure.
#' @param family A [family].
#' @param x Outcomes, or `NULL` for the enumerated support.
#' @return Numeric vector of log densities.
#' @keywords internal
induced_log_density <- new_generic(
  "induced_log_density",
  c("mixing", "family"),
  function(mixing, family, x = NULL) {
    S7::S7_dispatch()
  }
)

#' Draw outcomes from the induced marginal
#' @param mixing A [mixing] measure.
#' @param family A [family].
#' @param n_obs Number of draws.
#' @return `(n_obs, d)` numeric matrix of draws.
#' @keywords internal
induced_draw <- new_generic(
  "induced_draw",
  c("mixing", "family"),
  function(mixing, family, n_obs) {
    S7::S7_dispatch()
  }
)

method(induced_log_density, list(point_mixing, family)) <- function(
  mixing,
  family,
  x = NULL
) {
  log_density(family, mixing@theta_star, x)
}

method(induced_draw, list(point_mixing, family)) <- function(
  mixing,
  family,
  n_obs
) {
  simulate(family, mixing@theta_star, n_obs)
}

method(induced_log_density, list(finite_mixing, family)) <- function(
  mixing,
  family,
  x = NULL
) {
  row_logsumexp(sweep(
    log_density_batch(family, mixing@components, x),
    2L,
    log(mixing@weights),
    "+"
  ))
}

method(induced_draw, list(finite_mixing, family)) <- function(
  mixing,
  family,
  n_obs
) {
  idx <- sample.int(
    length(mixing@weights),
    n_obs,
    replace = TRUE,
    prob = mixing@weights
  )
  do.call(
    rbind,
    lapply(idx, function(c_i) {
      simulate(family, mixing@components[, c_i], 1L)
    })
  )
}
