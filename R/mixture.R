#' @include mixing_measure.R family.R
NULL

#' Distribution over a family's sample space (abstract)
#'
#' The sample-space branch of [distribution]: a law over outcomes, and the type
#' the RIPr problem consumes as the alternative \eqn{Q}{Q}. A [mixture] is one
#' implementation; a closed-form outcome law can subclass this directly with no
#' mixing measure behind it.
#' @examples
#' # `outcome_distribution` is abstract; [mixture()] subclasses it, e.g.
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' Q <- mixture(point_mixing(c(0.5, 0.3, 0.2)), fam)
#' S7::S7_inherits(Q, outcome_distribution)
#' @export
outcome_distribution <- new_class(
  "outcome_distribution",
  parent = distribution,
  abstract = TRUE
)


#' Log density of an outcome distribution
#' @param dist An [outcome_distribution].
#' @param x `(M, K)` matrix of outcomes, or a length-`K` vector for one outcome.
#' @return Length-`M` numeric vector.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' Q <- mixture(point_mixing(c(0.5, 0.3, 0.2)), fam)
#' dist_log_density(Q, c(2L, 1L, 1L))
#' @export
dist_log_density <- new_generic("dist_log_density", "dist", function(dist, x) {
  S7::S7_dispatch()
})


#' Draw outcomes from a distribution
#' @param dist An [outcome_distribution].
#' @param n_obs Number of draws.
#' @return `(n_obs, K)` numeric matrix.
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' Q <- mixture(point_mixing(c(0.5, 0.3, 0.2)), fam)
#' dist_draw(Q, n_obs = 3L)
#' @export
dist_draw <- new_generic("dist_draw", "dist", function(dist, n_obs) {
  S7::S7_dispatch()
})


#' The family a distribution lives over, if it has one
#'
#' Returns `NULL` when there is none. \eqn{Q}{Q} is fixed a priori and need not
#' be a mixture over the parameter space, so the family cannot in general be
#' recovered from it.
#' @param dist An [outcome_distribution].
#' @return A [sampling_family], or `NULL`.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' Q <- mixture(point_mixing(c(0.5, 0.3, 0.2)), fam)
#' dist_family(Q)
#' @export
dist_family <- new_generic("dist_family", "dist", function(dist) {
  S7::S7_dispatch()
})


method(dist_family, outcome_distribution) <- function(dist) NULL


#' Mixture `P_W` induced by a mixing measure and a family
#'
#' The law of `X` when \eqn{\theta \sim W}{theta ~ W} and
#' \eqn{X \sim p_\theta}{X ~ p_theta}. Unlike a bare [mixing_measure] it carries
#' its family, so it is evaluable on data. Both \eqn{Q = P_{W_1}}{Q = P_W1} and
#' \eqn{\widehat{P}^* = P_{\widehat{W}_0}}{P_star_hat} are mixtures.
#'
#' @param mixing A [mixing_measure] over the parameter space.
#' @param family The [sampling_family] whose kernel is mixed.
#' @return A `mixture`.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' mixture(point_mixing(c(0.5, 0.3, 0.2)), fam)
#' @export
mixture <- new_class(
  "mixture",
  parent = outcome_distribution,
  properties = list(mixing = mixing_measure, family = sampling_family)
)


method(dist_family, mixture) <- function(dist) dist@family


method(dist_log_density, mixture) <- function(dist, x) {
  induced_log_density(dist@mixing, dist@family, x)
}


method(dist_draw, mixture) <- function(dist, n_obs) {
  induced_draw(dist@mixing, dist@family, n_obs)
}


# --- Induced-mixture formulas: double dispatch on (mixing, family) -----------

#' Induced log density `log int p_theta(x) dW(theta)`
#' @param mixing A [mixing_measure].
#' @param family A [sampling_family].
#' @param x `(M, K)` matrix of outcomes.
#' @return Length-`M` numeric vector.
#' @keywords internal
induced_log_density <- new_generic(
  "induced_log_density",
  c("mixing", "family"),
  function(mixing, family, x) S7::S7_dispatch()
)


#' Draw from the induced mixture
#' @param mixing A [mixing_measure].
#' @param family A [sampling_family].
#' @param n_obs Number of draws.
#' @return `(n_obs, K)` numeric matrix.
#' @keywords internal
induced_draw <- new_generic(
  "induced_draw",
  c("mixing", "family"),
  function(mixing, family, n_obs) S7::S7_dispatch()
)


method(induced_log_density, list(point_mixing, sampling_family)) <- function(
  mixing,
  family,
  x
) {
  log_density(family, mixing@theta_star, x)
}


method(induced_draw, list(point_mixing, sampling_family)) <- function(
  mixing,
  family,
  n_obs
) {
  draw(family, mixing@theta_star, n_obs)
}


method(induced_log_density, list(finite_mixing, sampling_family)) <- function(
  mixing,
  family,
  x
) {
  row_logsumexp(add_by_col(
    log_density_batch(family, mixing@components, x),
    log(mixing@weights)
  ))
}


method(induced_draw, list(finite_mixing, sampling_family)) <- function(
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
  # One call per distinct component, indexed back into place: samplers are
  # vectorised over draws but not over parameters.
  k <- ncol(draw(family, mixing@components[, idx[1L]], 1L))
  out <- matrix(NA_real_, nrow = n_obs, ncol = k)
  for (c_i in unique(idx)) {
    rows <- which(idx == c_i)
    out[rows, ] <- draw(family, mixing@components[, c_i], length(rows))
  }
  out
}
