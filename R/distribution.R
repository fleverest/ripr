#' @include space.R family.R
NULL

#' Distributions over a space
#'
#' A `distribution` is a law over a [space]: the type the RIPr problem
#' consumes as the alternative \eqn{Q}{Q}, and the type it returns as
#' \eqn{\widehat{P}^*}{P_star_hat}. It answers two questions, [log_density()]
#' and [draw()], over the points its `sample_space` admits.
#'
#' \eqn{Q}{Q} is fixed a priori and need not have come from a family at all: a
#' and is usable everywhere an [induced_distribution()] is.
#' @param sample_space The [space] this is a law over.
#' @examples
#' # `distribution` is abstract; induced_distribution() subclasses it, e.g.
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' Q <- fam(c(0.5, 0.3, 0.2))
#' S7::S7_inherits(Q, distribution)
#' @export
distribution <- new_class(
  "distribution",
  abstract = TRUE,
  properties = list(sample_space = space)
)


#' Log density of a distribution
#'
#' The density/mass of a distribution at a point in the sample space. Contrast
#' [compile_loglik()], which takes a family and varies \eqn{\theta}{theta} at
#' fixed outcomes: that is the shape the optimiser needs, and this is the shape
#' a reader needs.
#' @param dist A [distribution].
#' @param x `(M, K)` matrix of outcomes, or a length-`K` vector for one outcome.
#' @return Length-`M` numeric vector.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' log_density(fam(c(0.5, 0.3, 0.2)), c(2L, 1L, 1L))
#' @export
log_density <- new_generic("log_density", "dist", function(dist, x) {
  S7::S7_dispatch()
})


#' Draw outcomes from a distribution
#' @param dist A [distribution].
#' @param n_obs Number of draws.
#' @return `(n_obs, K)` numeric matrix.
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' draw(fam(c(0.5, 0.3, 0.2)), n_obs = 3L)
#' @export
draw <- new_generic("draw", "dist", function(dist, n_obs) {
  S7::S7_dispatch()
})


#' The family a distribution was induced from, if any
#'
#' Returns `NULL` when there is none. \eqn{Q}{Q} is fixed a priori and need not
#' be induced from a family at all, so a family cannot always be recovered.
#' @param object A [distribution].
#' @param ... Ignored.
#' @return A [parametric_family], or `NULL`.
#' @name family.distribution
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' family(fam(c(0.5, 0.3, 0.2)))
NULL


method(family, distribution) <- function(object, ...) NULL


#' The distribution `P_W` induced by a mixing measure `W` and a family
#'
#' The law of `X` when \eqn{\theta \sim W}{theta ~ W} and
#' \eqn{X \sim P_\theta}{X ~ P_theta}. Both the alternative
#' \eqn{Q = P_{W_1}}{Q = P_W1} and the fitted
#' \eqn{\widehat{P}^* = P_{\widehat{W}_0}}{P_star_hat} are of this form.
#'
#' A kernel extends canonically from points to measures, so the degenerate case
#' is not a separate type: a bare parameter vector is taken as a point mass, and
#' `induced_distribution(fam, theta)` is the family at `theta`. Calling a family
#' directly is the shorthand for both, and usually reads better:
#' `fam(theta)` and `fam(W)`.
#'
#' @param family The [parametric_family] whose kernel is pushed forward.
#' @param mixing A [distribution] over its parameter space, or a parameter
#'   vector for a point mass.
#' @return An `induced_distribution`.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' induced_distribution(fam, c(0.5, 0.3, 0.2))
#' induced_distribution(fam, finite_dist(
#'   components = cbind(c(0.6, 0.2, 0.2), c(0.2, 0.6, 0.2)),
#'   weights = c(0.5, 0.5)
#' ))
#' @export
induced_distribution <- new_class(
  "induced_distribution",
  parent = distribution,
  properties = list(family = parametric_family, mixing = distribution),
  constructor = function(family, mixing) {
    if (!S7_inherits(mixing, distribution)) {
      mixing <- dirac(theta = as.numeric(mixing))
    }
    new_object(
      S7_object(),
      sample_space = family@sample_space,
      family = family,
      mixing = mixing
    )
  }
)


method(family, induced_distribution) <- function(object, ...) object@family


method(log_density, induced_distribution) <- function(dist, x) {
  induced_log_density(dist@mixing, dist@family, x)
}


method(draw, induced_distribution) <- function(dist, n_obs) {
  induced_draw(dist@mixing, dist@family, n_obs)
}


#' @rdname induced_distribution
#' @usage NULL
method(print, induced_distribution) <- function(x, ...) {
  cat("<induced_distribution>", format(x), "\n")
  invisible(x)
}


#' @description `format()` names the family and how it was mixed, distinguishing
#'   a point mass from a genuine mixture.
#' @rdname induced_distribution
#' @usage NULL
method(format, induced_distribution) <- function(x, ...) {
  name <- attr(S7_class(x@family), "name")
  n <- n_atoms(x@mixing)
  detail <- if (S7_inherits(x@mixing, dirac)) {
    paste0(
      "at theta = (",
      paste(signif(x@mixing@theta, 3L), collapse = ", "),
      ")"
    )
  } else if (is.na(n)) {
    paste0("mixed over ", attr(S7_class(x@mixing), "name"))
  } else {
    paste0("mixed over ", n, " atoms")
  }
  paste(name, detail)
}


# --- Induced-mixture formulas: double dispatch on (mixing, family) -----------

#' Induced log density `log int dP_theta(x) dW(theta)`
#' @param mixing A [distribution] over the parameter space.
#' @param family A [parametric_family].
#' @param x `(M, K)` matrix of outcomes.
#' @return Length-`M` numeric vector.
#' @keywords internal
induced_log_density <- new_generic(
  "induced_log_density",
  c("mixing", "family"),
  function(mixing, family, x) S7::S7_dispatch()
)


#' Draw from the induced mixture
#' @param mixing A [distribution] over the parameter space.
#' @param family A [parametric_family].
#' @param n_obs Number of draws.
#' @return `(n_obs, K)` numeric matrix.
#' @keywords internal
induced_draw <- new_generic(
  "induced_draw",
  c("mixing", "family"),
  function(mixing, family, n_obs) S7::S7_dispatch()
)


#' @description Sampling from a mixture can be done easily by sampling once
#' from the mixing distribution, then sampling from the family at that parameter
#' value.
#' @rdname induced_draw
#' @usage NULL
method(induced_draw, list(distribution, parametric_family)) <- function(
  mixing,
  family,
  n_obs
) {
  # `draw()` returns rows, so we transpose
  kernel_draw(family, t(draw(mixing, n_obs)))
}
