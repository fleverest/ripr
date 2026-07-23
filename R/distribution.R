#' Distribution over a family's sample space
#'
#' A `distribution` is any distribution over the outcomes of a sampling
#' [family()]: it can evaluate its own log-density over outcomes and (for Monte
#' Carlo engines) draw from itself. It is deliberately role-neutral -- the same
#' object plays the *alternative* Q when passed to [ripr_problem()], and the
#' fitted reverse information projection P* is itself a [mixture_dist()]. A
#' `distribution` knows nothing about nulls, projections, or algorithms.
#' @export
distribution <- new_class("distribution", abstract = TRUE)

#' Log density of a distribution over outcomes
#'
#' Contract: dispatches on the `(distribution, family)` pair; `x` is an `(N, d)`
#' matrix of outcomes, or `NULL` for the family's enumerated support. Returns a
#' length-N numeric vector of `log p(x)`. For mixing-measure distributions this
#' is the induced marginal over outcomes, not the mixing density.
#' @param dist A `distribution`.
#' @param family A `family`.
#' @param x Outcomes, or `NULL` for the enumerated support.
#' @return Numeric vector of log densities.
#' @export
dist_log_density <- new_generic(
  "dist_log_density",
  c("dist", "family"),
  function(dist, family, x = NULL) {
    S7::S7_dispatch()
  }
)

#' Draw outcomes from a distribution
#'
#' Contract: returns an `(n_obs, d)` numeric matrix of draws from the induced
#' outcome distribution. With a non-NULL `seed` the result is a pure function of
#' `(dist, family, n_obs, seed)` and the caller's RNG state is left untouched.
#' This is the common-random-numbers source for Monte Carlo engines.
#' @param dist A `distribution`.
#' @param family A `family`.
#' @param n_obs Number of draws.
#' @param seed Optional integer seed.
#' @return `(n_obs, d)` numeric matrix of draws.
#' @export
dist_sample <- new_generic(
  "dist_sample",
  c("dist", "family"),
  function(dist, family, n_obs, seed = NULL) {
    S7::S7_dispatch()
  }
)

#' Point distribution `p_{theta*}` (the family at a single parameter)
#'
#' @param theta_star Numeric parameter vector.
#' @return A `point_dist`.
#' @export
point_dist <- new_class(
  "point_dist",
  parent = distribution,
  properties = list(theta_star = class_numeric)
)

method(dist_log_density, list(point_dist, family)) <- function(
  dist,
  family,
  x = NULL
) {
  log_density(family, dist@theta_star, x)
}

method(dist_sample, list(point_dist, family)) <- function(
  dist,
  family,
  n_obs,
  seed = NULL
) {
  simulate(family, dist@theta_star, n_obs, seed = seed)
}

#' Finite mixture distribution `sum_c w_c p_{theta_c}`
#'
#' A mixture of the family over parameter atoms. This is both a convenient
#' alternative Q and the shape of the fitted projection P* returned by
#' [run_ripr()]; [prune_mixture()] trims it to its effective support.
#'
#' @param components `(d, C)` numeric matrix; one parameter vector per column.
#' @param weights Length-`C` numeric vector summing to 1.
#' @return A `mixture_dist`.
#' @export
mixture_dist <- new_class(
  "mixture_dist",
  parent = distribution,
  properties = list(
    components = class_any,
    weights = class_numeric
  ),
  validator = function(self) {
    if (!is.matrix(self@components)) {
      return("`components` must be a matrix with one component per column")
    }
    if (ncol(self@components) != length(self@weights)) {
      return("`weights` needs one entry per component column")
    }
    if (abs(sum(self@weights) - 1) > 1e-9) {
      return("`weights` must sum to 1")
    }
    NULL
  }
)

method(dist_log_density, list(mixture_dist, family)) <- function(
  dist,
  family,
  x = NULL
) {
  row_logsumexp(sweep(
    log_density_batch(family, dist@components, x),
    2L,
    log(dist@weights),
    "+"
  ))
}

method(dist_sample, list(mixture_dist, family)) <- function(
  dist,
  family,
  n_obs,
  seed = NULL
) {
  with_rng_seed(seed, {
    idx <- sample.int(
      length(dist@weights),
      n_obs,
      replace = TRUE,
      prob = dist@weights
    )
    do.call(rbind, lapply(idx, function(c_i) {
      simulate(family, dist@components[, c_i], 1L)
    }))
  })
}

#' Plug-in (adaptive) distribution -- interface stub
#'
#' Placeholder for a distribution updated between batches by a supplied
#' `update_fn`. Only the interface exists: resolve it to a concrete distribution
#' before constructing an engine; the density and sampler generics error.
#'
#' @param update_fn Function of `(state, data)` returning a new distribution.
#' @return A `plugin_dist`.
#' @export
plugin_dist <- new_class(
  "plugin_dist",
  parent = distribution,
  properties = list(update_fn = class_function)
)

method(dist_log_density, list(plugin_dist, family)) <- function(
  dist,
  family,
  x = NULL
) {
  stop(
    "plugin_dist is an interface stub: resolve it to a concrete distribution ",
    "(via its update_fn) before constructing an engine"
  )
}

method(dist_sample, list(plugin_dist, family)) <- function(
  dist,
  family,
  n_obs,
  seed = NULL
) {
  stop(
    "plugin_dist is an interface stub: resolve it to a concrete distribution ",
    "(via its update_fn) before sampling"
  )
}

#' Normalise a distribution argument
#'
#' Accepts a `distribution` (returned unchanged) and errors otherwise. The
#' single normalisation seam so problem constructors need not special-case the
#' distribution they are handed.
#' @param x A `distribution`.
#' @return A `distribution`.
#' @keywords internal
as_distribution <- function(x) {
  if (S7_inherits(x, distribution)) {
    return(x)
  }
  stop("cannot interpret this object as a distribution")
}

#' Prune a finite mixture to its effective support
#'
#' Drops the components with weight at or below `threshold` and renormalises the
#' survivors. Components are pruned by weight only, never merged, so several
#' surviving components may still coincide in parameter space.
#'
#' @param x A `mixture_dist`.
#' @param threshold Weight threshold; components with weight `<= threshold` are
#'   dropped. Default `1e-6`.
#' @return A `mixture_dist` over the surviving components.
#' @export
prune_mixture <- new_generic("prune_mixture", "x", function(x, threshold = 1e-6) {
  S7::S7_dispatch()
})

method(prune_mixture, mixture_dist) <- function(x, threshold = 1e-6) {
  keep <- x@weights > threshold
  if (!any(keep)) {
    stop("no component exceeds `threshold`; lower it to keep some support")
  }
  w <- x@weights[keep]
  mixture_dist(
    components = x@components[, keep, drop = FALSE],
    weights = w / sum(w)
  )
}
