#' Alternative hypothesis distribution Q
#'
#' An `alternative` specifies the numerator distribution the RIPr projects: its
#' density over outcomes (given a sampling family) and a sampler for Monte Carlo
#' engines. Alternatives know nothing about nulls or algorithms.
#' @export
alternative <- new_class("alternative", abstract = TRUE)

#' Log density of Q over outcomes
#'
#' Contract: dispatches on the `(alternative, family)` pair; `x` is an `(N, d)`
#' matrix of outcomes, or `NULL` for the family's enumerated support. Returns a
#' length-N numeric vector of `log q(x)`. For mixing-measure alternatives this
#' is the induced marginal over outcomes, not the mixing density.
#' @param alt An `alternative`.
#' @param family A `family`.
#' @param x Outcomes, or `NULL` for the enumerated support.
#' @return Numeric vector of log densities.
#' @export
q_log_density <- new_generic(
  "q_log_density",
  c("alt", "family"),
  function(alt, family, x = NULL) {
    S7::S7_dispatch()
  }
)

#' Draw outcomes from Q
#'
#' Contract: returns an `(n_obs, d)` numeric matrix of draws from the induced
#' outcome distribution. With a non-NULL `seed` the result is a pure function of
#' `(alt, family, n_obs, seed)` and the caller's RNG state is left untouched.
#' This is the common-random-numbers source for Monte Carlo engines.
#' @param alt An `alternative`.
#' @param family A `family`.
#' @param n_obs Number of draws.
#' @param seed Optional integer seed.
#' @return `(n_obs, d)` numeric matrix of draws.
#' @export
q_sample <- new_generic(
  "q_sample",
  c("alt", "family"),
  function(alt, family, n_obs, seed = NULL) {
    S7::S7_dispatch()
  }
)

#' Point alternative `Q = p_{theta*}`
#'
#' @param theta_star Numeric parameter vector of the single alternative.
#' @return A `point_alt`.
#' @export
point_alt <- new_class(
  "point_alt",
  parent = alternative,
  properties = list(theta_star = class_numeric)
)

method(q_log_density, list(point_alt, family)) <- function(
  alt,
  family,
  x = NULL
) {
  log_density(family, alt@theta_star, x)
}

method(q_sample, list(point_alt, family)) <- function(
  alt,
  family,
  n_obs,
  seed = NULL
) {
  simulate(family, alt@theta_star, n_obs, seed = seed)
}

#' Finite mixture alternative `Q = sum_c w_c p_{theta_c}`
#'
#' @param components `(d, C)` numeric matrix; one parameter vector per column.
#' @param weights Length-`C` numeric vector summing to 1.
#' @return A `mixture_alt`.
#' @export
mixture_alt <- new_class(
  "mixture_alt",
  parent = alternative,
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

method(q_log_density, list(mixture_alt, family)) <- function(
  alt,
  family,
  x = NULL
) {
  row_logsumexp(sweep(
    log_density_batch(family, alt@components, x),
    2L,
    log(alt@weights),
    "+"
  ))
}

method(q_sample, list(mixture_alt, family)) <- function(
  alt,
  family,
  n_obs,
  seed = NULL
) {
  with_rng_seed(seed, {
    idx <- sample.int(
      length(alt@weights),
      n_obs,
      replace = TRUE,
      prob = alt@weights
    )
    do.call(rbind, lapply(idx, function(c_i) {
      simulate(family, alt@components[, c_i], 1L)
    }))
  })
}

#' Plug-in (adaptive) alternative -- interface stub
#'
#' Placeholder for alternatives updated between batches by a supplied
#' `update_fn`. Only the interface exists: resolve it to a concrete alternative
#' before constructing an engine; the density and sampler generics error.
#'
#' @param update_fn Function of `(state, data)` returning a new alternative.
#' @return A `plugin_alt`.
#' @export
plugin_alt <- new_class(
  "plugin_alt",
  parent = alternative,
  properties = list(update_fn = class_function)
)

method(q_log_density, list(plugin_alt, family)) <- function(
  alt,
  family,
  x = NULL
) {
  stop(
    "plugin_alt is an interface stub: resolve it to a concrete alternative ",
    "(via its update_fn) before constructing an engine"
  )
}

method(q_sample, list(plugin_alt, family)) <- function(
  alt,
  family,
  n_obs,
  seed = NULL
) {
  stop(
    "plugin_alt is an interface stub: resolve it to a concrete alternative ",
    "(via its update_fn) before sampling"
  )
}

#' Normalise an alternative argument
#'
#' Accepts an `alternative` (returned unchanged) and errors otherwise. The
#' single normalisation seam so problem constructors need not special-case
#' their `alternative` argument.
#' @param x An `alternative`.
#' @return An `alternative`.
#' @keywords internal
as_alternative <- function(x) {
  if (S7_inherits(x, alternative)) {
    return(x)
  }
  stop("cannot interpret this object as an alternative distribution")
}
