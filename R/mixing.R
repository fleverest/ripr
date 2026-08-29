#' @include distribution.R region.R
NULL


# Distributions for producing mixtures
#
# This file is an extension of distribution.R, including distributions for
# producing useful mixtures.

#' Distributions with finite support
#'
#' A [distribution] defined by finitely many atoms with probabilities
#' represented by weights on them.  [dirac()] and [finite_dist()] are the two
#' instances.
#'
#' @examples
#' S7::S7_inherits(dirac(theta = c(0.5, 0.5)), discrete_dist)
#' S7::S7_inherits(
#'   finite_dist(components = cbind(c(0.6, 0.4)), weights = 1),
#'   discrete_dist
#' )
#' @param sample_space The [space] this is a law over. Inherited from
#'   [distribution]; a concrete subclass derives it rather than taking it, so
#'   it is never passed by a caller.
#' @export
discrete_dist <- new_class(
  "discrete_dist",
  parent = distribution,
  abstract = TRUE
)


#' A point mass at a single parameter value
#'
#' The degenerate [distribution]: all of its mass at `theta`. Used as a mixing
#' measure it is the case where no mixing happens at all, so `fam(dirac(theta))`
#' is equivalent to `fam(theta)`.
#' @param theta Numeric parameter vector.
#' @return A `dirac`.
#' @examples
#' dirac(theta = c(0.4, 0.35, 0.25))
#' @export
dirac <- new_class(
  "dirac",
  parent = discrete_dist,
  properties = list(
    theta = class_numeric,
    sample_space = new_property(
      space,
      getter = function(self) real_region(length(self@theta))
    )
  )
)


#' A distribution on finitely many atoms, `sum_c w_c delta_{theta_c}`
#'
#' A discrete measure on finitely many parameter atoms; the shape of
#' \eqn{\widehat{W}_0}{W0_hat}.
#'
#' @param components `(K, C)` numeric matrix, one parameter vector per column.
#' @param weights Length-`C` numeric vector summing to 1.
#' @return A `finite_dist`.
#' @examples
#' d <- finite_dist(
#'   components = cbind(c(0.6, 0.2, 0.2), c(0.2, 0.6, 0.2)),
#'   weights = c(0.5, 0.5)
#' )
#' n_atoms(d)
#' atoms(d)
#' weights(d)
#' @export
finite_dist <- new_class(
  "finite_dist",
  parent = discrete_dist,
  properties = list(
    components = class_any,
    weights = class_numeric,
    sample_space = new_property(
      space,
      getter = function(self) real_region(nrow(self@components))
    )
  ),
  validator = function(self) {
    if (!is.matrix(self@components)) {
      return("`components` must be a matrix with one component per column")
    }
    if (ncol(self@components) != length(self@weights)) {
      return("`weights` needs one entry per column of `components`")
    }
    if (any(self@weights < 0)) {
      return("`weights` must be non-negative")
    }
    if (abs(sum(self@weights) - 1) > 1e-9) {
      return("`weights` must sum to 1")
    }
    NULL
  }
)


#' Continuous distributions
#'
#' A [distribution] that is continuous rather than discrete, so it has a
#' density whereas [finite_dist]s have atoms and weights. See [gaussian_dist()],
#' [dirichlet()] and [truncated_dirichlet()].
#'
#' @examples
#' # `continuous_dist` is abstract; dirichlet() subclasses it, e.g.
#' S7::S7_inherits(dirichlet(alpha = c(2, 1, 1)), continuous_dist)
#' n_atoms(dirichlet(alpha = c(2, 1, 1)))
#' @param sample_space The [space] this is a law over. Inherited from
#'   [distribution].
#' @export
continuous_dist <- new_class(
  "continuous_dist",
  parent = distribution,
  abstract = TRUE
)


#' @description Check that a [discrete_dist] is supported in a space. Returns
#' `TRUE` is all of `dist`s atoms are contained in `space`.
#' @rdname supported_in
#' @usage NULL
method(supported_in, discrete_dist) <- function(dist, space) {
  tryCatch(
    all(apply(atoms(dist), 2L, function(theta) contains(space, theta))),
    error = function(e) FALSE
  )
}


#' Number of atoms in a distribution
#' @param x A [distribution] over a parameter space.
#' @return Integer.
#' @examples
#' n_atoms(dirac(theta = c(0.5, 0.5)))
#' n_atoms(finite_dist(components = cbind(c(0.6, 0.4), c(0.2, 0.8)), weights = c(0.5, 0.5)))
#' @export
n_atoms <- new_generic("n_atoms", "x", function(x) S7::S7_dispatch())


method(n_atoms, dirac) <- function(x) 1L


method(n_atoms, finite_dist) <- function(x) ncol(x@components)


#' @description A continuous measure has no atoms to count, which is `NA`
#'   rather than `0`: zero would say the measure was empty.
#' @rdname n_atoms
#' @usage NULL
method(n_atoms, continuous_dist) <- function(x) NA_integer_


#' Atoms of a distribution
#' @param x A [distribution] over a parameter space.
#' @return `(K, C)` numeric matrix.
#' @examples
#' atoms(finite_dist(components = cbind(c(0.6, 0.4), c(0.2, 0.8)), weights = c(0.5, 0.5)))
#' @export
atoms <- new_generic("atoms", "x", function(x) S7::S7_dispatch())


method(atoms, dirac) <- function(x) matrix(x@theta, ncol = 1L)


method(atoms, finite_dist) <- function(x) x@components


#' The refusal a continuous measure owes both support accessors
#' @keywords internal
#' @noRd
refuse_continuous <- function(x, what) {
  stop(
    "`",
    what,
    "()` is not defined for a `",
    attr(S7_class(x), "name"),
    "`: a continuous distribution has a density rather than a support to ",
    "list. Use `draw()` to sample it, or `reference_point()` for the ",
    "point it concentrates on.",
    call. = FALSE
  )
}


#' @rdname atoms
#' @usage NULL
method(atoms, continuous_dist) <- function(x) refuse_continuous(x, "atoms")


#' @description A point mass draws the same parameter every time.
#' @rdname draw
#' @usage NULL
method(draw, dirac) <- function(dist, n_obs) {
  matrix(dist@theta, nrow = n_obs, ncol = length(dist@theta), byrow = TRUE)
}


#' @description A finite distribution draws its atoms with probability equal to
#'   their weights. Repeats stay in place rather than being grouped:
#'   [kernel_draw()] is vectorised over parameters, so a repeated row costs
#'   nothing.
#' @rdname draw
#' @usage NULL
method(draw, finite_dist) <- function(dist, n_obs) {
  idx <- sample.int(
    length(dist@weights),
    n_obs,
    replace = TRUE,
    prob = dist@weights
  )
  t(dist@components[, idx, drop = FALSE])
}


#' Weights of a distribution
#' @param object A [distribution] over a parameter space.
#' @param ... Ignored.
#' @return Numeric vector summing to 1.
#' @name weights.distribution
#' @examples
#' weights(finite_dist(components = cbind(c(0.5, 0.5)), weights = 1))
NULL


method(weights, dirac) <- function(object, ...) 1


method(weights, finite_dist) <- function(object, ...) object@weights


#' @rdname weights.distribution
#' @usage NULL
method(weights, continuous_dist) <- function(object, ...) {
  refuse_continuous(object, "weights")
}


#' Replace a distribution with a finite one comprised of atoms drawn from it
#'
#' Constructs a [finite_dist] comprising `n` draws from `dist` as equally
#' weighted atoms. It is the empirical distribution based on a sample drawn
#' from it. Converges to `dist` as `n` grows.
#'
#' This can be used as a convenient way to approximate a mixture distribution
#' if no closed-form exists. See examples.
#'
#' @param dist A [distribution] to sample.
#' @param n Number of draws.
#' @return A [finite_dist] on `n` equally weighted atoms.
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 6L, k = 3L)
#'
#' W <- dirichlet(alpha = c(4, 3, 2))
#' approx <- discretise(W, 500L)
#' n_atoms(approx)
#'
#' x <- enumerate_space(fam@sample_space)
#' log_density(fam(W), c(2,2,2))
#' log_density(fam(approx), c(2,2,2))
#' @export
discretise <- function(dist, n) {
  rlang::check_number_whole(n, min = 1, max = 2147483647)
  n <- as.integer(n)
  finite_dist(
    components = t(draw(dist, n)),
    weights = rep(1 / n, n)
  )
}


#' Drop atoms below a weight threshold and renormalise
#'
#' Atoms are pruned by weight only, never merged, so survivors may still
#' coincide in parameter space.
#'
#' @param x A [finite_dist].
#' @param threshold Atoms with weight `<= threshold` are dropped.
#' @return A [finite_dist] over the survivors.
#' @examples
#' w <- finite_dist(
#'   components = cbind(c(0.6, 0.4), c(0.2, 0.8), c(0.5, 0.5)),
#'   weights = c(0.98, 0.01, 0.01)
#' )
#' prune(w, threshold = 0.05)
#' @export
prune <- new_generic("prune", "x", function(x, threshold = 1e-8) {
  S7::S7_dispatch()
})


method(prune, finite_dist) <- function(x, threshold = 1e-8) {
  keep <- x@weights > threshold
  if (!any(keep)) {
    stop(
      "no atom has weight above `threshold` (",
      threshold,
      "); the largest is ",
      signif(max(x@weights), 3),
      ".",
      call. = FALSE
    )
  }
  w <- x@weights[keep]
  finite_dist(
    components = x@components[, keep, drop = FALSE],
    weights = w / sum(w)
  )
}


#' A representative point of a distribution or a parameter space
#'
#' One point that stands for the whole, used to seed the RIPr optimiser's
#' starting atoms. The guarantee is that it lies **in the support**: it is a
#' point the optimiser may legally start from.
#'
#' It is deliberately not called a mode. For most of these it is one -- the
#' heaviest atom of a [finite_dist], the mean of a [gaussian_dist] -- but a
#' [dirichlet()] with any concentration at or below 1 has its mode on the
#' boundary or not at all, and falls back to the mean. Naming it for the mode
#' would make the name a lie in exactly the case where the fallback matters.
#'
#' Defined for a [parametric_family] as well: when the alternative is not a
#' mixture there is nothing to take a point from, and the family's own space
#' answers instead with the point closest to the origin -- the centroid of a
#' simplex, the origin itself for an unconstrained space.
#' @param x A [distribution] over a parameter space, or a [parametric_family].
#' @return Numeric vector of a parameter's length, inside the support.
#' @examples
#' w <- finite_dist(
#'   components = cbind(c(0.6, 0.4), c(0.2, 0.8)),
#'   weights = c(0.3, 0.7)
#' )
#' reference_point(w)
#' reference_point(multinomial_family(n_trials = 4L, k = 3L))
#' @export
reference_point <- new_generic("reference_point", "x", \(x) S7::S7_dispatch())


#' @description A family answers with the point of its parameter space closest
#'   to the origin: the centroid of a simplex, the origin itself for an
#'   unconstrained space. This is the fallback when the alternative is not a
#'   mixture and so has no point of its own to offer.
#' @rdname reference_point
#' @usage NULL
method(reference_point, parametric_family) <- function(x) {
  family <- x
  space <- family@parameter_space
  project(space, rep(0, space_dim(space)))
}


method(reference_point, dirac) <- function(x) x@theta


method(reference_point, finite_dist) <- function(x) {
  x@components[, which.max(x@weights)]
}


method(mixture_log_density, list(dirac, parametric_family)) <- function(
  mixing,
  family,
  x
) {
  kernel_loglik(family, mixing@theta, x)
}


method(mixture_log_density, list(finite_dist, parametric_family)) <- function(
  mixing,
  family,
  x
) {
  row_logsumexp(add_by_col(
    kernel_loglik_batch(family, mixing@components, x),
    log(mixing@weights)
  ))
}


method(mixture_log_density, list(continuous_dist, parametric_family)) <-
  function(mixing, family, x) {
    stop(
      "no induced density is implemented for a `",
      attr(S7_class(mixing), "name"),
      "` over a `",
      attr(S7_class(family), "name"),
      "`. Mixing a continuous measure through a kernel is an integral, and ",
      "only some pairings have one in closed or quadrature form.\n",
      "This is not approximated by Monte Carlo by default, see ",
      "`discretise(mixing, n)` if you would like to approximate the ",
      "mixture by sampling atoms for a mixing measure.",
      call. = FALSE
    )
  }
