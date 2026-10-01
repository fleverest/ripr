#' @include distribution.R region.R
NULL


# Mixing measures: distributions over parameter spaces that induce mixtures.

#' Distributions with finite support
#'
#' A [distribution] on finitely many weighted atoms. [finite_dist()] is the
#' standard instance, [dirac()] builds a one-atom [finite_dist()].
#'
#' @examples
#' S7::S7_inherits(
#'   finite_dist(atoms = rbind(c(0.6, 0.4)), weights = 1),
#'   discrete_dist
#' )
#' @param sample_space The [space] this is a law over. Derived by concrete
#'   subclasses, never passed by a caller.
#' @export
discrete_dist <- new_class(
  "discrete_dist",
  parent = distribution,
  abstract = TRUE
)


#' A distribution on finitely many atoms, `sum_c w_c delta_{theta_c}`
#'
#' A discrete measure on finitely many parameter atoms; the shape of
#' \eqn{\widehat{W}_0}{W0_hat}.
#'
#' @param atoms `(C, d)` numeric matrix, one parameter vector (atom) per row.
#' @param weights Length-`C` numeric vector summing to 1.
#' @return A `finite_dist`.
#' @examples
#' d <- finite_dist(
#'   atoms = rbind(c(0.6, 0.2, 0.2), c(0.2, 0.6, 0.2)),
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
    atoms = class_any,
    weights = class_numeric,
    sample_space = new_property(
      space,
      getter = function(self) real_region(ncol(self@atoms))
    )
  ),
  validator = function(self) {
    if (!is.matrix(self@atoms)) {
      return("`atoms` must be a matrix with one atom per row")
    }
    if (nrow(self@atoms) != length(self@weights)) {
      return("`weights` needs one entry per row of `atoms`")
    }
    if (any(self@weights < 0)) {
      return("`weights` must be non-negative")
    }
    if (abs(sum(self@weights) - 1) > 1e-8) {
      return("`weights` must sum to 1")
    }
    NULL
  }
)


#' A point mass at a single value
#'
#' A degenerate mixing measure such that `fam(dirac(theta))` is equivalent to
#' `fam(theta)`. It is a [finite_dist] with one atom, not a class of its own.
#' Drawing samples from [dirac()] does not change the state of the pseudorandom
#' number generator.
#'
#' @param theta Numeric parameter vector.
#' @return A [finite_dist] with one atom, at `theta`, of weight 1.
#' @examples
#' d <- dirac(theta = c(0.4, 0.35, 0.25))
#' d
#' atoms(d)
#' @export
dirac <- function(theta) {
  if (!is.numeric(theta) || length(theta) == 0L || anyNA(theta)) {
    stop(
      "`theta` must be a non-empty numeric vector without NAs.",
      call. = FALSE
    )
  }
  finite_dist(atoms = matrix(as.numeric(theta), nrow = 1L), weights = 1)
}


#' @rdname finite_dist
#' @usage NULL
method(format, finite_dist) <- function(x, ...) {
  sprintf(
    "finite_dist: %s in R^%d",
    count_label(nrow(x@atoms), "atom"),
    ncol(x@atoms)
  )
}


#' @description `print()` shows the atoms and their weights as a table when
#'   there are at most eight, and the atom count with the heaviest atom
#'   otherwise; [atoms()] and [weights()] give the full support either way.
#' @rdname finite_dist
#' @usage NULL
method(print, finite_dist) <- function(x, ...) {
  cat("<finite_dist>\n")
  cat(
    "  ",
    count_label(nrow(x@atoms), "atom"),
    " in R^",
    ncol(x@atoms),
    "\n",
    sep = ""
  )
  if (nrow(x@atoms) <= 8L) {
    m <- x@atoms
    if (is.null(colnames(m))) {
      colnames(m) <- paste0("theta", seq_len(ncol(m)))
    }
    print(signif(cbind(m, weight = x@weights), 4L))
  } else {
    i <- which.max(x@weights)
    cat(
      "  heaviest atom ",
      theta_label(x@atoms[i, ]),
      " with weight ",
      signif(x@weights[i], 3L),
      "\n",
      sep = ""
    )
  }
  invisible(x)
}


#' Continuous distributions
#'
#' A [distribution] with a density rather than atoms and weights. See for
#' instance [gaussian_dist()], [dirichlet()] and [truncated_dirichlet()].
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


#' @description A [discrete_dist] is supported in `space` when all its atoms are.
#' @rdname supported_in
#' @usage NULL
method(supported_in, discrete_dist) <- function(dist, space) {
  tryCatch(
    all(apply(atoms(dist), 1L, function(theta) contains(space, theta))),
    error = function(e) FALSE
  )
}


#' Number of atoms in a distribution
#' @param x A [distribution] over a parameter space.
#' @return Integer.
#' @examples
#' n_atoms(dirac(theta = c(0.5, 0.5)))
#' n_atoms(finite_dist(atoms = rbind(c(0.6, 0.4), c(0.2, 0.8)), weights = c(0.5, 0.5)))
#' @export
n_atoms <- new_generic("n_atoms", "x", function(x) S7::S7_dispatch())


method(n_atoms, finite_dist) <- function(x) nrow(x@atoms)


#' @description `NA` for a continuous measure.
#' @rdname n_atoms
#' @usage NULL
method(n_atoms, continuous_dist) <- function(x) NA_integer_


#' Atoms of a distribution
#' @param x A [distribution] over a parameter space.
#' @return `(C, d)` numeric matrix, one atom per row.
#' @examples
#' atoms(finite_dist(atoms = rbind(c(0.6, 0.4), c(0.2, 0.8)), weights = c(0.5, 0.5)))
#' @export
atoms <- new_generic("atoms", "x", function(x) S7::S7_dispatch())


method(atoms, finite_dist) <- function(x) x@atoms


#' The error [atoms()] and [weights()] raise for a continuous measure.
#' @keywords internal
#' @noRd
refuse_continuous <- function(x, what) {
  stop(
    "`",
    what,
    "()` is not defined for a `",
    class_name(x),
    "`: a continuous distribution has a density rather than a support to ",
    "list. Use `draw()` to sample it, or `reference_point()` for the ",
    "point it concentrates on.",
    call. = FALSE
  )
}


#' @rdname atoms
#' @usage NULL
method(atoms, continuous_dist) <- function(x) refuse_continuous(x, "atoms")


#' @description A finite distribution draws its atoms with probability equal to
#'   their weights; a [dirac()] does so without changing the state of the
#'    pseudorandom number generator.
#' @rdname draw
#' @usage NULL
method(draw, finite_dist) <- function(dist, n) {
  if (nrow(dist@atoms) == 1L) {
    return(matrix(
      dist@atoms[1L, ],
      nrow = n,
      ncol = ncol(dist@atoms),
      byrow = TRUE
    ))
  }
  idx <- sample.int(
    length(dist@weights),
    n,
    replace = TRUE,
    prob = dist@weights
  )
  dist@atoms[idx, , drop = FALSE]
}


#' Weights of a distribution
#' @param object A [distribution] over a parameter space.
#' @param ... Ignored.
#' @return Numeric vector summing to 1.
#' @name weights.distribution
#' @examples
#' weights(finite_dist(atoms = rbind(c(0.5, 0.5)), weights = 1))
NULL


method(weights, finite_dist) <- function(object, ...) object@weights


#' @rdname weights.distribution
#' @usage NULL
method(weights, continuous_dist) <- function(object, ...) {
  refuse_continuous(object, "weights")
}


#' Replace a distribution with the empirical distribution of draws from it
#'
#' Results in a [finite_dist] of `n` equally weighted draws from `dist`,
#' converging to `dist` as `n` grows. Useful for approximating a mixture that
#' has no closed form, or for approximating a distribution for which it is
#' difficult to compute expectations under.
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
    atoms = draw(dist, n),
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
#'   atoms = rbind(c(0.6, 0.4), c(0.2, 0.8), c(0.5, 0.5)),
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
    atoms = x@atoms[keep, , drop = FALSE],
    weights = w / sum(w)
  )
}


#' A point in the support of `x`, used only to seed the optimiser's starting
#' atoms. A family, when the alternative is not a mixture, answers with its
#' parameter space's point nearest the origin.
#' @keywords internal
#' @noRd
reference_point <- new_generic("reference_point", "x", \(x) S7::S7_dispatch())


method(reference_point, parametric_family) <- function(x) {
  family <- x
  space <- family@parameter_space
  project(space, rep(0, space_dim(space)))
}


method(reference_point, finite_dist) <- function(x) {
  x@atoms[which.max(x@weights), ]
}


method(mixture_log_density, list(finite_dist, parametric_family)) <- function(
  mixing,
  family,
  x
) {
  row_logsumexp(add_by_col(
    kernel_loglik_batch(family, mixing@atoms, x),
    log(mixing@weights)
  ))
}


method(mixture_log_density, list(continuous_dist, parametric_family)) <-
  function(mixing, family, x) {
    stop(
      "no induced density is implemented for a `",
      class_name(mixing),
      "` over a `",
      class_name(family),
      "`. Mixing a continuous measure through a kernel is an integral, and ",
      "only some pairings have one in closed or quadrature form.\n",
      "This is not approximated by Monte Carlo by default, see ",
      "`discretise(mixing, n)` if you would like to approximate the ",
      "mixture by sampling atoms for a mixing measure.",
      call. = FALSE
    )
  }
