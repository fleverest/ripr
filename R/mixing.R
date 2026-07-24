#' Mixing measure over a family's parameter space
#'
#' A `mixing` is the parameter-space branch of [distribution]: a law over the
#' parameter set Theta of a sampling [family()] -- a point mass ([point_mixing()]),
#' a finite set of weighted atoms ([finite_mixing()]), or a continuous prior
#' (e.g. [gaussian_mixing()]). On its own it is a law over *parameters* and knows
#' nothing of a sample space; paired with a family via [as_marginal()] it induces
#' a [marginal] over *outcomes* -- the role it plays as the alternative Q and as
#' the fitted reverse information projection P*. A `mixing` knows nothing about
#' nulls, projections, or algorithms.
#' @export
mixing <- new_class("mixing", parent = distribution, abstract = TRUE)

#' Point mixing measure: a mass at a single parameter `theta*`
#'
#' Its induced marginal is the family at `theta*`, `p_{theta*}`.
#' @param theta_star Numeric parameter vector.
#' @return A `point_mixing`.
#' @export
point_mixing <- new_class(
  "point_mixing",
  parent = mixing,
  properties = list(theta_star = class_numeric)
)

#' Finite mixing measure `sum_c w_c delta_{theta_c}`
#'
#' A discrete measure on a finite set of parameter atoms. Its induced marginal is
#' the mixture `sum_c w_c p_{theta_c}` -- both a convenient alternative Q and the
#' shape of the fitted projection P* returned by [run_ripr()]. [prune_mixture()]
#' trims it to its effective support.
#'
#' @param components `(d, C)` numeric matrix; one parameter vector per column.
#' @param weights Length-`C` numeric vector summing to 1.
#' @return A `finite_mixing`.
#' @export
finite_mixing <- new_class(
  "finite_mixing",
  parent = mixing,
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

#' Normalise a mixing-measure argument
#'
#' Accepts a `mixing` (returned unchanged) and errors otherwise. The single
#' normalisation seam so problem constructors need not special-case the mixing
#' measure they are handed.
#' @param x A `mixing`.
#' @return A `mixing`.
#' @keywords internal
as_mixing <- function(x) {
  if (S7_inherits(x, mixing)) {
    return(x)
  }
  stop("cannot interpret this object as a mixing measure")
}

#' Prune a finite mixing measure to its effective support
#'
#' Drops the atoms with weight at or below `threshold` and renormalises the
#' survivors. Atoms are pruned by weight only, never merged, so several surviving
#' atoms may still coincide in parameter space.
#'
#' @param x A `finite_mixing`.
#' @param threshold Weight threshold; atoms with weight `<= threshold` are
#'   dropped. Default `1e-6`.
#' @return A `finite_mixing` over the surviving atoms.
#' @export
prune_mixture <- new_generic(
  "prune_mixture",
  "x",
  function(x, threshold = 1e-6) {
    S7::S7_dispatch()
  }
)

method(prune_mixture, finite_mixing) <- function(x, threshold = 1e-6) {
  keep <- x@weights > threshold
  if (!any(keep)) {
    stop("no atom exceeds `threshold`; lower it to keep some support")
  }
  w <- x@weights[keep]
  finite_mixing(
    components = x@components[, keep, drop = FALSE],
    weights = w / sum(w)
  )
}
