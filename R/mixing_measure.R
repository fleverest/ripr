#' Probability distribution (abstract root)
#'
#' Two branches: a [mixing_measure] over a family's *parameter* space, and an
#' [outcome_distribution] over its *sample* space. A [mixture()] bridges them.
#' The `distribution` class has no interface of its own.
#' @examples
#' # `distribution` is abstract; [mixing_measure] and [outcome_distribution]
#' # subclass it, e.g.
#' S7::S7_inherits(point_mixing(theta_star = c(0.5, 0.5)), distribution)
#' @export
distribution <- new_class("distribution", abstract = TRUE)


#' Mixing measure over a family's parameter space
#'
#' A law \eqn{W}{W} over the parameter set of a [parametric_family]: a point mass
#' ([point_mixing()]) or a finite weighted set of atoms ([finite_mixing()]).
#' Paired with a family via [mixture()] it induces the law \eqn{P_W}{P_W} over
#' outcomes.
#'
#' The mixing-measure/mixture distinction follows the mixture-model literature:
#' \eqn{W}{W} is the `mixing_measure`, \eqn{P_W}{P_W} is the [mixture].
#' @references
#'   \insertRef{Lindsay1995}{ripr}
#' @examples
#' # `mixing_measure` is abstract; [point_mixing()] and [finite_mixing()]
#' # subclass it, e.g.
#' S7::S7_inherits(point_mixing(theta_star = c(0.5, 0.5)), mixing_measure)
#' @export
mixing_measure <- new_class(
  "mixing_measure",
  parent = distribution,
  abstract = TRUE
)


#' Point mixing measure: a mass at a single parameter
#'
#' Its induced [mixture] is the family at `theta_star`.
#' @param theta_star Numeric parameter vector.
#' @return A `point_mixing`.
#' @examples
#' point_mixing(theta_star = c(0.4, 0.35, 0.25))
#' @export
point_mixing <- new_class(
  "point_mixing",
  parent = mixing_measure,
  properties = list(theta_star = class_numeric)
)


#' Finite mixing measure `sum_c w_c delta_{theta_c}`
#'
#' A discrete measure on finitely many parameter atoms; the shape of
#' \eqn{\widehat{W}_0}{W0_hat}.
#'
#' @param components `(K, C)` numeric matrix, one parameter vector per column.
#' @param weights Length-`C` numeric vector summing to 1.
#' @return A `finite_mixing`.
#' @examples
#' finite_mixing(
#'   components = cbind(c(0.6, 0.2, 0.2), c(0.2, 0.6, 0.2)),
#'   weights = c(0.5, 0.5)
#' )
#' @export
finite_mixing <- new_class(
  "finite_mixing",
  parent = mixing_measure,
  properties = list(components = class_any, weights = class_numeric),
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


#' Number of atoms in a mixing measure
#' @param x A [mixing_measure].
#' @return Integer.
#' @examples
#' n_atoms(point_mixing(theta_star = c(0.5, 0.5)))
#' n_atoms(finite_mixing(components = cbind(c(0.6, 0.4), c(0.2, 0.8)), weights = c(0.5, 0.5)))
#' @export
n_atoms <- new_generic("n_atoms", "x", function(x) S7::S7_dispatch())


method(n_atoms, point_mixing) <- function(x) 1L


method(n_atoms, finite_mixing) <- function(x) ncol(x@components)


#' Parameter atoms of a mixing measure
#' @param x A [mixing_measure].
#' @return `(K, C)` numeric matrix.
#' @examples
#' atoms(finite_mixing(components = cbind(c(0.6, 0.4), c(0.2, 0.8)), weights = c(0.5, 0.5)))
#' @export
atoms <- new_generic("atoms", "x", function(x) S7::S7_dispatch())


method(atoms, point_mixing) <- function(x) matrix(x@theta_star, ncol = 1L)


method(atoms, finite_mixing) <- function(x) x@components


#' Weights of a mixing measure
#' @param object A [mixing_measure].
#' @param ... Ignored.
#' @return Numeric vector summing to 1.
#' @name weights.mixing_measure
#' @examples
#' weights(finite_mixing(components = cbind(c(0.5, 0.5)), weights = 1))
NULL


method(weights, point_mixing) <- function(object, ...) 1


method(weights, finite_mixing) <- function(object, ...) object@weights


#' Drop atoms below a weight threshold and renormalise
#'
#' Atoms are pruned by weight only, never merged, so survivors may still
#' coincide in parameter space.
#'
#' @param x A [finite_mixing].
#' @param threshold Atoms with weight `<= threshold` are dropped.
#' @return A [finite_mixing] over the survivors.
#' @examples
#' w <- finite_mixing(
#'   components = cbind(c(0.6, 0.4), c(0.2, 0.8), c(0.5, 0.5)),
#'   weights = c(0.98, 0.01, 0.01)
#' )
#' prune(w, threshold = 0.05)
#' @export
prune <- new_generic("prune", "x", function(x, threshold = 1e-8) {
  S7::S7_dispatch()
})


method(prune, finite_mixing) <- function(x, threshold = 1e-8) {
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
  finite_mixing(
    components = x@components[, keep, drop = FALSE],
    weights = w / sum(w)
  )
}


#' Modal parameter of a mixing measure
#'
#' The point the measure places the highest mass or density. Used to seed the
#' starting point for the RIPr optimiser.
#' @param x A [mixing_measure].
#' @return Numeric vector of a parameter's length.
#' @examples
#' w <- finite_mixing(
#'   components = cbind(c(0.6, 0.4), c(0.2, 0.8)),
#'   weights = c(0.3, 0.7)
#' )
#' mode_parameter(w)
#' @export
mode_parameter <- new_generic("mode_parameter", "x", \(x) S7::S7_dispatch())


method(mode_parameter, point_mixing) <- function(x) x@theta_star


method(mode_parameter, finite_mixing) <- function(x) {
  x@components[, which.max(x@weights)]
}
