#' @include polyhedron_region.R
NULL

# --- Unconstrained region -----------------------------------------------------

#' The whole of `R^d` as a region
#'
#' The unconstrained case: \eqn{\Theta = \mathbb{R}^d}{Theta = R^d}, with the
#' identity chart. This is the parameter space of a [gaussian_family()], and the
#' way to say that a null places no constraint at all.
#'
#' Being unbounded it has no vertices, so like [halfspace_region()] it admits no
#' certified gap bound.
#'
#' The name says what the region is, rather than what it is made of: a
#' [real_space] is a *sample* space, and the two were too easy to confuse while
#' this one carried the same `real_` prefix.
#'
#' @param d Integer dimension.
#' @return An `unconstrained_region`.
#' @examples
#' unconstrained_region(2L)
#' project(unconstrained_region(2L), c(3, -1))
#' @export
unconstrained_region <- new_class(
  "unconstrained_region",
  parent = polyhedron_region,
  properties = list(
    n_dim = new_property(
      class_numeric,
      getter = function(self) nrow(self@generators$v)
    )
  ),
  constructor = function(d) {
    d <- as.integer(d)
    stopifnot(
      "`d` must be a single positive integer" = length(d) == 1L &&
        !is.na(d) &&
        d >= 1L
    )
    new_object(
      polyhedron_region(
        .hv = hv_from_h(
          # A caller reading facets wants to be told there are none.
          # `q_hrep()` states the same region to cddlib as the trivially true
          # row `0 . x <= 1`, via `as_hmatrix()`'s zero-row branch, because
          # cddlib needs a row to work with.
          h = list(
            a = matrix(numeric(0), nrow = 0L, ncol = d),
            b = numeric(0),
            eq = logical(0)
          ),
          v = make_generators(matrix(0, nrow = d, ncol = 1L), NULL, diag(d))
        )
      )
    )
  }
)


method(project, unconstrained_region) <- function(space, theta) {
  as.numeric(theta)
}


method(contains, unconstrained_region) <- function(space, theta, tol = 1e-8) {
  length(theta) == as.integer(space@n_dim) && all(is.finite(theta))
}


method(region_phrase, unconstrained_region) <- function(space) {
  sprintf("all of R^%d", space_dim(space))
}
