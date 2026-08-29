#' @include polyhedron_region.R
NULL

# --- Real region --------------------------------------------------------------

#' The whole of `R^d` as a region
#'
#' All of \eqn{\mathbb{R}^d}{R^d}: \eqn{\Theta = \mathbb{R}^d}{Theta = R^d}, with the
#' identity chart. This is the parameter space of a [gaussian_family()], and the
#' way to say that a null places no constraint at all.
#'
#' Being unbounded it has no vertices, so like [halfspace_region()] it admits no
#' certified gap bound.
#'
#' @param d Integer dimension.
#' @return A `real_region`.
#' @examples
#' real_region(2L)
#' project(real_region(2L), c(3, -1))
#' @export
real_region <- new_class(
  "real_region",
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


method(project, real_region) <- function(space, theta) {
  as.numeric(theta)
}


method(contains, real_region) <- function(space, theta, tol = 1e-8) {
  length(theta) == as.integer(space@n_dim) && all(is.finite(theta))
}


method(region_phrase, real_region) <- function(space) {
  sprintf("all of R^%d", space_dim(space))
}


#' @description All of \eqn{\mathbb{R}^d}{R^d} is admissible, so the inherited
#'   checks (numeric, correct shape, no missing values) are almost all that is
#'   needed. Only finiteness has to be added.
#' @rdname validate_outcome
#' @usage NULL
method(validate_outcome, real_region) <- function(space, x) {
  x <- check_outcome_shape(x, space_dim(space))
  if (any(!is.finite(x))) {
    stop("outcomes must be finite.", call. = FALSE)
  }
  x
}
