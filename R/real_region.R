#' @include polyhedron_region.R
NULL

# --- Real region --------------------------------------------------------------

#' The whole of `R^d` as a region
#'
#' \eqn{\Theta = \mathbb{R}^d}{Theta = R^d} with the identity chart: the
#' parameter space of a [gaussian_family()], or an unconstrained null. Being
#' unbounded, it cannot be certified.
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
    d = new_property(
      class_numeric,
      getter = function(self) ncol(self@generators$v)
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
        .hv = hv_fill(
          # No facets; `as_hmatrix()` gives cddlib the trivial row `0 . x <= 1`
          # since it needs at least one.
          h = list(
            a = matrix(numeric(0), nrow = 0L, ncol = d),
            b = numeric(0),
            eq = logical(0)
          ),
          v = make_generators(matrix(0, nrow = 1L, ncol = d), NULL, diag(d))
        )
      )
    )
  }
)


method(project, real_region) <- function(space, theta) {
  as.numeric(theta)
}


method(contains, real_region) <- function(space, theta, tol = 1e-8) {
  length(theta) == as.integer(space@d) && all(is.finite(theta))
}


method(region_phrase, real_region) <- function(space) {
  sprintf("all of R^%d", space_dim(space))
}


#' @description Adds a finiteness check to the shared shape checks.
#' @noRd
method(validate_outcome, real_region) <- function(space, x) {
  x <- check_outcome_shape(x, space_dim(space))
  if (any(!is.finite(x))) {
    stop("outcomes must be finite.", call. = FALSE)
  }
  x
}
