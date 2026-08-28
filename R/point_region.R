#' @include polytope_region.R
NULL

# --- Point region --------------------------------------------------------

#' A single parameter point
#'
#' The degenerate convex set \eqn{\{\theta\}}{{theta}}. Its role is to let a
#' likelihood ratio state the null it is valid for: \eqn{Q / P_\theta}{Q / P_theta}
#' is an e-variable for \eqn{\{P_\theta\}}{{P_theta}}.
#'
#' It is technically a simplex, but it requires no validation so we have a
#' separate class for it. [certify()] treats it separately too: the supremum
#' over one parameter is the expectation at that point, so a point null is
#' certified by direct evaluation rather than by any enclosing method.
#'
#' @param theta The parameter vector.
#' @return A `point_region`.
#' @examples
#' point_region(theta = c(0.5, 0.3, 0.2))
#' @export
point_region <- new_class(
  "point_region",
  parent = polytope_region,
  properties = list(
    theta = new_property(
      class_numeric,
      getter = function(self) as.numeric(self@generators$v[, 1L])
    )
  ),
  constructor = function(theta) {
    theta <- as.numeric(theta)
    if (length(theta) == 0L || !all(is.finite(theta))) {
      stop("`theta` must be a finite numeric vector.", call. = FALSE)
    }
    d <- length(theta)
    # The exact facets `theta_i == b_i`, rather than letting the polytope
    # constructor derive rows that describe the same point less directly.
    new_object(polytope_region(
      .hv = hv_from_h(
        h = list(a = diag(d), b = theta, eq = rep(TRUE, d)),
        v = make_generators(matrix(theta, ncol = 1L), NULL, NULL)
      )
    ))
  },
  validator = function(self) {
    if (ncol(self@generators$v) != 1L) {
      return("a point region holds exactly one vertex")
    }
    NULL
  }
)


method(project, point_region) <- function(space, theta) space@theta

method(contains, point_region) <- function(space, theta, tol = 1e-8) {
  max(abs(space@theta - theta)) <= tol
}


method(region_phrase, point_region) <- function(space) {
  sprintf("the point (%s)", toString(signif(space@theta, 4L)))
}
