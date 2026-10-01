#' @include simplex_region.R
NULL

# --- Point region --------------------------------------------------------

#' A single parameter point
#'
#' The degenerate convex set \eqn{\{\theta\}}{{theta}}, e.g. the null for
#' which a likelihood ratio \eqn{Q / P_\theta}{Q / P_theta} is an e-variable.
#' It is the 0-simplex; [certify()] handles it by direct evaluation.
#'
#' @param theta The parameter vector.
#' @return A `point_region`, which is also a [simplex_region()].
#' @examples
#' point_region(theta = c(0.5, 0.3, 0.2))
#' @export
point_region <- new_class(
  "point_region",
  parent = simplex_region,
  properties = list(
    theta = new_property(
      class_numeric,
      getter = function(self) as.numeric(self@generators$v[1L, ])
    )
  ),
  constructor = function(theta) {
    theta <- as.numeric(theta)
    if (length(theta) == 0L || !all(is.finite(theta))) {
      stop("`theta` must be a finite numeric vector.", call. = FALSE)
    }
    d <- length(theta)
    new_object(simplex_region(
      .hv = hv_fill(
        h = list(a = diag(d), b = theta, eq = rep(TRUE, d)),
        v = make_generators(matrix(theta, nrow = 1L), NULL, NULL)
      )
    ))
  },
  validator = function(self) {
    if (nrow(self@generators$v) != 1L) {
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
