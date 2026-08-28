#' @include polyhedron_region.R
NULL

# --- Halfspace region --------------------------------------------------------

#' Halfspace given by a normal and an offset
#'
#' The unbounded case:
#'   \eqn{\{\theta : a^\top \theta \le b\}}{{theta : a'theta <= b}}
#' (H-representation). The plurality region
#' \eqn{\{\theta : \theta_1 \le \theta_j\}}{{theta : theta_1 <= theta_j}} is
#' `normal = e_1 - e_j`, `offset = 0`.
#'
#' Coordinates are `(z, c)`: `z` positions a point on the bounding hyperplane
#' in an orthonormal basis, and `c >= 0` is the distance inward.
#'
#' Having no vertices, a halfspace admits no certified gap bound. That is a
#' property of the representation rather than a missing feature -- see
#' [polytope_region()] for the bounded alternative.
#'
#' @param normal Normal vector `a`; must be non-zero.
#' @param offset Offset `b`.
#' @return A `halfspace_region`.
#' @examples
#' # `{theta : theta_1 <= theta_2}`
#' halfspace_region(normal = c(1, -1), offset = 0)
#' @export
halfspace_region <- new_class(
  "halfspace_region",
  parent = polyhedron_region,
  properties = list(
    normal = class_numeric,
    offset = class_numeric,
    unit = class_numeric,
    anchor = class_numeric,
    basis = class_any
  ),
  constructor = function(normal, offset = 0) {
    normal <- as.numeric(normal)
    if (all(normal == 0)) {
      stop("`normal` must be a non-zero vector.", call. = FALSE)
    }
    offset <- as.numeric(offset)
    d <- length(normal)
    nrm <- sqrt(sum(normal^2))
    unit <- normal / nrm
    anchor <- offset * normal / nrm^2
    basis <- if (d >= 2L) {
      qr.Q(qr(cbind(unit, diag(d))))[, 2:d, drop = FALSE]
    } else {
      NULL
    }
    new_object(
      polyhedron_region(
        .hv = hv_from_h(
          # Computed from the callers input interpreted as a H-rep.
          # Not interpreting as V-rep because of potential lossy conversions.
          h = list(a = matrix(normal, nrow = 1L), b = offset, eq = FALSE),
          # We override the V-rep double representation because mathematically
          # the H- and V-reps share equivalent numbers, but floating-point
          # magic can perturb.
          v = make_generators(
            vertices = matrix(anchor, ncol = 1L),
            rays = matrix(-unit, ncol = 1L),
            lines = if (d >= 2L) basis else NULL
          )
        )
      ),
      normal = normal,
      offset = offset,
      unit = unit,
      anchor = anchor,
      basis = basis
    )
  }
)


method(project, halfspace_region) <- function(space, theta) {
  slack <- sum(space@normal * theta) - space@offset
  if (slack <= 0) {
    return(theta)
  }
  theta - (slack / sum(space@normal^2)) * space@normal
}


method(contains, halfspace_region) <- function(space, theta, tol = 1e-8) {
  sum(space@normal * theta) <= space@offset + tol * sqrt(sum(space@normal^2))
}


method(region_phrase, halfspace_region) <- function(space) {
  sprintf(
    "the halfspace normal . theta <= %s, normal (%s)",
    format(signif(space@offset, 4L)),
    toString(signif(space@normal, 4L))
  )
}
