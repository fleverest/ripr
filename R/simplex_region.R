#' @include polytope_region.R
NULL

# --- Simplex region ----------------------------------------------------------

#' A simplex: an affinely independent vertex set
#'
#' The [polytope_region()] whose vertices are affinely independent, so the hull
#' is a simplex of dimension `ncol(vertices) - 1`.
#'
#' Being a simplex is not on its own enough to certify. The Bernstein enclosure
#' asks for two further things, both checked by [certify()] rather than here:
#' one vertex per coordinate, and every vertex inside the standard simplex. A
#' lower-dimensional cell fails the first; the tetrahedron in the examples below
#' fails the second, being a full-dimensional simplex of `R^3` rather than of
#' \eqn{\Delta}{Delta}. Either way [certify()] refuses, naming the condition.
#'
#' Neither is a defect in the region. Both still chart, project and fit as
#' usual, and `sup_lb()` still searches them; they simply have no implemented
#' bounding method, as a [gaussian_family()] null already does not.
#'
#' @inheritParams polytope_region
#' @return A `simplex_region`, which is also a [polytope_region()].
#' @examples
#' # The 2-simplex in R^3, e.g. the entire multinomial parameter space:
#' simplex_region(vertices = diag(3))
#'
#' # The plurality region `{theta : theta_1 <= theta_2}` within it:
#' simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1)))
#'
#' # A tetrahedron in R^3, e.g. a piece of a triangulated Gaussian null:
#' simplex_region(
#'   vertices = cbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
#' )
#' @export
simplex_region <- new_class(
  "simplex_region",
  parent = polytope_region,
  validator = function(self) {
    vertices <- self@vertices
    if (!all(is.finite(vertices))) {
      return("every vertex coordinate must be finite")
    }
    n_v <- ncol(vertices)
    d <- nrow(vertices)
    if (n_v > d + 1L) {
      return(paste0(
        "at most ",
        d + 1L,
        " points can be affinely independent in ",
        d,
        " dimensions; got ",
        n_v,
        ". A hull of more vertices is a `polytope_region`"
      ))
    }
    if (n_v > 1L) {
      # Affine independence
      eq_rows <- if (!is.null(self@facets)) {
        sum(self@facets$eq)
      } else {
        sum(q_hrep(self)[, 1L] == "1")
      }
      hull_dim <- d - eq_rows
      if (hull_dim != n_v - 1L) {
        return(paste0(
          "the vertices are affinely dependent: the ",
          n_v,
          " vertices span a hull of dimension ",
          hull_dim,
          " rather than ",
          n_v - 1L
        ))
      }
    }
    NULL
  },
  constructor = function(vertices = NULL, .hv = NULL) {
    if (!is.null(.hv)) {
      return(new_object(polytope_region(.hv = .hv)))
    }
    new_object(polytope_region(vertices = vertices))
  }
)


method(contains, simplex_region) <- function(space, theta, tol = 1e-8) {
  # Barycentric membership test.
  #
  # The vertices are affinely independent so every point of the hull has unique
  # weights. A point being inside means each solved weight is non-negative and
  # the least-squares residual is zero (up to some small tolerance).
  v <- space@vertices
  a <- qr.coef(qr(rbind(v, 1)), c(theta, 1))
  if (anyNA(a)) {
    return(contains(S7::super(space, to = polyhedron_region), theta, tol))
  }
  all(a >= -tol) && max(abs(v %*% a - theta)) <= tol
}
