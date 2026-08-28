#' @include polyhedron_region.R
NULL

# --- Polytope region ----------------------------------------------------------

#' Convex hull of a set of vertices
#'
#' The bounded case: a polytope given by its vertices (a V-representation),
#' parametrised by convex combinations of them.
#'
#' For instance, the plurality region in the standard simplex
#' \eqn{\{\theta : \theta_1 \le \theta_j\}}{{theta : theta_1 <= theta_j}} is of
#' this form.
#'
#' Holding vertices rather than constraints is what a certified bound needs: de
#' Casteljau subdivision works on a vertex set natively. [certify()] takes only
#' the [simplex_region()] cells described below.
#'
#' [simplex_region()] is a special case of [polytope_region()] with affinely
#' independent vertices, which is what the Bernstein enclosure in [certify()]
#' needs upstream.
#'
#' Use [simplex_region()] when that is what you mean, and use [polytope_region()]
#' when it is not.
#'
#' [cells()] bridges the two. A polytope that is not a simplex decomposes into
#' simplices by a vertex fan triangulation. A square gives two triangles; a
#' pentagon, three.
#' @param vertices `(d, V)` numeric matrix, one vertex per column.
#' @param .hv Internal only: the underlying dual-representation record in both
#'   double and rational representation.
#' @return A `polytope_region`, which is also a [polyhedron_region()] with
#'   empty ray and lineality blocks.
#' @examples
#' # A square in R^2
#' polytope_region(vertices = cbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)))
#' @export
polytope_region <- new_class(
  "polytope_region",
  parent = polyhedron_region,
  properties = list(
    vertices = new_property(
      class_any,
      getter = function(self) self@generators$v
    ),
    n_vertices = new_property(
      class_numeric,
      getter = function(self) ncol(self@generators$v)
    )
  ),
  constructor = function(vertices = NULL, .hv = NULL) {
    if (!is.null(.hv)) {
      return(new_object(polyhedron_region(.hv = .hv)))
    }
    if (!is.matrix(vertices) || ncol(vertices) == 0L) {
      stop(
        "`vertices` must be a matrix with one vertex per column.",
        call. = FALSE
      )
    }
    if (!all(is.finite(vertices))) {
      stop("`vertices` must all be finite.", call. = FALSE)
    }
    new_object(polyhedron_region(vertices = vertices))
  },
  validator = function(self) {
    if (ncol(self@generators$r) > 0L || ncol(self@generators$l) > 0L) {
      return("a polytope is bounded, so rays and lines must be empty")
    }
    NULL
  }
)
