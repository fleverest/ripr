#' @include polyhedron_region.R
NULL

# --- Polytope region ----------------------------------------------------------

#' Convex hull of a set of vertices
#'
#' A bounded polytope given by its vertices, parametrised by convex
#' combinations of them; e.g. the plurality region
#' \eqn{\{\theta : \theta_1 \le \theta_j\}}{{theta : theta_1 <= theta_j}}
#' in the standard simplex. Use [simplex_region()] when the vertices are
#' affinely independent. [cells()] triangulates a polytope into simplices,
#' which is what [certify()] works on.
#' @param vertices `(V, d)` numeric matrix, one vertex per row.
#' @param .hv Internal only: the underlying dual-representation record in both
#'   double and rational representation.
#' @return A `polytope_region`, which is also a [polyhedron_region()] with
#'   empty ray and lineality blocks.
#' @examples
#' # A square in R^2
#' polytope_region(vertices = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)))
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
      getter = function(self) nrow(self@generators$v)
    )
  ),
  constructor = function(vertices = NULL, .hv = NULL) {
    if (!is.null(.hv)) {
      return(new_object(polyhedron_region(.hv = .hv)))
    }
    if (!is.matrix(vertices) || nrow(vertices) == 0L) {
      stop(
        "`vertices` must be a matrix with one vertex per row.",
        call. = FALSE
      )
    }
    if (!all(is.finite(vertices))) {
      stop("`vertices` must all be finite.", call. = FALSE)
    }
    new_object(polyhedron_region(vertices = vertices))
  },
  validator = function(self) {
    if (nrow(self@generators$r) > 0L || nrow(self@generators$l) > 0L) {
      return("a polytope is bounded, so rays and lines must be empty")
    }
    NULL
  }
)
