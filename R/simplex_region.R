#' @include polytope_region.R
NULL

# --- Simplex region ----------------------------------------------------------

#' A simplex: an affinely independent vertex set
#'
#' The [polytope_region()] whose vertices are affinely independent, so the hull
#' is a simplex of dimension `nrow(vertices) - 1`.
#'
#' [certify()] additionally needs every vertex inside the standard simplex (the
#' tetrahedron below fails this). Such simplices still chart, project and fit;
#' they just cannot be certified.
#'
#' @inheritParams polytope_region
#' @return A `simplex_region`, which is also a [polytope_region()].
#' @examples
#' # The 2-simplex in R^3, e.g. the entire multinomial parameter space:
#' simplex_region(vertices = diag(3))
#'
#' # The plurality region `{theta : theta_1 <= theta_2}` within it:
#' simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1)))
#'
#' # A tetrahedron in R^3, e.g. a piece of a triangulated Gaussian null:
#' simplex_region(
#'   vertices = rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
#' )
#' @export
simplex_region <- new_class(
  "simplex_region",
  parent = polytope_region,
  validator = function(self) simplex_defect(self@hv),
  constructor = function(vertices = NULL, .hv = NULL) {
    if (!is.null(.hv)) {
      return(new_object(polytope_region(.hv = .hv)))
    }
    new_object(polytope_region(vertices = vertices))
  }
)


#' Why a vertex set is not a simplex, or `NULL` if it is
#'
#' Exact: `scdd()` gives the affine hull as equality rows, so hull dimension is
#' `d - sum(eq)`. Also used by `region_from_qh()` to pick a cell's class.
#' @keywords internal
#' @noRd
simplex_defect <- function(hv) {
  n_v <- nrow(hv@v$v)
  d <- ncol(hv@v$v)
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
  hull_dim <- d - sum(hv@h$eq)
  if (n_v > 1L && hull_dim != n_v - 1L) {
    return(paste0(
      "the vertices are affinely dependent: the ",
      n_v,
      " vertices span a hull of dimension ",
      hull_dim,
      " rather than ",
      n_v - 1L
    ))
  }
  NULL
}


method(contains, simplex_region) <- function(space, theta, tol = 1e-8) {
  # Barycentric weights `a` (unique, by affine independence) solve
  # `a V = theta`, `sum(a) = 1`; inside iff all `a >= 0` and the residual is 0.
  v <- space@vertices
  a <- qr.coef(qr(rbind(t(v), 1)), c(theta, 1))
  if (anyNA(a)) {
    return(contains(S7::super(space, to = polyhedron_region), theta, tol))
  }
  all(a >= -tol) && max(abs(drop(a %*% v) - theta)) <= tol
}
