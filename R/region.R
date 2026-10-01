#' @include space.R
#' @include polyhedra.R
NULL

# --- The region hierarchy -----------------------------------------------------

#' Regions of a parameter space
#'
#' A `region` is a subset of \eqn{\mathbb{R}^d}{R^d} with convex geometry: a
#' null hypothesis, a family's parameter space, or a truncated prior's support.
#' Unlike a [count_space], it can be charted, projected onto and combined with
#' the set algebra.
#'
#' `region` is abstract. A [convex_region] carries the geometry
#' ([space_dim()], [contains()], [project()], [chart()]); a [union_region] is a
#' finite union of convex regions and need not be convex. Every region answers
#' [parts()] (the convex regions it was declared as) and [cells()] (the convex
#' regions the algorithms decompose it into); these differ once a region is
#' triangulated.
#' @examples
#' # Every convex_region is a region:
#' s <- simplex_region(vertices = diag(3))
#' S7::S7_inherits(s, region)
#' S7::S7_inherits(s, convex_region)
#'
#' # A union of them is a region, but not a convex one:
#' u <- union_region(s, halfspace_region(normal = c(1, -1, 0)))
#' S7::S7_inherits(u, region)
#' S7::S7_inherits(u, convex_region)
#' @export
region <- new_class("region", parent = space, abstract = TRUE)


#' The convex regions a region was declared as
#'
#' What was originally constructed, unchanged; a convex region is its own only
#' part. Use [cells()] instead when feeding an optimiser or for visualisation.
#' @param space A [region].
#' @return A list of [convex_region] objects.
#' @examples
#' s <- simplex_region(vertices = diag(3))
#' parts(s)
#' parts(union_region(s, halfspace_region(normal = c(1, -1, 0))))
#' @export
parts <- new_generic("parts", "space", function(space) S7::S7_dispatch())


method(parts, region) <- function(space) list(space)


#' Number of convex regions a region was declared as
#' @keywords internal
#' @noRd
n_parts <- function(space) length(parts(space))


#' The convex regions a region decomposes into for optimisation
#'
#' A union returns its pieces and a triangulable geometry its triangulation;
#' any other region is its own only cell. Cells may overlap.
#' Contrast [parts()], which is what the region was declared as.
#' @param space A [region].
#' @param max_cells The number of simplices triangulation may produce before
#'   giving up. Defaults to the `ripr.max_cells` option, or `1000L` when that
#'   is unset.
#' @param .budget Internal: the shared budget a union hands its parts. Leave
#'   it `NULL`.
#' @return A list of [convex_region] objects whose union is `space`.
#' @examples
#' cells(simplex_region(vertices = diag(3)))
#' @export
cells <- new_generic(
  "cells",
  "space",
  function(
    space,
    max_cells = getOption("ripr.max_cells", 1000L),
    .budget = NULL
  ) {
    S7::S7_dispatch()
  }
)


#' The simplex budget one `cells()` call runs under
#' @keywords internal
#' @noRd
cell_budget <- function(max_cells) {
  budget <- new.env(parent = emptyenv())
  budget$left <- max_cells
  budget$max_cells <- max_cells
  budget
}


method(cells, region) <- function(
  space,
  max_cells = getOption("ripr.max_cells", 1000L),
  .budget = NULL
) {
  list(space)
}


# --- Convex regions -----------------------------------------------------------

#' Convex regions
#'
#' A convex subset of a parameter space: a family's own \eqn{\Theta}{Theta},
#' or one (possibly overlapping) piece \eqn{\Theta_{0i}}{Theta_0i} of a null
#' \eqn{\Theta_0 = \bigcup_i \Theta_{0i}}{Theta_0 = union_i Theta_0i}. It
#' provides dimension, membership, projection and a [chart()] of constrained
#' coordinates for optimisers.
#' @examples
#' # `convex_region` is abstract; polytope_region(), simplex_region(),
#' # halfspace_region(), point_region() and real_region() subclass it:
#' s <- simplex_region(vertices = diag(3))
#' S7::S7_inherits(s, convex_region)
#' space_dim(s)
#' @export
convex_region <- new_class("convex_region", parent = region, abstract = TRUE)


#' Coordinate chart for a parameter space
#'
#' Closures mapping between the parameter space and a coordinate space, with
#' the coordinate constraints stated so a constrained optimiser (SLSQP) can run
#' on the set.
#' @param space A [convex_region].
#' @return A list comprising:
#' \describe{
#'   \item{`n_par`}{dimension of the coordinate space.}
#'   \item{`to_theta(u)`}{maps coordinates to the parameter vector.}
#'   \item{`to_theta_batch(u_mat)`}{`(N, n_par)` coordinates, one point per
#'   row, to `(N, d)` parameters.}
#'   \item{`from_theta(theta)`}{Coordinates for a point in the parameter space,
#'   satisfying the constraints.}
#'   \item{`jacobian(u)`}{`(d, n_par)` derivative of `to_theta` at `u`.}
#'   \item{`lower`}{length-`n_par` lower bounds on the coordinates; `-Inf`
#'   where a coordinate is free.}
#'   \item{`heq(u)`, `heqjac(u)`}{equality constraint (`heq(u) = 0` on the
#'   feasible set) and its Jacobian, or `NULL` when there is none.}
#'   \item{`seed(n)`}{`(n, n_par)` random feasible coordinates, one point per
#'   row, for a multi-start search, drawn to suit the space's own geometry.}
#' }
#' @examples
#' # One part of the K = 3 plurality null: `theta_1 <= theta_2` in the simplex.
#' s <- simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1)))
#' ch <- chart(s)
#' ch$n_par
#'
#' # Coordinates are barycentric weights over the three vertices ...
#' ch$to_theta(c(1 / 2, 1 / 4, 1 / 4))
#' ch$from_theta(c(0.1, 0.6, 0.3))
#'
#' # ... constrained to the simplex.
#' ch$lower
#' ch$heq(c(1 / 2, 1 / 4, 1 / 4))
#'
#' # An unbounded region also has lineality and cone coordinates. The halfspace
#' # `{theta_1 <= theta_2}` has one on its bounding hyperplane and one for the
#' # distance inward, the latter bounded below:
#' ch <- chart(halfspace_region(normal = c(1, -1), offset = 0))
#' ch$lower
#' ch$to_theta(c(0, sqrt(2)))
#' @export
chart <- new_generic("chart", "space", function(space) S7::S7_dispatch())


#' Euclidean projection onto a parameter space
#'
#' The closest point of the space to `theta`. Idempotent up to tolerance, and
#' its output always satisfies [contains()].
#' @param space A [convex_region].
#' @param theta Parameter vector.
#' @return A parameter vector in the space.
#' @examples
#' s <- simplex_region(vertices = diag(3))
#' project(s, c(2, -1, 0))
#' @export
project <- new_generic(
  "project",
  "space",
  function(space, theta) S7::S7_dispatch()
)


#' A starting atom in a parameter space
#'
#' Defaults to projecting a reference point, e.g. the alternative's mean.
#' @keywords internal
#' @noRd
init_point <- new_generic(
  "init_point",
  "space",
  function(space, ref) S7::S7_dispatch()
)


method(init_point, convex_region) <- function(space, ref) project(space, ref)


# --- Predicates ---------------------------------------------------------------

#' Is a region empty?
#'
#' Decided by an exact rational LP feasibility check (`rcdd`). A [union_region]
#' is empty when all its `parts()` are.
#' @param space A [region].
#' @return `TRUE` or `FALSE`.
#' @examples
#' # The K = 3 plurality cell `{theta_1 <= theta_2}` is not empty:
#' is_empty(simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))))
#'
#' # Nor is a single point, or the whole space:
#' is_empty(point_region(theta = c(0.5, 0.3, 0.2)))
#' is_empty(real_region(3L))
#' @export
is_empty <- new_generic(
  "is_empty",
  "space",
  function(space) S7::S7_dispatch()
)


#' Is a region bounded?
#'
#' Bounded means no rays and no lineality; a [union_region] is bounded when
#' all its `parts()` are.
#' @param space A [region].
#' @return `TRUE` or `FALSE`.
#' @examples
#' is_bounded(simplex_region(vertices = diag(3)))
#' is_bounded(halfspace_region(normal = c(1, -1, 0)))
#' @export
is_bounded <- new_generic("is_bounded", "space", function(space) {
  S7::S7_dispatch()
})
