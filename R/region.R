#' @include sample_space.R
#' @include polyhedra.R
NULL

# --- The region hierarchy -----------------------------------------------------

#' Regions of a parameter space
#'
#' A `region` is a subset of a family's parameter space: the set a null
#' hypothesis is stated over, or the support of a truncated prior. The `_region`
#' suffix marks the parameter space side of the package throughout, as `_space`
#' marks the sample space side.
#'
#' `region` is abstract and splits in two. A [convex_region] is one that is
#' convex, and carries the geometry: [space_dim()], [contains()], [project()],
#' and a [chart()] to optimise in. A [union_region] is a finite union of convex
#' regions, and need not be convex.
#'
#' Every region answers [parts()] and [cells()]. `parts()` gives the convex
#' regions it was *declared* as; `cells()` gives the convex regions the
#' algorithms *decompose* it into. They agree on every geometry the package
#' currently has, but if one is triangulated they will differ.
#'
#' Part and cell are roles, not types. The same [simplex_region()] is a part
#' when a user declares it as a piece of a plurality null, and a cell when a
#' triangulation produces it from a [polytope_region()]. That is why neither
#' word appears in a class name.
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
region <- new_class("region", abstract = TRUE)


#' The convex regions a region was declared as
#'
#' What the caller asked for, unchanged. A convex region is its own only part.
#'
#' Contrast [cells()], which is what the algorithms decompose a region into.
#' Use `parts()` when reporting what was declared, and `cells()` when feeding an
#' optimiser or an enclosure, or for visualisation.
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
#' @inheritParams parts
#' @return Integer.
#' @examples
#' n_parts(simplex_region(vertices = diag(3)))
#' @export
n_parts <- function(space) length(parts(space))


# Every region is list-like over its parts, so the results of the set algebra
# read uniformly whether one cell survived (a bare convex region), several (a
# union), or none (an empty_region): `length()` is 1, n, or 0, and `x[[i]]`,
# `x[i]`, `lapply()` work on all of them.

#' @rdname parts
#' @usage NULL
#' @export
method(length, region) <- function(x) length(parts(x))


#' @rdname parts
#' @usage NULL
#' @export
method(as.list, region) <- function(x, ...) parts(x)


#' @rdname parts
#' @usage NULL
#' @export
method(`[[`, region) <- function(x, i, ...) parts(x)[[i]]


#' @description Subsetting with `[` returns a region again: the union of the
#'   selected parts, the lone part itself, or an `empty_region` for an empty
#'   selection.
#' @rdname parts
#' @usage NULL
#' @export
method(`[`, region) <- function(x, i, ...) {
  chosen <- parts(x)[i]
  if (length(chosen) == 0L) {
    return(empty_region())
  }
  union_region(chosen)
}


#' The convex regions a region decomposes into for optimisation
#'
#' Every region is its own only cell unless it says otherwise. A region that is
#' a union of convex pieces returns those pieces here, and a geometry that can
#' be triangulated will return its triangulation. Things that search or enclose
#' then run on each convex cell, without needing to know which region it came
#' from.
#'
#' Cells may overlap. This is fine for our case because a supremum over a union
#' is a supremum over any cover. Overlaps may cost time but the results will
#' still be valid.
#'
#' Contrast [parts()], which is what the region was declared as.
#' @param space A [region].
#' @param max_cells The number of simplices triangulation may produce before
#'   giving up (default `1000L`).
#' @param .budget Internal: the shared budget a union hands its parts. Leave
#'   it `NULL`.
#' @return A list of [convex_region] objects whose union is `space`.
#' @examples
#' cells(simplex_region(vertices = diag(3)))
#' @export
cells <- new_generic(
  "cells",
  "space",
  function(space, max_cells = 1000L, .budget = NULL) S7::S7_dispatch()
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


method(cells, region) <- function(space, max_cells = 1000L, .budget = NULL) {
  list(space)
}


#' Number of convex regions a region decomposes into
#' @inheritParams cells
#' @return Integer.
#' @examples
#' n_cells(simplex_region(vertices = diag(3)))
#' @export
n_cells <- function(space) length(cells(space))


# --- Convex regions -----------------------------------------------------------

#' Convex regions
#'
#' A `convex_region` is a convex subset of a family's parameter space. The same
#' type serves two roles. A family's own \eqn{\Theta}{Theta} is one, and so is
#' each convex piece \eqn{\Theta_{0i}}{Theta_0i} of a null hypothesis
#' \eqn{\Theta_0 = \bigcup_i \Theta_{0i}}{Theta_0 = union_i Theta_0i}. A
#' null's pieces may overlap.
#'
#' A `convex_region` object encodes the geometry: it encodes dimension,
#' membership checking, projection, and a [chart()] that maps constrained
#' generator coordinates to the region.
#'
#' Not to be confused with [sample_space]. Outcomes from a sample space are
#' only validated, but in this package parameters need coordinates for a
#' optimiser to search over, which is what [chart()] is for. There may be
#' null geometries that do not permit a [chart()], but these are currently
#' beyond the scope of this package.
#' @examples
#' # `convex_region` is abstract; polytope_region(), simplex_region(),
#' # halfspace_region(), point_region() and unconstrained_region() subclass it:
#' s <- simplex_region(vertices = diag(3))
#' S7::S7_inherits(s, convex_region)
#' space_dim(s)
#' @export
convex_region <- new_class("convex_region", parent = region, abstract = TRUE)


#' Coordinate chart for a parameter space
#'
#' Returns a list of closures that define mappings between the parameter space
#' and a coordinate space, together with the coordinate constraints stated
#' explicitly, which is what lets a constrained optimiser such as SLSQP run
#' on the set.
#' @param space A [convex_region].
#' @return A list comprising:
#' \describe{
#'   \item{`n_par`}{dimension of the coordinate space.}
#'   \item{`to_theta(u)`}{maps coordinates to the parameter vector.}
#'   \item{`to_theta_batch(u_mat)`}{`(n_par, N)` coordinates to `(d, N)`
#'   parameters.}
#'   \item{`from_theta(theta)`}{Coordinates for a point in the parameter space,
#'   satisfying the constraints.}
#'   \item{`jacobian(u)`}{`(d, n_par)` derivative of `to_theta` at `u`.}
#'   \item{`lower`}{length-`n_par` lower bounds on the coordinates; `-Inf`
#'   where a coordinate is free.}
#'   \item{`heq(u)`, `heqjac(u)`}{equality constraint (`heq(u) = 0` on the
#'   feasible set) and its Jacobian, or `NULL` when there is none.}
#'   \item{`seed(n)`}{`(n_par, n)` random feasible coordinates for a
#'   multi-start search, drawn to suit the space's own geometry.}
#' }
#' @examples
#' # One part of the K = 3 plurality null: `theta_1 <= theta_2` in the simplex.
#' s <- simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1)))
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


#' Does a parameter vector belong to the space?
#' @param space A [convex_region].
#' @param theta Parameter vector.
#' @param tol Tolerance.
#' @return `TRUE` or `FALSE`.
#' @examples
#' s <- simplex_region(vertices = diag(3))
#' contains(s, c(1 / 3, 1 / 3, 1 / 3))
#' contains(s, c(2, -1, 0))
#' @export
contains <- new_generic(
  "contains",
  "space",
  function(space, theta, tol = 1e-8) S7::S7_dispatch()
)


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
#' @param space A [convex_region].
#' @param ref Reference parameter vector.
#' @return A parameter vector in the space.
#' @examples
#' s <- simplex_region(vertices = diag(3))
#' init_point(s, c(1, 0, 0))
#' @export
init_point <- new_generic(
  "init_point",
  "space",
  function(space, ref) S7::S7_dispatch()
)


method(init_point, convex_region) <- function(space, ref) project(space, ref)

# --- Representations ----------------------------------------------------------

#' The half-space description of a region
#'
#' Returns the inequalities and equalities cutting a [convex_region] out of its
#' ambient space, with rcdd's flag columns stripped and its sign convention
#' undone, so that `a %*% theta <= b` reads directly.
#'
#' `eq` flags equality to `rcdd`, e.g. for defining the sum(theta)==1 constraint
#' for probabilities.
#'
#' @param space A [convex_region].
#' @return A list with `a` (an `m` by `d` matrix), `b` (length `m`) and `eq`
#'   (logical, length `m`).
#' @keywords internal
#' @noRd
h_rep <- new_generic("h_rep", "space", function(space) S7::S7_dispatch())


#' The generator description of a region
#'
#' Returns the points, rays and lines whose combination is the region, with
#' rcdd's flag columns stripped. Generators are columns, as vertices are
#' everywhere else in the package.
#'
#' `r` and `l` are empty for a bounded region, which is what [is_bounded()]
#' tests. A ray direction is only determined modulo the lineality space, so
#' there can be multiple correct solutions.
#'
#' @param space A [convex_region].
#' @return A list with `v` (vertices), `r` (rays) and `l` (lineality basis),
#'   each a `d` by `n` matrix.
#' @keywords internal
#' @noRd
v_rep <- new_generic("v_rep", "space", function(space) S7::S7_dispatch())


# Both representations are the rational one, converted once at the boundary.
# Every class differs only in how it states itself to cddlib, which is `q_hrep()`
# and `q_vrep()` above.
method(h_rep, convex_region) <- function(space) from_hmatrix(q_hrep(space))


method(v_rep, convex_region) <- function(space) from_vmatrix(q_vrep(space))


#' The rational H- and V-representations of a region
#'
#' The exact form underneath `h_rep()` and `v_rep()`, as an rcdd matrix with its
#' flag columns intact. Composing operations here rather than on the double
#' form is what keeps a chain of them exact.
#'
#' `q_hrep()` of an [unconstrained_region] is the trivially true row
#' `0 . x <= 1` rather than no rows at all, because cddlib needs a row to work
#' with. `h_rep()` reports the same region as zero rows, which is the honest
#' answer for a caller reading facets.
#'
#' @param space A [convex_region].
#' @return An rcdd H- or V-representation matrix of rationals.
#' @keywords internal
#' @noRd
q_hrep <- new_generic("q_hrep", "space", function(space) S7::S7_dispatch())


#' @keywords internal
#' @noRd
q_vrep <- new_generic("q_vrep", "space", function(space) S7::S7_dispatch())


#' Is a region empty?
#'
#' `rcdd` solves a linear feasibility program (a zero objective over the region's
#' half-space description) in exact rational arithmetic. A region is empty if
#' its constraints have no common solution, which will probably happen for cells
#' produced by intersections or complements rather than by declaring an empty
#' set.
#'
#' A [union_region] is empty when every one of its `parts()` is, so this is
#' defined on [region] rather than on [convex_region].
#'
#' @param space A [region].
#' @return `TRUE` or `FALSE`.
#' @examples
#' # The K = 3 plurality cell `{theta_1 <= theta_2}` is not empty:
#' is_empty(simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))))
#'
#' # Nor is a single point, or the whole space:
#' is_empty(point_region(theta = c(0.5, 0.3, 0.2)))
#' is_empty(unconstrained_region(3L))
#' @export
is_empty <- new_generic(
  "is_empty",
  "space",
  function(space) S7::S7_dispatch()
)


method(is_empty, convex_region) <- function(space) {
  q_is_empty(q_hrep(space))
}


#' Is a region bounded?
#'
#' A region is bounded when its generator description has no rays and no
#' lineality.
#'
#' A [union_region] is bounded when every one of its `parts()` is.
#'
#' @param space A [region].
#' @return `TRUE` or `FALSE`.
#' @examples
#' is_bounded(simplex_region(vertices = diag(3)))
#' is_bounded(halfspace_region(normal = c(1, -1, 0)))
#' @export
is_bounded <- new_generic("is_bounded", "space", function(space) {
  S7::S7_dispatch()
})


method(is_bounded, convex_region) <- function(space) {
  v <- v_rep(space)
  ncol(v$r) == 0L && ncol(v$l) == 0L
}


#' The dimension of a region's affine hull
#'
#' An [empty_region()] has dimension -1, a point has dimension 0, a segment 1,
#' and so on.
#'
#' Contrast [space_dim()], which is the number of coordinates a point of the
#' region carries. That ambient dimension belongs to the family's parameter
#' space; the affine dimension belongs to the region's own geometry. It is
#' most the ambient dimension, but may be lower when the region imposes extra
#' constraints (e.g. a lower-dimension null may set theta1 = theta2).
#'
#' @param space A [region].
#' @return An integer between -1 and `space_dim(space)`.
#' @examples
#' region_dim(point_region(theta = c(0.5, 0.3, 0.2)))
#' region_dim(simplex_region(vertices = diag(3)))
#' region_dim(unconstrained_region(3L))
#' @export
region_dim <- new_generic("region_dim", "space", function(space) {
  S7::S7_dispatch()
})


method(region_dim, convex_region) <- function(space) {
  as.integer(q_dim(q_vrep(space)))
}
