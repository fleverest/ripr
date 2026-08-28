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
# union), or none (NULL): `length()` is 1, n, or 0, and `x[[i]]`, `x[i]`,
# `lapply()` work on all of them.

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
#'   selected parts, the lone part itself, or `NULL` for an empty selection.
#' @rdname parts
#' @usage NULL
#' @export
method(`[`, region) <- function(x, i, ...) {
  chosen <- parts(x)[i]
  if (length(chosen) == 0L) {
    return(NULL)
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
#' @param ... Passed on to the method. The triangulating method for
#'   [polyhedron_region()] takes `max_cells` (default `1000L`), the point at
#'   which it gives up.
#' @return A list of [convex_region] objects whose union is `space`.
#' @examples
#' cells(simplex_region(vertices = diag(3)))
#' @export
cells <- new_generic("cells", "space", function(space, ...) S7::S7_dispatch())


method(cells, region) <- function(space, ...) list(space)


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


#' Bundle an objective for [maximise_over()]
#'
#' @param value Function of a parameter vector returning a scalar.
#' @param grad Function of a parameter vector returning the gradient, length `d`.
#' @param value_batch Optional function of a `(d, N)` matrix returning `N`
#'   values. Used to score the seed grid; defaults to applying `value` per
#'   column, which is correct but slower.
#' @return A list to pass to [maximise_over()].
#' @keywords internal
objective <- function(value, grad, value_batch = NULL) {
  if (is.null(value_batch)) {
    value_batch <- function(theta_mat) {
      vapply(seq_len(ncol(theta_mat)), \(i) value(theta_mat[, i]), numeric(1))
    }
  }
  list(value = value, grad = grad, value_batch = value_batch)
}


#' Maximise an objective over a parameter space
#'
#' Multi-start SLSQP in the space's own [chart()]: seed coordinates are scored
#' in batch, the best `n_restarts` are refined under the chart's declared
#' constraints, and the best refinement wins. Written once and shared by all
#' geometries, since only the charts should differ.
#'
#' **The result is a lower bound on the true supremum, not the supremum.** The
#' objective is generally non-convex, so this search is a heuristic. Restarts
#' will converge to local maxima. Anything needing an upper bound must obtain
#' it some other way, e.g. via certification.
#'
#' `seeds` should always include the current mixture's atoms. Without them the
#' returned value can fall below `max_j G(theta_j)`, which the mixture already
#' guarantees is at least 1, and a duality gap computed from it would come out
#' spuriously negative.
#'
#' @param space A [convex_region].
#' @param obj An [objective()].
#' @param seeds Optional `(d, m)` matrix of parameter-space points to seed from.
#' @param n_seeds Random seeds drawn from the chart.
#' @param n_restarts How many of the best seeds to refine.
#' @return `list(theta = , value = )` with `theta` in the space.
#' @keywords internal
maximise_over <- function(
  space,
  obj,
  seeds = NULL,
  n_seeds = 200L,
  n_restarts = 25L
) {
  ch <- chart(space)

  # A zero-dimensional chart has nothing to search: the space is a point.
  if (ch$n_par == 0L) {
    theta <- ch$to_theta(numeric(0))
    return(list(theta = theta, value = obj$value(theta)))
  }

  # slsqp() calls fn and gr separately at the same point, so cache the pair.
  last_u <- NULL
  last <- NULL
  fn_gr <- function(u) {
    if (!is.null(last_u) && identical(u, last_u)) {
      return(last)
    }
    theta <- ch$to_theta(u)
    last_u <<- u
    last <<- list(
      value = -obj$value(theta),
      gradient = -as.numeric(obj$grad(theta) %*% ch$jacobian(u))
    )
    last
  }

  refine <- function(u0, fallback) {
    res <- tryCatch(
      {
        fit <- nloptr::slsqp(
          u0,
          fn = \(u) fn_gr(u)$value,
          gr = \(u) fn_gr(u)$gradient,
          lower = ch$lower,
          heq = ch$heq,
          heqjac = ch$heqjac,
          control = list(xtol_rel = 1e-8, maxeval = 1000L)
        )
        list(par = fit$par, value = fit$value)
      },
      error = function(e) list(par = u0, value = fallback)
    )
    if (!is.finite(res$value)) {
      res <- list(par = u0, value = fallback)
    }
    res
  }

  starts <- ch$seed(n_seeds)
  if (!is.null(seeds)) {
    seeds <- as.matrix(seeds)
    coords <- vapply(
      seq_len(ncol(seeds)),
      \(i) ch$from_theta(seeds[, i]),
      numeric(ch$n_par)
    )
    starts <- cbind(matrix(coords, nrow = ch$n_par), starts)
  }

  scores <- obj$value_batch(ch$to_theta_batch(starts))
  top <- order(scores, decreasing = TRUE)[seq_len(min(
    n_restarts,
    length(scores)
  ))]

  best <- list(par = starts[, top[1L]], value = Inf)
  for (i in top) {
    fallback <- if (is.finite(scores[i])) -scores[i] else Inf
    res <- refine(starts[, i], fallback = fallback)
    if (is.finite(res$value) && res$value < best$value) best <- res
  }
  list(theta = ch$to_theta(best$par), value = -best$value)
}


# --- Shared coordinate helpers ------------------------------------------------

#' Euclidean projection onto the probability simplex
#'
#' Duchi et al. (2008); exact up to floating point, not iterative.
#' @keywords internal
#' @noRd
project_simplex <- function(y) {
  u <- sort(y, decreasing = TRUE)
  css <- cumsum(u)
  rho <- max(which(u + (1 - css) / seq_along(u) > 0))
  pmax(y + (1 - css[rho]) / rho, 0)
}


# --- Polyhedron region --------------------------------------------------------

#' Assemble and check a generator triple
#'
#' Shapes constructor input into the `list(v, r, l)` the class stores: `NULL`
#' blocks become empty, and a cone with no vertex is anchored at the origin by
#' `with_origin_vertex()`. Content validation (finiteness, agreeing dimensions)
#' is the validator's job; this only refuses what it cannot shape.
#' @keywords internal
#' @noRd
make_generators <- function(vertices, rays, lines) {
  given <- Filter(Negate(is.null), list(vertices, rays, lines))
  if (length(given) == 0L) {
    stop(
      "`vertices`, `rays` and `lines` cannot all be NULL.",
      call. = FALSE
    )
  }
  if (!all(vapply(given, is.matrix, logical(1)))) {
    stop(
      "generators must be matrices, one generator per column.",
      call. = FALSE
    )
  }
  d <- nrow(given[[1L]])
  as_block <- function(x) {
    if (is.null(x)) no_generators(d) else x
  }
  with_origin_vertex(
    list(v = as_block(vertices), r = as_block(rays), l = as_block(lines))
  )
}


#' Derive the exact facet matrix at construction, within a size guard
#'
#' Facet count is combinatorial in the generator count in the worst case, and
#' the construction-time measurements stop at 16 vertices in R^5. The
#' constructor leaves `@facets` NULL and `h_rep()` derives on demand instead.
#' Returns the rational matrix directly. The record keeps it, so the exact work
#' is done once rather than recomputed.
#' @keywords internal
#' @noRd
derive_qfacets <- function(qv, n_generators, guard = 100L) {
  if (n_generators > guard) {
    return(NULL)
  }
  q_scdd(qv)
}


#' A convex polyhedron given by its generators
#'
#' The concrete base of every convex region in the package: the Minkowski--Weul
#' form \eqn{\mathrm{conv}(V) + \mathrm{cone}(R) + \mathrm{span}(L)}{conv(V) + cone(R) + span(L)}.
#' [polytope_region()], [simplex_region()], [halfspace_region()],
#' [point_region()] and [unconstrained_region()] are all special cases that
#' have friendlier constructors and do extra validation; use this one when
#' you know the generators themselves, or [h_region()] if you have the
#' H-representation.
#'
#' The [chart()] defines the generator map:
#' \deqn{\theta(u) = V a + L z + R c}{theta(u) = V a + L z + R c}
#' with coordinates `u = (a, z, c)` ordered vertex weights, then lineality,
#' then rays, and the constraints `a >= 0`, `sum(a) = 1` and `c >= 0` declared
#' to the optimiser rather than substituted away. A lone vertex contributes no
#' coordinate, so `n_par = nv + nl + nr` when `nv > 1` and `nl + nr`
#' otherwise. That every polyhedron admits such a generator form, dually to
#' its H-representation, is the Minkowski--Weyl theorem (Ziegler 1995,
#' Theorem 1.2).
#'
#' @param vertices `(d, nv)` numeric matrix, one vertex per column, or `NULL`
#'   for a cone anchored at the origin.
#' @param rays `(d, nr)` numeric matrix of recession directions, or `NULL`.
#' @param lines `(d, nl)` numeric matrix spanning the lineality space, or
#'   `NULL`.
#' @param .hv Internal use only.
#' @return A `polyhedron_region`.
#' @section Properties:
#' \describe{
#'   \item{`generators`}{`list(v, r, l)` of numeric matrices, one generator
#'   per column. Read-only: derived from the static underlying record.}
#'   \item{`facets`}{The half-space description `(A, B, eq)`, read as
#'   `a %*% theta <= b` for each row, with `eq` flagging when equality
#'   should hold rather than inequality for each row. `facets` may be `NULL`
#'   if it has not yet been derived. Read-only.}
#'   \item{`hv`}{The internal dual-representation record the views above
#'   read from, holding both descriptions in both double and rational forms.
#'   This is computed once at construction so that downstream computations
#'   don't start from a rounding.
#' }
#' @references
#'   \insertRef{Ziegler1995}{ripr}
#'
#'   \insertRef{BeckTeboulle2009}{ripr}
#'
#'   \insertRef{DuchiShalevShwartz2008}{ripr}
#'
#'   \insertRef{ODonoghueCandes2015}{ripr}
#' @examples
#' # The halfspace `{theta_1 <= 0}` in R^2, by hand:
#' polyhedron_region(
#'   vertices = matrix(c(0, 0), ncol = 1),
#'   rays = matrix(c(-1, 0), ncol = 1),
#'   lines = matrix(c(0, 1), ncol = 1)
#' )
#' @export
polyhedron_region <- new_class(
  "polyhedron_region",
  parent = convex_region,
  properties = list(
    hv = hv_region,
    generators = new_property(
      class_list,
      getter = function(self) self@hv@v
    ),
    facets = new_property(
      class_any,
      getter = function(self) self@hv@h
    )
  ),
  constructor = function(
    vertices = NULL,
    rays = NULL,
    lines = NULL,
    .hv = NULL
  ) {
    if (!is.null(.hv)) {
      if (!is.null(vertices) || !is.null(rays) || !is.null(lines)) {
        stop(
          "generators and `.hv` cannot both be given: the record already ",
          "carries the generators.",
          call. = FALSE
        )
      }
      return(new_object(S7_object(), hv = .hv))
    }
    # Every region carries its exact rational representations, generated at
    # construction. This way set operations stay exact.
    g <- make_generators(vertices, rays, lines)
    new_object(S7_object(), hv = hv_from_v(g))
  },
  validator = function(self) {
    g <- self@generators
    if (!identical(names(g), c("v", "r", "l"))) {
      return("`generators` must be a list with elements `v`, `r` and `l`")
    }
    if (!all(vapply(g, \(x) is.matrix(x) && is.numeric(x), logical(1)))) {
      return("every generator block must be a numeric matrix")
    }
    dims <- vapply(g, nrow, integer(1))
    if (length(unique(dims)) > 1L) {
      return(paste0(
        "generator blocks disagree on the ambient dimension: ",
        paste(dims, collapse = ", ")
      ))
    }
    if (!all(vapply(g, \(x) all(is.finite(x)), logical(1)))) {
      return("every generator coordinate must be finite")
    }
    if (ncol(g$v) == 0L) {
      return("`generators$v` must hold at least one point")
    }
    f <- self@facets
    if (!is.null(f)) {
      if (!is.list(f) || !all(c("a", "b", "eq") %in% names(f))) {
        return("`facets` must be NULL or a list with `a`, `b` and `eq`")
      }
      if (!is.matrix(f$a) || ncol(f$a) != nrow(g$v)) {
        return("`facets$a` must have one column per ambient dimension")
      }
      if (length(f$b) != nrow(f$a) || length(f$eq) != nrow(f$a)) {
        return("`facets` must have one `b` and one `eq` entry per row of `a`")
      }
    }
    NULL
  }
)


#' A convex polyhedron given by its half-space description
#'
#' The dual constructor to [polyhedron_region()]: the set
#' `{theta : a %*% theta <= b}`, with `eq` flagging rows that hold with
#' equality (e.g. `sum(theta) == 1`). The generators are derived by one exact
#' double-description step, and the rows given here are kept as the region's
#' facets exactly as declared.
#'
#' @param a `(m, d)` numeric matrix of facet normals, one constraint per row.
#' @param b Numeric right-hand side, length `m`.
#' @param eq Logical, length `m` or recycled; `TRUE` marks an equality row.
#' @return A [polyhedron_region()].
#' @examples
#' # The unit square in R^2:
#' h_region(a = rbind(diag(2), -diag(2)), b = c(1, 1, 0, 0))
#'
#' # The halfspace `{theta_1 <= theta_2}`:
#' h_region(a = matrix(c(1, -1), nrow = 1), b = 0)
#' @export
h_region <- function(a, b, eq = FALSE) {
  if (!is.matrix(a) || !is.numeric(a) || nrow(a) == 0L) {
    stop(
      "`a` must be a numeric matrix with one constraint per row.",
      call. = FALSE
    )
  }
  b <- as.numeric(b)
  if (length(b) != nrow(a)) {
    stop("`b` must have one entry per row of `a`.", call. = FALSE)
  }
  eq <- rep_len(as.logical(eq), nrow(a))
  if (anyNA(eq)) {
    stop("`eq` must be TRUE or FALSE for every row.", call. = FALSE)
  }
  h <- list(a = a, b = b, eq = eq)
  if (h_is_empty(h)) {
    stop(
      "the constraints have no common solution: the region would be empty.",
      call. = FALSE
    )
  }
  polyhedron_region(.hv = hv_from_h(h))
}


method(space_dim, polyhedron_region) <- function(space) {
  nrow(space@generators$v)
}


method(v_rep, polyhedron_region) <- function(space) {
  # The generators as declared, not as cddlib would return them: a redundant
  # vertex stays.
  space@generators
}


method(q_vrep, polyhedron_region) <- function(space) {
  space@hv@qv
}


method(h_rep, polyhedron_region) <- function(space) {
  f <- space@facets
  if (!is.null(f)) {
    return(f)
  }
  v_to_h(space@generators)
}


# The record already holds the exact H whenever one was derived; only above
# the facet guard is it absent, and a facet *derived* from a vertex rep in
# general position is an exact rational a double cannot hold, so the on-demand
# path re-runs the exact conversion rather than re-rationalising `@facets`.
method(q_hrep, polyhedron_region) <- function(space) {
  if (!is.null(space@hv@qh)) {
    return(space@hv@qh)
  }
  # Only above the facet guard: the exact H was never derived, so derive on
  # demand without keeping it (S7 value semantics leave nowhere to put it).
  q_scdd(q_vrep(space))
}


method(is_bounded, polyhedron_region) <- function(space) {
  ncol(space@generators$r) == 0L && ncol(space@generators$l) == 0L
}


method(contains, polyhedron_region) <- function(space, theta, tol = 1e-8) {
  # Facet violation, normalised by row norm so the tolerance means the same
  # thing on every facet. Equality rows describe the affine hull and are
  # tested two-sided.
  h <- h_rep(space)
  slack <- as.numeric(h$a %*% theta) - h$b
  scale <- sqrt(rowSums(h$a^2))
  all(ifelse(h$eq, abs(slack) <= tol * scale, slack <= tol * scale))
}


#' Solve for a point's least-squares weights over a generator triple
#'
#' Minimises `|| V a + L z + R c - theta ||^2` over `a` in the simplex,
#' `c >= 0`, `z` free. The minimiser's image `V a + L z + R c` is the
#' Euclidean projection of `theta` onto the polyhedron, so one solve serves
#' both [project()] and the chart's `from_theta()`.
#'
#' Two paths. The direct path substitutes the simplex's affine constraint out,
#' `a = (1 - sum(b), b)`, and solves the remaining least squares by a min-norm
#' SVD solve; the solution is kept whenever it lands in the constraint set,
#' which it does for every point of a region whose vertices are affinely
#' independent -- every simplex cell, in particular -- making the common case
#' exact and non-iterative. Otherwise (a redundant vertex set, or `theta`
#' outside the region) fall through to accelerated projected gradient with
#' per-block proximal steps: simplex projection on `a` (Duchi et al. 2008), a
#' non-negative clamp on `c`, nothing on `z`. FISTA acceleration (Beck and
#' Teboulle 2009) with the momentum reset whenever it points uphill
#' (O'Donoghue and Candes 2015), exactly the scheme of the polytope projection
#' this generalises.
#'
#' @param g `list(v, r, l)`, one generator per column.
#' @param theta Parameter vector.
#' @param tol The tolerance for whether or not we accept the unconstrained sol.
#' @param max_it Maximum number of iterations for proj. grad. descent
#' @return `list(a, z, c, theta_hat)`: the three weight blocks and their image.
#' @keywords internal
#' @noRd
generator_weights <- function(g, theta, tol = 1e-9, max_it = 20000L) {
  # Setup vertex, lineality and ray matrices
  v <- g$v
  l <- g$l
  r <- g$r
  n_v <- ncol(v)
  n_l <- ncol(l)
  n_r <- ncol(r)

  # Maps (a, z, c) |-> V a + L z + R c = theta^*, the projection of theta
  image <- function(a, z, cc) as.numeric(v %*% a + l %*% z + r %*% cc)

  # Now we minimise || theta^* - theta ||^2 subject to the constraints:
  # 1. sum(a) = 1
  # 2. a >= 0     : Because a are barycentric coordinates in V
  # 3. c >= 0    : Because c is a coordinate along a ray, and rays are one-sided

  # First we ignore inequalities and just check if they hold anyway.
  # The constrained problem is harder, so this is easy to check.

  # We force the constraint on a by setting a = (1-sum(a_2,...) a_2 ...)
  # writing b = (a_2 ...), a = (1-sum(b), b). Then V a = v_1 + V^- b, with
  # V^- being all but the first vertex subtracting v_1. Then we solve regular
  # least-squares and check the constraint.
  m <- cbind(v[, -1L, drop = FALSE] - v[, 1L], l, r)
  x <- min_norm_solve(m, theta - v[, 1L])
  b <- x[seq_len(n_v - 1L)]
  a <- c(1 - sum(b), b)
  z <- x[(n_v - 1L) + seq_len(n_l)]
  cc <- x[(n_v - 1L) + n_l + seq_len(n_r)]
  # If the constraints hold (or close enough)
  if (all(a >= -tol) && all(cc >= -tol)) {
    a <- pmax(a, 0)
    a <- a / sum(a)
    cc <- pmax(cc, 0)
    return(list(a = a, z = z, c = cc, theta_hat = image(a, z, cc)))
  }

  # If any weight in a or c is genuinely negative, then we do a
  # constrained solve by projected gradient descent.

  # Setup:
  m <- cbind(v, l, r)
  i_v <- seq_len(n_v)
  i_l <- n_v + seq_len(n_l)
  i_r <- n_v + n_l + seq_len(n_r)

  # Goal: solve argmix_x || M x - theta ||^2 = argmin_x 1/2 || M x - theta ||^2.
  # Gradient is M'(M x - theta), so solve for x:
  # M'M x = M' theta
  mtm <- crossprod(m)
  mtt <- as.numeric(crossprod(m, theta))
  # Lipschitz constant of the gradient of 1/2 || Mx - theta ||^2
  lip <- max(svd(m)$d)^2
  # Project a onto simplex, c onto R+ to satisfy constraints
  prox <- function(x) {
    x[i_v] <- project_simplex(x[i_v])
    x[i_r] <- pmax(x[i_r], 0)
    x
  }

  # Now solve for x starting from default at center of simplex and 0 for c, z.
  x <- c(rep(1 / n_v, n_v), rep(0, n_l), rep(0, n_r))
  y <- x # FISTA extrapolated point
  t_k <- 1 # Momentum counter
  for (i in seq_len(max_it)) {
    new_x <- prox(y - (as.numeric(mtm %*% y) - mtt) / lip)
    if (sum((y - new_x) * (new_x - x)) > 0) {
      y <- new_x
      t_k <- 1
    } else {
      t_new <- (1 + sqrt(1 + 4 * t_k^2)) / 2
      y <- new_x + ((t_k - 1) / t_new) * (new_x - x)
      t_k <- t_new
    }
    converged <- max(abs(new_x - x)) < 1e-14
    x <- new_x
    if (converged) break
  }
  list(
    a = x[i_v],
    z = x[i_l],
    c = x[i_r],
    theta_hat = image(x[i_v], x[i_l], x[i_r])
  )
}


#' Minimum-norm least-squares solution of `m x = rhs`
#'
#' A SVD solve with small singular values dropped.
#' @keywords internal
#' @noRd
min_norm_solve <- function(m, rhs) {
  if (ncol(m) == 0L) {
    return(numeric(0))
  }
  sv <- svd(m)
  keep <- sv$d > max(dim(m)) * .Machine$double.eps * max(sv$d)
  as.numeric(
    sv$v[, keep, drop = FALSE] %*%
      (crossprod(sv$u[, keep, drop = FALSE], rhs) / sv$d[keep])
  )
}


method(project, polyhedron_region) <- function(space, theta) {
  generator_weights(space@generators, theta)$theta_hat
}


method(chart, polyhedron_region) <- function(space) {
  g <- space@generators
  v <- g$v
  l <- g$l
  r <- g$r
  n_v <- ncol(v)
  n_l <- ncol(l)
  n_r <- ncol(r)
  # Direct generator coordinates `u = (a, z, c)`: barycentric weights over the
  # vertices (constrained to the simplex), free lineality coordinates, and
  # non-negative ray coefficients. A lone vertex v_1 contributes no free
  # coordinate: `theta = v_1 + L z + R c`.
  free_v <- if (n_v > 1L) n_v else 0L
  i_v <- seq_len(free_v)
  base <- if (free_v == 0L) as.numeric(v[, 1L]) else numeric(nrow(v))
  jac <- cbind(v[, i_v, drop = FALSE], l, r)
  n_par <- free_v + n_l + n_r

  list(
    n_par = n_par,
    to_theta = function(u) base + as.numeric(jac %*% u),
    to_theta_batch = function(u_mat) jac %*% u_mat + base,
    from_theta = function(theta) {
      w <- generator_weights(g, theta)
      c(if (free_v > 0L) w$a, w$z, w$c)
    },
    jacobian = function(u) jac,
    # Bounds and constraints for a constrained optimiser (SLSQP): the
    # barycentric block sums to one and is non-negative, the ray block is
    # non-negative, the lineality block is free.
    lower = c(rep(0, free_v), rep(-Inf, n_l), rep(0, n_r)),
    heq = if (free_v > 0L) {
      function(u) sum(u[i_v]) - 1
    },
    heqjac = if (free_v > 0L) {
      function(u) matrix(c(rep(1, free_v), rep(0, n_l + n_r)), nrow = 1L)
    },
    # Seeding draws random coefficients for each of a, z and c:
    # a: a uniform Dirichlet over vertices V
    # z: standard normals along the lineality space
    # c: Exp(1) on the rays, boundary-biased with mean 1
    seed = function(n) {
      u_v <- if (free_v > 0L) {
        gm <- matrix(stats::rgamma(n * n_v, shape = 1), nrow = n_v)
        div_by_col(gm, colSums(gm))
      } else {
        matrix(numeric(0), nrow = 0L, ncol = n)
      }
      rbind(
        u_v,
        matrix(stats::rnorm(n * n_l), nrow = n_l, ncol = n),
        matrix(stats::rgamma(n * n_r, shape = 1), nrow = n_r, ncol = n)
      )
    }
  )
}


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

# --- Unconstrained region -----------------------------------------------------

#' The whole of `R^d` as a region
#'
#' The unconstrained case: \eqn{\Theta = \mathbb{R}^d}{Theta = R^d}, with the
#' identity chart. This is the parameter space of a [gaussian_family()], and the
#' way to say that a null places no constraint at all.
#'
#' Being unbounded it has no vertices, so like [halfspace_region()] it admits no
#' certified gap bound.
#'
#' The name says what the region is, rather than what it is made of: a
#' [real_space] is a *sample* space, and the two were too easy to confuse while
#' this one carried the same `real_` prefix.
#'
#' @param d Integer dimension.
#' @return An `unconstrained_region`.
#' @examples
#' unconstrained_region(2L)
#' project(unconstrained_region(2L), c(3, -1))
#' @export
unconstrained_region <- new_class(
  "unconstrained_region",
  parent = polyhedron_region,
  properties = list(
    n_dim = new_property(
      class_numeric,
      getter = function(self) nrow(self@generators$v)
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
        .hv = hv_from_h(
          # A caller reading facets wants to be told there are none.
          # `q_hrep()` states the same region to cddlib as the trivially true
          # row `0 . x <= 1`, via `as_hmatrix()`'s zero-row branch, because
          # cddlib needs a row to work with.
          h = list(
            a = matrix(numeric(0), nrow = 0L, ncol = d),
            b = numeric(0),
            eq = logical(0)
          ),
          v = make_generators(matrix(0, nrow = d, ncol = 1L), NULL, diag(d))
        )
      )
    )
  }
)


method(project, unconstrained_region) <- function(space, theta) {
  as.numeric(theta)
}


method(contains, unconstrained_region) <- function(space, theta, tol = 1e-8) {
  length(theta) == as.integer(space@n_dim) && all(is.finite(theta))
}


# --- Printing -----------------------------------------------------------------

#' A count with its noun, pluralised
#' @keywords internal
#' @noRd
count_label <- function(n, singular, plural = paste0(singular, "s")) {
  sprintf("%d %s", n, if (n == 1L) singular else plural)
}


#' One-line summary of a convex region
#'
#' `format()` appends this to the class name and `print()` puts under the
#' banner. The base counts the generator blocks; classes with a closed form
#' say it directly.
#' @keywords internal
#' @noRd
region_phrase <- new_generic("region_phrase", "space", function(space) {
  S7::S7_dispatch()
})


method(region_phrase, polyhedron_region) <- function(space) {
  g <- space@generators
  counts <- c(
    count_label(ncol(g$v), "vertex", "vertices"),
    if (ncol(g$r) > 0L) count_label(ncol(g$r), "ray"),
    if (ncol(g$l) > 0L) count_label(ncol(g$l), "line")
  )
  sprintf("%s in R^%d", paste(counts, collapse = ", "), space_dim(space))
}


method(region_phrase, point_region) <- function(space) {
  sprintf("the point (%s)", toString(signif(space@theta, 4L)))
}


method(region_phrase, halfspace_region) <- function(space) {
  sprintf(
    "the halfspace normal . theta <= %s, normal (%s)",
    format(signif(space@offset, 4L)),
    toString(signif(space@normal, 4L))
  )
}


method(region_phrase, unconstrained_region) <- function(space) {
  sprintf("all of R^%d", space_dim(space))
}


#' @rdname polyhedron_region
#' @usage NULL
#' @export
method(format, polyhedron_region) <- function(x, ...) {
  sprintf("%s: %s", attr(S7_class(x), "name"), region_phrase(x))
}


#' @description `print()` summarises the geometry -- the generator blocks, the
#'   facet counts, and for a small polytope the vertices themselves -- rather
#'   than dumping the properties.
#' @rdname polyhedron_region
#' @usage NULL
#' @export
method(print, polyhedron_region) <- function(x, ...) {
  cat("<", attr(S7_class(x), "name"), ">\n", sep = "")
  cat(
    "  ",
    region_phrase(x),
    if (!is_bounded(x) && !S7_inherits(x, unconstrained_region)) ", unbounded",
    "\n",
    sep = ""
  )
  f <- x@facets
  if (is.null(f)) {
    cat("  facets not derived: generator count is above the size guard\n")
  } else if (length(f$eq) == 0L) {
    cat("  facets: none\n")
  } else {
    n_eq <- sum(f$eq)
    n_ineq <- length(f$eq) - n_eq
    line <- paste(
      c(
        if (n_ineq > 0L) count_label(n_ineq, "inequality", "inequalities"),
        if (n_eq > 0L) count_label(n_eq, "equality", "equalities")
      ),
      collapse = ", "
    )
    cat("  facets: ", line, "\n", sep = "")
  }
  if (
    S7_inherits(x, polytope_region) &&
      !S7_inherits(x, point_region) &&
      ncol(x@vertices) <= 8L
  ) {
    cat("  vertices:\n")
    print(x@vertices)
  }
  invisible(x)
}
