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
#' @return A list of [convex_region] objects whose union is `space`.
#' @examples
#' cells(simplex_region(vertices = diag(3)))
#' @export
cells <- new_generic("cells", "space", function(space) S7::S7_dispatch())


method(cells, region) <- function(space) list(space)


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
#' membership checking, projection, and a [chart()] that maps unconstrained
#' coordinates to the region.
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


#' Unconstrained coordinate chart for a parameter space
#'
#' Returns a list of closures that define mappings between the parameter space
#' and an unconstrained coordinate space, which is what lets BFGS run on the
#' constrained set.
#'
#' Charts for a compact space generally cover only the relative interior, so
#' an optimiser never exactly solves a maximum attained at a vertex or at
#' infinity, though at this point we are in the realm of numerical precision
#' anyway.
#' @param space A [convex_region].
#' @return A list of closures comprising:
#' \describe{
#'   \item{`n_par`}{dimension of the coordinate space.}
#'   \item{`to_theta(u)`}{coordinates to a parameter vector in the space.}
#'   \item{`to_theta_batch(u_mat)`}{`(n_par, N)` coordinates to `(d, N)`
#'   parameters.}
#'   \item{`from_theta(theta)`}{Reparemetrised coordinates for a point in the
#'   space.}
#'   \item{`jacobian(u)`}{`(d, n_par)` derivative of `to_theta` at `u`.}
#'   \item{`seed(n)`}{`(n_par, n)` random coordinates for a multi-start search,
#'   drawn to suit the space's own geometry.}
#' }
#' @examples
#' s <- simplex_region(vertices = diag(3))
#' ch <- chart(s)
#' ch$n_par
#' ch$to_theta(c(0, 0))
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
#' Multi-start BFGS in the space's own [chart()]: seed coordinates are scored in
#' batch, the best `n_restarts` are refined, and the best refinement wins.
#' Written once and shared by every geometry, since only the chart differs.
#'
#' **The result is a lower bound on the true supremum, not the supremum.** The
#' objective is generally non-convex, the search is heuristic, and charts cover
#' only the relative interior, so a maximum at a vertex is approached and never
#' attained. Anything needing an upper bound must obtain it elsewhere.
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

  # optim() calls fn and gr separately at the same point, so cache the pair.
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
    tryCatch(
      stats::optim(
        u0,
        fn = \(u) fn_gr(u)$value,
        gr = \(u) fn_gr(u)$gradient,
        method = "BFGS"
      ),
      error = function(e) list(par = u0, value = fallback)
    )
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
    res <- refine(starts[, i], fallback = -scores[i])
    if (res$value < best$value) best <- res
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


#' Softmax of `c(0, u)`, a bijection from `R^{n-1}` to the simplex interior
#' @keywords internal
#' @noRd
softmax0 <- function(u) {
  e <- exp(c(0, u) - max(0, u))
  e / sum(e)
}


#' Inverse of softmax0, with a guard against exact zeros
#' @keywords internal
#' @noRd
softmax0_inv <- function(alpha, eps = 1e-12) {
  log(alpha[-1L] + eps) - log(alpha[1L] + eps)
}


#' Jacobian `d alpha / d u` of softmax0, `(n, n-1)`
#' @keywords internal
#' @noRd
softmax0_jacobian <- function(alpha) {
  (diag(alpha) - outer(alpha, alpha))[, -1L, drop = FALSE]
}


#' @keywords internal
#' @noRd
softplus <- function(s) log1p(exp(-abs(s))) + pmax(s, 0)


#' @keywords internal
#' @noRd
softplus_inv <- function(t) log(expm1(pmax(t, 1e-12)))


#' @keywords internal
#' @noRd
sigmoid <- function(s) 1 / (1 + exp(-s))


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


#' Derive a facet description at construction, within a size guard
#'
#' Facet count is combinatorial in the generator count in the worst case, and
#' the construction-time measurements stop at 16 vertices in R^5. Above
#' `guard` generators -- columns: one generator per column, as everywhere in
#' the package -- the constructor leaves `@facets` NULL and `h_rep()` derives
#' on demand instead.
#' @keywords internal
#' @noRd
derive_facets <- function(g, guard = 100L) {
  if (ncol(g$v) + ncol(g$r) + ncol(g$l) > guard) {
    return(NULL)
  }
  v_to_h(g)
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
#' The [chart()] is based directly on the generators:
#' \deqn{\theta(u) = V\,\mathrm{softmax0}(u_v) + L u_l + R\,\mathrm{softplus}(u_r)}{
#' theta(u) = V softmax0(u_v) + L u_l + R softplus(u_r)}
#' with `n_par = (nv - 1) + nl + nr` coordinates ordered `(u_v, u_l, u_r)`:
#' vertex weights, then lineality, then rays. That every polyhedron admits
#' such a generator form, dually to its H-representation, is the Minkowski--Weyl
#' theorem (Ziegler 1995, Theorem 1.2).
#'
#' @param vertices `(d, nv)` numeric matrix, one vertex per column, or `NULL`
#'   for a cone anchored at the origin.
#' @param rays `(d, nr)` numeric matrix of recession directions, or `NULL`.
#' @param lines `(d, nl)` numeric matrix spanning the lineality space, or
#'   `NULL`.
#' @param facets Optionally, the half-space description already known to the
#'   caller: a list with `a` (`(m, d)` matrix), `b` (length `m`) and `eq`
#'   (logical, length `m`), read as `a %*% theta <= b` with `eq` flagging
#'   equality rows.
#' @return A `polyhedron_region`.
#' @section Properties:
#' \describe{
#'   \item{`generators`}{`list(v, r, l)` of numeric matrices, one generator
#'   per column.}
#'   \item{`facets`}{The half-space description as above, or `NULL` when it has
#'   not been derived.}
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
    generators = class_list,
    facets = class_any
  ),
  constructor = function(
    vertices = NULL,
    rays = NULL,
    lines = NULL,
    facets = NULL
  ) {
    g <- make_generators(vertices, rays, lines)
    if (is.null(facets)) {
      facets <- derive_facets(g)
    }
    new_object(S7_object(), generators = g, facets = facets)
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
  g <- h_to_v(h)
  polyhedron_region(
    vertices = g$v,
    rays = g$r,
    lines = g$l,
    facets = h
  )
}


method(space_dim, polyhedron_region) <- function(space) {
  nrow(space@generators$v)
}


method(v_rep, polyhedron_region) <- function(space) {
  # The generators as declared, not as cddlib would return them: a redundant
  # vertex stays.
  space@generators
}


method(q_vrep, polyhedron_region) <- function(space) as_vmatrix(v_rep(space))


method(h_rep, polyhedron_region) <- function(space) {
  f <- space@facets
  if (!is.null(f)) {
    return(f)
  }
  v_to_h(space@generators)
}


# A facet *derived* from vertex rep in general position is an exact rational a
# double cannot hold, so the rational H-representation always re-runs the
# exact conversion rather than re-rationalising the stored double `@facets`.
# Can override this with `as_hmatrix(h_rep(space))` when the facets are,
# declared directly by the user, for them it is exact and cheaper.
method(q_hrep, polyhedron_region) <- function(space) q_scdd(q_vrep(space))


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
  i_v <- seq_len(n_v - 1L)
  i_l <- (n_v - 1L) + seq_len(n_l)
  i_r <- (n_v - 1L) + n_l + seq_len(n_r)

  to_theta <- function(u) {
    as.numeric(v %*% softmax0(u[i_v]) + l %*% u[i_l] + r %*% softplus(u[i_r]))
  }

  list(
    n_par = (n_v - 1L) + n_l + n_r,
    to_theta = to_theta,
    to_theta_batch = function(u_mat) {
      alpha <- vapply(
        seq_len(ncol(u_mat)),
        \(i) softmax0(u_mat[i_v, i]),
        numeric(n_v)
      )
      out <- v %*% matrix(alpha, nrow = n_v)
      if (n_l > 0L) {
        out <- out + l %*% u_mat[i_l, , drop = FALSE]
      }
      if (n_r > 0L) {
        out <- out + r %*% softplus(u_mat[i_r, , drop = FALSE])
      }
      out
    },
    from_theta = function(theta) {
      w <- generator_weights(g, theta)
      c(softmax0_inv(w$a), w$z, softplus_inv(w$c))
    },
    jacobian = function(u) {
      cbind(
        v %*% softmax0_jacobian(softmax0(u[i_v])),
        l,
        r %*% diag(sigmoid(u[i_r]), nrow = n_r)
      )
    },
    # Seeding draws random coefficients for each of a, z and c:
    # a: a uniform Dirichlet over vertices V
    # z: standard normals along the lineality space
    # c: boundary-biased normals on the rays
    # It is really just a heuristic for the search space.
    # The entire chart is heuristic, softmax/softplus parametrisation
    # may not be ideal, but it permits unconstrained optimisers like BFGS to run
    # without much thought.
    seed = function(n) {
      u_v <- if (n_v > 1L) {
        # Dirichlet draws
        gm <- matrix(stats::rgamma(n * n_v, shape = 1), nrow = n_v)
        alpha <- div_by_col(gm, colSums(gm))
        matrix(
          vapply(seq_len(n), \(i) softmax0_inv(alpha[, i]), numeric(n_v - 1L)),
          nrow = n_v - 1L
        )
      } else {
        matrix(numeric(0), nrow = 0L, ncol = n)
      }
      rbind(
        # Dirichlet for V coeffs
        u_v,
        # Std normal for lineality space
        matrix(stats::rnorm(n * n_l), nrow = n_l, ncol = n),
        # normal with mean -1, stdev 2 for rays
        # (slightly biased toward 0 through softplus)
        matrix(stats::rnorm(n * n_r, mean = -1, sd = 2), nrow = n_r, ncol = n)
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
#' @param vertices `(d, V)` numeric matrix, one vertex per column.
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
  constructor = function(vertices) {
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
      eq_rows <- if (is.null(self@facets)) {
        sum(q_hrep(self)[, 1L] == "1")
      } else {
        sum(self@facets$eq)
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
  constructor = function(vertices) {
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
#' Coordinates are `(z, s)`: `z` positions a point on the bounding hyperplane in
#' an orthonormal basis, and `softplus(s) >= 0` is the distance inward.
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
        # The foot of the normal on the bounding hyperplane, the outward
        # direction as the one ray, and an orthonormal basis of the hyperplane,
        # which is exactly the lineality space.
        vertices = matrix(anchor, ncol = 1L),
        rays = matrix(-unit, ncol = 1L),
        lines = if (d >= 2L) basis else NULL,
        # The caller's own inequality, exactly: deriving it back from the
        # QR-computed basis through cddlib would be strictly less exact.
        facets = list(a = matrix(normal, nrow = 1L), b = offset, eq = FALSE)
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

# `@facets` holds the declared inequality exactly, so rationalising it via
# as_hmatrix is exact.
method(q_hrep, halfspace_region) <- function(space) as_hmatrix(h_rep(space))


# --- Point region --------------------------------------------------------

#' A single parameter point
#'
#' The degenerate convex set \eqn{\{\theta\}}{{theta}}. Its role is to let a
#' likelihood ratio state the null it is valid for: \eqn{Q / P_\theta}{Q / P_theta}
#' is an e-variable for \eqn{\{P_\theta\}}{{P_theta}}.
#'
#' It is technically a simplex, but it requires no validation so we have a
#' separate class for it.
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
    # The exact facets `theta_i == b_i` replace the ones the polytope
    # constructor derives through cddlib, which describe the same point less
    # directly.
    new_object(
      polytope_region(vertices = matrix(theta, ncol = 1L)),
      facets = list(a = diag(d), b = theta, eq = rep(TRUE, d))
    )
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

# `@facets` holds `theta` itself as equality rows, exact under rationalisation.
method(q_hrep, point_region) <- function(space) as_hmatrix(h_rep(space))


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
        vertices = matrix(0, nrow = d, ncol = 1L),
        lines = diag(d),
        # A caller reading facets wants to be told there are none. `q_hrep()`
        # states the same region to cddlib as the trivially true row
        # `0 . x <= 1`, via `as_hmatrix()`'s zero-row branch, because cddlib
        # needs a row to work with.
        facets = list(
          a = matrix(numeric(0), nrow = 0L, ncol = d),
          b = numeric(0),
          eq = logical(0)
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

# `@facets` is the exact zero-row description, so no re-derivation is needed.
method(q_hrep, unconstrained_region) <- function(space) as_hmatrix(h_rep(space))


# --- Union region -------------------------------------------------------------

#' A finite union of convex regions
#'
#' The union \eqn{\bigcup_i \Theta_{0i}}{union_i Theta_0i} of finitely many
#' [convex_region]s, which generally is not convex. A null hypothesis is one
#' such union, and so is the support of a truncated prior, so the union is
#' worth a type of its own rather than an untyped list passed around by
#' whoever happens to hold it.
#'
#' A `union_region` is a [region] but deliberately **not** a [convex_region].
#' [chart()], [project()], `maximise_over()` assume convexity, and a union of
#' convex sets is not convex. What this class does implement is [space_dim()]
#' [contains()], [parts()] and [cells()].
#'
#' Given exactly one convex region, `union_region()` returns it unchanged.
#'
#' @param ... [convex_region] objects, other `union_region` objects, and lists
#'   of either, in any combination and any nesting. A `union_region` argument
#'   flattens rather than nests.
#' @return A `union_region`, or the lone [convex_region] it was given.
#' @section Properties:
#' \describe{
#'   \item{`parts`}{The flat list of convex cells, as declared.}
#'   \item{`disjoint`}{`NULL`. A cache for a disjoint decomposition, filled by
#'   a later phase; nothing computes it yet.}
#'   \item{`triangulation`}{`NULL`. A cache for a simplicial decomposition, on
#'   the same terms.}
#' }
#' @examples
#' # The K = 3 plurality null: two overlapping sub-simplices.
#' union_region(
#'   simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'   simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#' )
#'
#' # Nesting is flattened, so these agree:
#' s <- simplex_region(vertices = diag(3))
#' h <- halfspace_region(normal = c(1, -1, 0))
#' n_parts(union_region(s, h))
#' n_parts(union_region(list(s, h)))
#' n_parts(union_region(union_region(s), list(h)))
#'
#' # One cell is already a region, so it is handed back as it came:
#' identical(union_region(s), s)
#' @export
union_region <- new_class(
  "union_region",
  parent = region,
  properties = list(
    parts = class_list,
    disjoint = class_any,
    triangulation = class_any
  ),
  constructor = function(...) {
    flat <- flatten_parts(list(...))
    if (length(flat) == 1L && S7_inherits(flat[[1L]], convex_region)) {
      return(flat[[1L]])
    }
    new_object(
      S7_object(),
      parts = flat,
      disjoint = NULL,
      triangulation = NULL
    )
  },
  validator = function(self) {
    if (length(self@parts) == 0L) {
      return("`parts` must be a non-empty list")
    }
    ok <- vapply(
      self@parts,
      \(p) S7_inherits(p, convex_region),
      logical(1)
    )
    if (!all(ok)) {
      return("every element of `parts` must be a `convex_region`")
    }
    # Ambient dimension only; shape and codimension are free. Cells of
    # differing ambient dimension have no common space to union in, and
    # comparing one against a parameter would silently recycle rather than
    # complain. Both dimensions are named, since neither is more wrong.
    dims <- vapply(self@parts, space_dim, integer(1))
    if (length(unique(dims)) > 1L) {
      return(paste0(
        "every element of `parts` must have the same ambient dimension; got ",
        paste(unique(dims), collapse = ", ")
      ))
    }
    NULL
  }
)


#' Flatten union-ish input into a list of convex parts
#'
#' Descends bare lists, unwraps unions into their own parts, and leaves anything
#' else alone as a leaf for the validator to name.
#' @keywords internal
#' @noRd
flatten_parts <- function(x) {
  if (S7_inherits(x, union_region)) {
    return(x@parts)
  }
  if (S7_inherits(x, convex_region)) {
    return(list(x))
  }
  if (is.list(x) && !S7_inherits(x)) {
    return(c(list(), unlist(lapply(x, flatten_parts), recursive = FALSE)))
  }
  list(x)
}


#' Coerce region-ish input to a [region]
#'
#' A [region] passes through untouched; a list becomes a [union_region] of its
#' elements, which for a one-element list is that element itself.
#' @param x A [region], or a list of them.
#' @return A [region].
#' @keywords internal
#' @noRd
as_region <- function(x) {
  if (S7_inherits(x, region)) x else union_region(x)
}


method(space_dim, union_region) <- function(space) {
  # The validator has already established that there is at least one part and
  # that they agree, so the first one speaks for all of them.
  space_dim(space@parts[[1L]])
}


method(contains, union_region) <- function(space, theta, tol = 1e-8) {
  any(vapply(space@parts, \(p) contains(p, theta, tol), logical(1)))
}

#' The refusal both representations owe a union
#' @keywords internal
#' @noRd
refuse_union <- function(what) {
  stop(
    "`",
    what,
    "()` is not defined for a `union_region`: a union is not an ",
    "intersection of half-spaces and has no single generator set. Take the ",
    "representation of each of `parts()` or `cells()` instead.",
    call. = FALSE
  )
}

method(h_rep, union_region) <- function(space) refuse_union("h_rep")
method(v_rep, union_region) <- function(space) refuse_union("v_rep")
method(q_hrep, union_region) <- function(space) refuse_union("q_hrep")
method(q_vrep, union_region) <- function(space) refuse_union("q_vrep")


method(is_empty, union_region) <- function(space) {
  all(vapply(space@parts, is_empty, logical(1)))
}


method(is_bounded, union_region) <- function(space) {
  all(vapply(space@parts, is_bounded, logical(1)))
}


method(parts, union_region) <- function(space) space@parts


#' @description A union's cells are its parts' cells, flattened: the parts are
#'   what was declared, the cells are what the algorithms run on.
#' @rdname cells
#' @usage NULL
method(cells, union_region) <- function(space) {
  unlist(lapply(space@parts, cells), recursive = FALSE)
}


#' The count of cells, as it should read in a message
#' @keywords internal
#' @noRd
parts_label <- function(n) sprintf("%d part%s", n, if (n == 1L) "" else "s")


#' @rdname union_region
#' @usage NULL
#' @export
method(print, union_region) <- function(x, ...) {
  n <- length(x@parts)
  cat("<", attr(S7_class(x), "name"), ">\n", sep = "")
  cat("  ", parts_label(n), ", dimension ", space_dim(x), "\n", sep = "")
  # The cells share a dimension, so the header has already said it and the
  # class name is all that is left to distinguish them. One line each while
  # that is readable, a tally beyond it: a triangulated null can hold hundreds
  # of cells, and listing them tells the reader nothing the tally does not.
  named <- vapply(x@parts, \(p) attr(S7_class(p), "name"), character(1))
  if (n <= 6L) {
    for (nm in named) {
      cat("    ", nm, "\n", sep = "")
    }
  } else {
    tally <- table(named)
    for (nm in names(tally)) {
      cat("    ", tally[[nm]], " x ", nm, "\n", sep = "")
    }
  }
  invisible(x)
}


#' @description `format()` gives the same summary on one line, without the class
#'   banner and the per-cell listing that `print()` adds.
#' @rdname union_region
#' @usage NULL
#' @export
method(format, union_region) <- function(x, ...) {
  sprintf(
    "%s: %s, dimension %d",
    attr(S7_class(x), "name"),
    parts_label(length(x@parts)),
    space_dim(x)
  )
}
