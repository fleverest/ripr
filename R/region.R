#' @include sample_space.R
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
#'   \item{`from_theta(theta)`}{coordinates for a point in the space. A cheap
#'   inverse, not a metric projection; compose with [project()] first if `theta`
#'   may lie outside.}
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
#'   Projected onto `space` before use, so points on other parts are fine.
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
      gradient = -as.vector(obj$grad(theta) %*% ch$jacobian(u))
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
    projected <- vapply(
      seq_len(ncol(seeds)),
      \(i) ch$from_theta(project(space, seeds[, i])),
      numeric(ch$n_par)
    )
    starts <- cbind(matrix(projected, nrow = ch$n_par), starts)
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


# --- Polytope region ----------------------------------------------------------

#' The left pseudo-inverse of a vertex matrix, for barycentric recovery
#' @keywords internal
#' @noRd
vertex_pinv <- function(vertices) {
  sv <- svd(vertices)
  keep <- sv$d > max(dim(vertices)) * .Machine$double.eps * max(sv$d)
  sv$v[, keep, drop = FALSE] %*% (t(sv$u[, keep, drop = FALSE]) / sv$d[keep])
}


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
#' @return A `polytope_region`.
#' @references
#'   \insertRef{BeckTeboulle2009}{ripr}
#'
#'   \insertRef{DuchiShalevShwartz2008}{ripr}
#'
#'   \insertRef{ODonoghueCandes2015}{ripr}
#' @examples
#' # A square in R^2
#' polytope_region(vertices = cbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)))
#' @export
polytope_region <- new_class(
  "polytope_region",
  parent = convex_region,
  properties = list(
    vertices = new_property(
      class_any,
      setter = function(self, value) {
        if (!is.matrix(value) || ncol(value) == 0L) {
          stop(
            "`vertices` must be a matrix with one vertex per column.",
            call. = FALSE
          )
        }
        if (!all(is.finite(value))) {
          stop("`vertices` must all be finite.", call. = FALSE)
        }
        attr(self, "vertices") <- value
        attr(self, "pinv") <- vertex_pinv(value)
        self
      }
    ),
    n_vertices = new_property(
      class_numeric,
      getter = function(self) ncol(self@vertices)
    ),
    pinv = class_any
  ),
  constructor = function(vertices) {
    new_object(S7_object(), vertices = vertices)
  }
)

method(space_dim, polytope_region) <- function(space) nrow(space@vertices)


method(chart, polytope_region) <- function(space) {
  vertices <- space@vertices
  pinv <- space@pinv
  n_v <- ncol(vertices)

  list(
    n_par = n_v - 1L,
    to_theta = function(u) as.vector(vertices %*% softmax0(u)),
    to_theta_batch = function(u_mat) {
      alpha <- vapply(
        seq_len(ncol(u_mat)),
        \(i) softmax0(u_mat[, i]),
        numeric(n_v)
      )
      vertices %*% matrix(alpha, nrow = n_v)
    },
    from_theta = function(theta) {
      alpha <- pmax(as.vector(pinv %*% theta), 0)
      softmax0_inv(alpha / sum(alpha))
    },
    jacobian = function(u) vertices %*% softmax0_jacobian(softmax0(u)),
    # Uniform Dirichlet over the vertices: already adapted to the geometry,
    # unlike Gaussian noise in an arbitrary coordinate system.
    seed = function(n) {
      g <- matrix(stats::rgamma(n * n_v, shape = 1), nrow = n_v)
      alpha <- div_by_col(g, colSums(g))
      matrix(
        vapply(seq_len(n), \(i) softmax0_inv(alpha[, i]), numeric(n_v - 1L)),
        nrow = n_v - 1L
      )
    }
  )
}


method(project, polytope_region) <- function(space, theta) {
  # Least squares over the vertex weights, constrained to the simplex:
  #   min_alpha ||V alpha - theta||^2  subject to  alpha in the simplex.
  #
  # Solved by FISTA (Beck and Teboulle, 2009): a projected gradient step,
  # accelerated by extrapolating along the previous step with the weight
  # (t_k - 1)/t_{k+1}. "Shrinkage-thresholding" refers to the proximal operator,
  # which here is the simplex projection. Acceleration makes the iterates
  # oscillate when the momentum term overshoots, so the momentum is reset
  # whenever the step points uphill (O'Donoghue and Candes, 2015).
  #
  # Chosen over the pseudo-inverse recovery in the chart because it handles
  # redundant vertex sets, and over a QP solver because the problem is tiny and
  # this needs no dependency.
  vertices <- space@vertices
  n_v <- ncol(vertices)
  lip <- max(svd(vertices)$d)^2
  vtv <- crossprod(vertices)
  vtt <- as.vector(crossprod(vertices, theta))

  alpha <- rep(1 / n_v, n_v)
  y <- alpha
  t_k <- 1
  for (i in seq_len(20000L)) {
    new_alpha <- project_simplex(y - (as.vector(vtv %*% y) - vtt) / lip)
    if (sum((y - new_alpha) * (new_alpha - alpha)) > 0) {
      y <- new_alpha
      t_k <- 1
    } else {
      t_new <- (1 + sqrt(1 + 4 * t_k^2)) / 2
      y <- new_alpha + ((t_k - 1) / t_new) * (new_alpha - alpha)
      t_k <- t_new
    }
    converged <- max(abs(new_alpha - alpha)) < 1e-14
    alpha <- new_alpha
    if (converged) break
  }
  as.vector(vertices %*% alpha)
}


method(contains, polytope_region) <- function(space, theta, tol = 1e-8) {
  max(abs(project(space, theta) - theta)) <= tol
}


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
    # `polytope_region`'s setter rejects non-finite vertices on every path that
    # assigns them, so this guards a hand-rolled `new_object()` only. It stays
    # first regardless: `NaN < 0` is `NA`, and `if (NA)` would surface as
    # "missing value where TRUE/FALSE needed" from inside S7.
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
      # Affine independence of the columns is linear independence of the edge
      # vectors from the first. Measured by the edge matrix's reciprocal
      # condition number rather than a determinant, which may not exist (the
      # vertex matrix need not be square) and would not be scale invariant if
      # it did: scaling every vertex by `c` scales `det` by `c^d`.
      #
      # The check above bounds `n_v - 1 <= d`, so `sv` has exactly `n_v - 1`
      # entries and its last is the one that decides rank.
      edges <- vertices[, -1L, drop = FALSE] - vertices[, 1L]
      sv <- svd(edges, nu = 0L, nv = 0L)$d
      if (sv[1L] <= 0) {
        return("every vertex is the same point, so they span nothing")
      }
      rcond <- sv[n_v - 1L] / sv[1L]
      if (rcond <= 1e-9) {
        return(paste0(
          "the vertices must span a non-degenerate simplex, but their ",
          "reciprocal condition number is ",
          format(rcond),
          ", so they are affinely dependent."
        ))
      }
    }
    NULL
  },
  constructor = function(vertices) {
    new_object(polytope_region(vertices = vertices))
  }
)


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
  parent = convex_region,
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
    d <- length(normal)
    nrm <- sqrt(sum(normal^2))
    basis <- if (d >= 2L) {
      qr.Q(qr(cbind(normal / nrm, diag(d))))[, 2:d, drop = FALSE]
    } else {
      NULL
    }
    new_object(
      S7_object(),
      normal = normal,
      offset = as.numeric(offset),
      unit = normal / nrm,
      anchor = offset * normal / nrm^2,
      basis = basis
    )
  }
)


method(space_dim, halfspace_region) <- function(space) {
  length(space@normal)
}


method(chart, halfspace_region) <- function(space) {
  d <- length(space@normal)
  anchor <- space@anchor
  basis <- space@basis
  unit <- space@unit

  on_plane <- function(z) {
    if (d >= 2L) anchor + as.vector(basis %*% z) else anchor
  }
  to_theta <- function(u) on_plane(u[-d]) - softplus(u[d]) * unit

  list(
    n_par = d,
    to_theta = to_theta,
    to_theta_batch = function(u_mat) {
      matrix(
        vapply(seq_len(ncol(u_mat)), \(i) to_theta(u_mat[, i]), numeric(d)),
        nrow = d
      )
    },
    from_theta = function(theta) {
      inward <- (space@offset - sum(space@normal * theta)) /
        sqrt(sum(space@normal^2))
      z <- if (d >= 2L) {
        as.vector(crossprod(basis, theta - anchor))
      } else {
        numeric(0L)
      }
      c(z, softplus_inv(inward))
    },
    jacobian = function(u) cbind(basis, -sigmoid(u[d]) * unit),
    # An unbounded set has no intrinsic scale, so seeds are Gaussian on the
    # hyperplane and biased toward the boundary, where the optimum usually sits.
    seed = function(n) {
      rbind(
        matrix(stats::rnorm(n * (d - 1L)), nrow = d - 1L, ncol = n),
        stats::rnorm(n, mean = -1, sd = 2)
      )
    }
  )
}


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
#' @param theta The parameter vector.
#' @return A `point_region`.
#' @examples
#' point_region(theta = c(0.5, 0.3, 0.2))
#' @export
point_region <- new_class(
  "point_region",
  parent = convex_region,
  properties = list(theta = class_numeric)
)

method(space_dim, point_region) <- function(space) length(space@theta)


method(chart, point_region) <- function(space) {
  theta <- space@theta
  list(
    n_par = 0L,
    to_theta = function(u) theta,
    to_theta_batch = function(u_mat) matrix(theta, ncol = 1L),
    from_theta = function(theta) numeric(0L),
    jacobian = function(u) matrix(0, length(theta), 0L),
    seed = function(n) matrix(numeric(0L), nrow = 0L, ncol = n)
  )
}

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
  parent = convex_region,
  properties = list(n_dim = class_numeric),
  constructor = function(d) {
    d <- as.integer(d)
    stopifnot(
      "`d` must be a single positive integer" = length(d) == 1L &&
        !is.na(d) &&
        d >= 1L
    )
    new_object(S7_object(), n_dim = d)
  }
)


method(space_dim, unconstrained_region) <- function(space) {
  as.integer(space@n_dim)
}


method(chart, unconstrained_region) <- function(space) {
  d <- as.integer(space@n_dim)
  list(
    n_par = d,
    to_theta = function(u) as.vector(u),
    to_theta_batch = function(u_mat) as.matrix(u_mat),
    from_theta = function(theta) as.vector(theta),
    jacobian = function(u) diag(d),
    # No intrinsic scale to adapt to, so standard normal is as good as anything.
    seed = function(n) matrix(stats::rnorm(n * d), nrow = d, ncol = n)
  )
}


method(project, unconstrained_region) <- function(space, theta) as.vector(theta)


method(contains, unconstrained_region) <- function(space, theta, tol = 1e-8) {
  length(theta) == as.integer(space@n_dim) && all(is.finite(theta))
}


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
