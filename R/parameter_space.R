#' @include sample_space.R
NULL

#' Parameter spaces
#'
#' A `parameter_space` is a convex set of parameter values. The same type serves
#' two roles. A family's own \eqn{\Theta}{Theta} is one, and so is each convex
#' piece \eqn{\Theta_{0i}}{Theta_0i} of a null hypothesis
#' \eqn{\Theta_0 = \bigcup_i \Theta_{0i}}{Theta_0 = union_i Theta_0i}. A null's
#' pieces may overlap.
#'
#' A `parameter_space` object encodes the geometry: it encodes dimension,
#' membership checking, projection, and a [chart()] that maps unconstrained
#' coordinates to the space. The oracle is written once against the chart, so
#' a new geometry needs no optimisation code.
#'
#' The contrast with [sample_space] is that these are searched over: outcomes
#' only ever need validating, whereas parameters need coordinates for a
#' optimiser to search over, which is what [chart()] supplies. There may be
#' null geometries that do not permit a [chart()], but these are currently
#' beyond the scope of this package.
#' @examples
#' # `parameter_space` is abstract; polytope_region(), simplex_region(),
#' # halfspace_region(), point_region() and real_region() subclass it, e.g.
#' s <- simplex_region(vertices = diag(3))
#' S7::S7_inherits(s, parameter_space)
#' space_dim(s)
#' @export
parameter_space <- new_class("parameter_space", abstract = TRUE)


#' Unconstrained coordinate chart for a parameter space
#'
#' Returns a list of closures that define mappings between the parameter space
#' and an unconstrained coordinate space, which is what lets BFGS run on a
#' constrained set.
#'
#' Charts for a compact space generally cover only the relative interior, so
#' an optimiser never exactly solves a maximum attained at a vertex or at
#' infinity, though at this point we are in the realm of numerical precision
#' anyway.
#' @param space A [parameter_space].
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
#' @param space A [parameter_space].
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
#' @param space A [parameter_space].
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
#' Defaults to projecting a reference point, typically the alternative's mean.
#' @param space A [parameter_space].
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


method(init_point, parameter_space) <- function(space, ref) project(space, ref)


#' The convex pieces a space decomposes into for optimisation
#'
#' Every space is its own only piece unless it says otherwise. A geometry that
#' is a union of convex cells (e.g. a triangulated polytope) returns those
#' cells here. Things that search or enclose then run on each convex piece,
#' without needing to know which geometry it came from.
#'
#' Pieces may overlap. This is fine for our case because a supremum over a
#' union is a supremum over any cover. Overlaps may cost time but the results
#' will still be valid.
#' @param space A [parameter_space].
#' @return A list of [parameter_space] objects whose union is `space`.
#' @examples
#' pieces(simplex_region(vertices = diag(3)))
#' @export
pieces <- new_generic("pieces", "space", function(space) S7::S7_dispatch())


method(pieces, parameter_space) <- function(space) list(space)


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
#' @param space A [parameter_space].
#' @param obj An [objective()].
#' @param seeds Optional `(d, m)` matrix of parameter-space points to seed from.
#'   Projected onto `space` before use, so points on other subnulls are fine.
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


# --- Polytope region ---------------------------------------------------------

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
  parent = parameter_space,
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
  parent = parameter_space,
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
  parent = parameter_space,
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


# --- Real region --------------------------------------------------------------

#' The whole of `R^d` as a parameter space
#'
#' The unconstrained case: \eqn{\Theta = \mathbb{R}^d}{Theta = R^d}, with the
#' identity chart. This is the parameter space of a [gaussian_family()], and the
#' way to say that a null places no constraint at all.
#'
#' Being unbounded it has no vertices, so like [halfspace_region()] it admits no
#' certified gap bound.
#'
#' @param d Integer dimension.
#' @return A `real_region`.
#' @examples
#' real_region(2L)
#' project(real_region(2L), c(3, -1))
#' @export
real_region <- new_class(
  "real_region",
  parent = parameter_space,
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


method(space_dim, real_region) <- function(space) as.integer(space@n_dim)


method(chart, real_region) <- function(space) {
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


method(project, real_region) <- function(space, theta) as.vector(theta)


method(contains, real_region) <- function(space, theta, tol = 1e-8) {
  length(theta) == as.integer(space@n_dim) && all(is.finite(theta))
}
