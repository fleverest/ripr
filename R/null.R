#' @include family.R
NULL

#' One convex piece of a null hypothesis
#'
#' The null is a finite union \eqn{\Theta_0 = \bigcup_i \Theta_{0i}}{Theta_0 = union_i Theta_0i}
#' of convex sets, and a `subnull` is one \eqn{\Theta_{0i}}{Theta_0i}. The pieces
#' may overlap; nothing assumes a partition.
#'
#' Not called a "face": a face of a polytope is the intersection with a
#' supporting hyperplane, hence a boundary object of codimension at least one,
#' whereas these are full-dimensional and may overlap. Each is a null hypothesis
#' in its own right, which is also how the per-piece e-values of a decomposition
#' are indexed.
#'
#' A subnull owns its geometry -- membership, projection, and a [chart()] giving
#' unconstrained coordinates -- and nothing else. The oracle is written once
#' against the chart, so a new subnull type needs no optimisation code.
#' @export
subnull <- new_class("subnull", abstract = TRUE)


#' Unconstrained coordinate chart for a subnull
#'
#' Returns a list of closures mapping between a subnull and an unconstrained
#' coordinate space, which is what lets BFGS run on a constrained set.
#'
#' Charts generally cover only the relative interior, so a maximum attained at a
#' vertex or at infinity is approached but never reached. That is why the search
#' below returns a lower bound on the true supremum.
#' @param subnull A [subnull].
#' @return A list of closures comprising:
#' \describe{
#'   \item{`n_par`}{dimension of the coordinate space.}
#'   \item{`to_theta(u)`}{coordinates to a parameter vector in the subnull.}
#'   \item{`to_theta_batch(u_mat)`}{`(n_par, N)` coordinates to `(d, N)`
#'   parameters.}
#'   \item{`from_theta(theta)`}{coordinates for a point in the subnull. A cheap
#'   inverse, not a metric projection; compose with [project()] first if `theta`
#'   may lie outside.}
#'   \item{`jacobian(u)`}{`(d, n_par)` derivative of `to_theta` at `u`.}
#'   \item{`seed(n)`}{`(n_par, n)` random coordinates for a multi-start search,
#'   drawn to suit the subnull's own geometry.}
#' }
#' @export
chart <- new_generic("chart", "subnull", function(subnull) S7::S7_dispatch())


#' Is a parameter in the subnull?
#' @param subnull A [subnull].
#' @param theta Parameter vector.
#' @param tol Tolerance.
#' @return `TRUE` or `FALSE`.
#' @export
contains <- new_generic(
  "contains",
  "subnull",
  function(subnull, theta, tol = 1e-8) S7::S7_dispatch()
)


#' Euclidean projection onto a subnull
#'
#' The closest point of the subnull to `theta`. Idempotent up to tolerance, and
#' its output always satisfies [contains()].
#' @param subnull A [subnull].
#' @param theta Parameter vector.
#' @return A parameter vector in the subnull.
#' @export
project <- new_generic(
  "project",
  "subnull",
  function(subnull, theta) S7::S7_dispatch()
)


#' A starting atom on the subnull
#'
#' Defaults to projecting a reference point, typically the alternative's mean.
#' @param subnull A [subnull].
#' @param ref Reference parameter vector.
#' @return A parameter vector in the subnull.
#' @export
init_point <- new_generic(
  "init_point",
  "subnull",
  function(subnull, ref) S7::S7_dispatch()
)


method(init_point, subnull) <- function(subnull, ref) project(subnull, ref)


#' Bundle an objective for [maximise_over()]
#'
#' @param value Function of a parameter vector returning a scalar.
#' @param grad Function of a parameter vector returning the gradient, length `d`.
#' @param value_batch Optional function of a `(d, N)` matrix returning `N`
#'   values. Used to score the seed grid; defaults to applying `value` per
#'   column, which is correct but slower.
#' @return A list to pass to [maximise_over()].
#' @export
objective <- function(value, grad, value_batch = NULL) {
  if (is.null(value_batch)) {
    value_batch <- function(theta_mat) {
      vapply(seq_len(ncol(theta_mat)), \(i) value(theta_mat[, i]), numeric(1))
    }
  }
  list(value = value, grad = grad, value_batch = value_batch)
}


#' Maximise an objective over a subnull
#'
#' Multi-start BFGS in the subnull's own [chart()]: seed coordinates are scored
#' in batch, the best `n_restarts` are refined, and the best refinement wins.
#' Written once and shared by every subnull type, since only the chart differs.
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
#' @param subnull A [subnull].
#' @param obj An [objective()].
#' @param seeds Optional `(d, m)` matrix of parameter-space points to seed from.
#'   Projected onto the subnull before use, so points on other subnulls are fine.
#' @param n_seeds Random seeds drawn from the chart.
#' @param n_restarts How many of the best seeds to refine.
#' @return `list(theta = , value = )` with `theta` in the subnull.
#' @export
maximise_over <- function(
  subnull,
  obj,
  seeds = NULL,
  n_seeds = 200L,
  n_restarts = 25L
) {
  ch <- chart(subnull)

  # A zero-dimensional chart has nothing to search: the subnull is a point.
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
      \(i) ch$from_theta(project(subnull, seeds[, i])),
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


# --- Simplex subnull ----------------------------------------------------------

#' Convex hull of a set of vertices
#'
#' The bounded case: a polytope given by its vertices (a V-representation),
#' parametrised by convex combinations of them. Every plurality subnull
#' \eqn{\{\theta : \theta_1 \le \theta_j\}}{{theta : theta_1 <= theta_j}} is of
#' this form.
#'
#' Holding vertices rather than constraints is what makes a certified bound
#' possible later: de Casteljau subdivision needs something to subdivide, and a
#' halfspace has no vertices.
#'
#' @param vertices `(d, V)` numeric matrix, one vertex per column.
#' @param pinv Optional precomputed `(V, d)` left pseudo-inverse.
#' @return A `simplex_null`.
#' @references
#' Beck, A. and Teboulle, M. (2009). A fast iterative shrinkage-thresholding
#' algorithm for linear inverse problems. *SIAM Journal on Imaging Sciences*
#' **2**(1), 183-202. \doi{10.1137/080716542}
#'
#' Duchi, J., Shalev-Shwartz, S., Singer, Y. and Chandra, T. (2008). Efficient
#' projections onto the l1-ball for learning in high dimensions. *Proceedings of
#' the 25th International Conference on Machine Learning*, 272-279.
#'
#' O'Donoghue, B. and Candes, E. (2015). Adaptive restart for accelerated
#' gradient schemes. *Foundations of Computational Mathematics* **15**, 715-732.
#' @export
simplex_null <- new_class(
  "simplex_null",
  parent = subnull,
  properties = list(
    vertices = class_any,
    n_vertices = class_numeric,
    pinv = class_any
  ),
  constructor = function(vertices, pinv = NULL) {
    if (!is.matrix(vertices) || ncol(vertices) == 0L) {
      stop(
        "`vertices` must be a matrix with one vertex per column.",
        call. = FALSE
      )
    }
    if (is.null(pinv)) {
      sv <- svd(vertices)
      keep <- sv$d > max(dim(vertices)) * .Machine$double.eps * max(sv$d)
      pinv <- sv$v[, keep, drop = FALSE] %*%
        (t(sv$u[, keep, drop = FALSE]) / sv$d[keep])
    }
    new_object(
      S7_object(),
      vertices = vertices,
      n_vertices = ncol(vertices),
      pinv = pinv
    )
  }
)


method(chart, simplex_null) <- function(subnull) {
  vertices <- subnull@vertices
  pinv <- subnull@pinv
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
      vapply(seq_len(n), \(i) softmax0_inv(alpha[, i]), numeric(n_v - 1L))
    }
  )
}


method(project, simplex_null) <- function(subnull, theta) {
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
  vertices <- subnull@vertices
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


method(contains, simplex_null) <- function(subnull, theta, tol = 1e-8) {
  max(abs(project(subnull, theta) - theta)) <= tol
}


# --- Halfspace subnull --------------------------------------------------------

#' Halfspace given by a normal and an offset
#'
#' The unbounded case:
#'   \eqn{\{\theta : a^\top \theta \le b\}}{{theta : a'theta <= b}}
#' (H-representation). The plurality subnull
#' \eqn{\{\theta : \theta_1 \le \theta_j\}}{{theta : theta_1 <= theta_j}} is
#' `normal = e_1 - e_j`, `offset = 0`.
#'
#' Coordinates are `(z, s)`: `z` positions a point on the bounding hyperplane in
#' an orthonormal basis, and `softplus(s) >= 0` is the distance inward.
#'
#' Having no vertices, a halfspace admits no certified gap bound. That is a
#' property of the representation rather than a missing feature -- see
#' [simplex_null()] for the bounded alternative.
#'
#' @param normal Normal vector `a`; must be non-zero.
#' @param offset Offset `b`.
#' @return A `halfspace_null`.
#' @export
halfspace_null <- new_class(
  "halfspace_null",
  parent = subnull,
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


method(chart, halfspace_null) <- function(subnull) {
  d <- length(subnull@normal)
  anchor <- subnull@anchor
  basis <- subnull@basis
  unit <- subnull@unit

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
      inward <- (subnull@offset - sum(subnull@normal * theta)) /
        sqrt(sum(subnull@normal^2))
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
        matrix(stats::rnorm(n * (d - 1L)), nrow = d - 1L),
        stats::rnorm(n, mean = -1, sd = 2)
      )
    }
  )
}


method(project, halfspace_null) <- function(subnull, theta) {
  slack <- sum(subnull@normal * theta) - subnull@offset
  if (slack <= 0) {
    return(theta)
  }
  theta - (slack / sum(subnull@normal^2)) * subnull@normal
}


method(contains, halfspace_null) <- function(subnull, theta, tol = 1e-8) {
  sum(subnull@normal * theta) <=
    subnull@offset + tol * sqrt(sum(subnull@normal^2))
}


# --- Singleton subnull --------------------------------------------------------

#' A single parameter point
#'
#' The degenerate convex set \eqn{\{\theta^*\}}{{theta*}}. Its role is to let a
#' likelihood ratio state the null it is valid for: \eqn{Q / P}{Q / P} is an
#' e-variable for \eqn{\{P\}}{{P}} and generally for nothing larger.
#'
#' @param theta The parameter vector.
#' @return A `singleton_null`.
#' @export
singleton_null <- new_class(
  "singleton_null",
  parent = subnull,
  properties = list(theta = class_numeric)
)

method(chart, singleton_null) <- function(subnull) {
  theta <- subnull@theta
  list(
    n_par = 0L,
    to_theta = function(u) theta,
    to_theta_batch = function(u_mat) matrix(theta, ncol = 1L),
    from_theta = function(theta) numeric(0L),
    jacobian = function(u) matrix(0, length(theta), 0L),
    seed = function(n) matrix(numeric(0L), nrow = 0L, ncol = n)
  )
}

method(project, singleton_null) <- function(subnull, theta) subnull@theta

method(contains, singleton_null) <- function(subnull, theta, tol = 1e-8) {
  max(abs(subnull@theta - theta)) <= tol
}


# --- The null hypothesis ------------------------------------------------------

#' A null hypothesis: a family together with its parameter region
#'
#' \eqn{H_0 = \{P_\theta : \theta \in \bigcup_i \Theta_{0i}\}}{H_0 = {P_theta : theta in union_i Theta_0i}}.
#' The null is a set of *distributions*, so it needs both the model and the
#' geometry; keeping them together means nothing downstream has to carry them as
#' separate arguments that could disagree.
#'
#' @param family A [sampling_family].
#' @param subnulls A non-empty list of [subnull] objects.
#' @return A `null_model`.
#' @export
null_model <- new_class(
  "null_model",
  properties = list(family = sampling_family, subnulls = class_list),
  validator = function(self) {
    if (length(self@subnulls) == 0L) {
      return("`subnulls` must be a non-empty list")
    }
    ok <- vapply(self@subnulls, \(s) S7_inherits(s, subnull), logical(1))
    if (!all(ok)) {
      return("every element of `subnulls` must be a `subnull`")
    }
    NULL
  }
)


#' Number of convex pieces in a null
#' @param null A [null_model].
#' @return Integer.
#' @export
n_subnulls <- function(null) length(null@subnulls)


#' Does any subnull contain this parameter?
#' @param null A [null_model].
#' @param theta Parameter vector.
#' @param tol Tolerance.
#' @return `TRUE` or `FALSE`.
#' @export
in_null <- function(null, theta, tol = 1e-8) {
  any(vapply(null@subnulls, \(s) contains(s, theta, tol), logical(1)))
}
