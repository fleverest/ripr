#' Face of a null-hypothesis region
#'
#' A `face` is one convex piece of a (generally non-convex, union-structured)
#' null hypothesis. Faces own their geometry (membership, projection,
#' parametrisation) and their optimisation oracle; algorithms only ever talk to
#' the generics below, so new face types drop in without touching algorithm
#' code.
#' @export
face <- new_class("face", abstract = TRUE)

#' Map face-local coordinates to a parameter vector on the face
#'
#' Contract: `alpha` lies in the face's coordinate domain (for polytope faces,
#' the `(V-1)`-simplex over vertices); the result is a parameter vector in the
#' face. The map must be affine with Jacobian [face_jacobian()].
#' @param face A `face`.
#' @param alpha Face-local coordinates.
#' @return Parameter vector on the face.
#' @keywords internal
parametrise <- new_generic("parametrise", "face", function(face, alpha) {
  S7::S7_dispatch()
})

#' Batched [parametrise()]: one input per row, one output per column
#'
#' Contract: `alpha_mat` is `(N, V)`; the result is `(d, N)`.
#' @param face A `face`.
#' @param alpha_mat `(N, V)` matrix of face-local coordinates.
#' @return `(d, N)` matrix of parameter vectors.
#' @keywords internal
parametrise_batch <- new_generic(
  "parametrise_batch",
  "face",
  function(face, alpha_mat) S7::S7_dispatch()
)

#' Constant Jacobian of [parametrise()]
#' @param face A `face`.
#' @return `(d, V)` Jacobian matrix.
#' @keywords internal
face_jacobian <- new_generic("face_jacobian", "face", function(face) {
  S7::S7_dispatch()
})

#' Recover face-local coordinates for a point on (or near) the face
#'
#' Contract: for `theta` on the face, returns coordinates `alpha` such that
#' `parametrise(face, alpha)` reproduces `theta` up to the recovery error of the
#' face's pseudo-inverse. This is the cheap seed-recovery map used by the EM
#' M-step; it is NOT a metric projection -- see [project()] for that.
#' @param face A `face`.
#' @param theta Parameter vector.
#' @return Face-local coordinates.
#' @keywords internal
face_coordinates <- new_generic(
  "face_coordinates",
  "face",
  function(face, theta) S7::S7_dispatch()
)

#' Metric projection onto the face
#'
#' Contract: returns the closest point of the face to `theta` (Euclidean
#' distance), computed in closed form where available. Idempotent up to
#' numerical tolerance, and its output always satisfies [contains()].
#' @param face A `face`.
#' @param theta Parameter vector.
#' @return The projected point, lying in the face.
#' @export
project <- new_generic("project", "face", function(face, theta) {
  S7::S7_dispatch()
})

#' Face membership test
#'
#' Contract: `TRUE` when `theta` lies within `tol` (sup-norm) of the face.
#' @param face A `face`.
#' @param theta Parameter vector.
#' @param tol Tolerance.
#' @return Logical scalar.
#' @export
contains <- new_generic(
  "contains",
  "face",
  function(face, theta, tol = 1e-8) {
    S7::S7_dispatch()
  }
)

#' Initial atom on the face given a reference point
#'
#' Contract: returns a point on the face suitable for initialising the first
#' atom, typically derived from the alternative's mean or mode.
#' @param face A `face`.
#' @param ref Reference point.
#' @return A point on the face.
#' @export
init_point <- new_generic("init_point", "face", function(face, ref) {
  S7::S7_dispatch()
})

#' Optimisation oracle over the face
#'
#' Contract: maximises the supplied objective over the face and returns
#' `list(theta = , value = )` with `theta` on the face. `objective` is a list of
#' closures over parameter vectors: `value(theta)`, `grad_theta(theta)`, and
#' `value_batch(theta_mat)`. Every implementation is a heuristic global search
#' over a generally non-convex objective, so the returned value is a *lower*
#' bound on the true face maximum, not the maximum itself. Callers that need an
#' upper bound on `sup G` -- certification, above all -- must account for that
#' separately; nothing here guarantees it.
#' @param face A `face`.
#' @param objective List of `value`, `grad_theta`, `value_batch` closures.
#' @param ... Ignored.
#' @param n_seeds Random seeds for the global search. Treated as a floor:
#'   implementations that stratify over `seed_centres` may draw more so that
#'   every centre receives a minimum share.
#' @param n_restarts Number of best seeds refined by BFGS.
#' @param seed_alpha Optional single seed for a local refinement (skips the
#'   global search).
#' @param seed_centres Optional `d x m` matrix of parameter-space points to
#'   centre the random search on -- in practice the current mixture's atoms.
#'   Points are projected onto the face before use, so centres lying on other
#'   faces are admissible. `NULL` leaves the choice to the implementation.
#'   Unbounded faces have no intrinsic scale, so without this hint their search
#'   is centred on an arbitrary anchor and can miss the optimum badly; bounded
#'   faces parametrised by their own vertices do not need it.
#' @return `list(theta = , value = )` with `theta` on the face.
#' @export
oracle <- new_generic(
  "oracle",
  "face",
  function(
    face,
    objective,
    ...,
    n_seeds = NULL,
    n_restarts = 25L,
    seed_alpha = NULL,
    seed_centres = NULL
  ) {
    S7::S7_dispatch()
  }
)

#' Union-structured null hypothesis region
#'
#' Holds the list of faces whose union is the null hypothesis. Every face stays
#' active on every oracle sweep.
#'
#' @param faces List of `face` objects.
#' @return A `null_region`.
#' @export
null_region <- new_class(
  "null_region",
  properties = list(
    faces = class_list
  ),
  validator = function(self) {
    if (length(self@faces) == 0L) {
      return("`faces` must be a non-empty list")
    }
    ok <- vapply(self@faces, function(f) S7_inherits(f, face), logical(1L))
    if (!all(ok)) {
      return("every element of `faces` must be a `face` object")
    }
    NULL
  }
)

# Euclidean projection of a vector onto the probability simplex (Duchi et al.
# 2008); deterministic and exact up to floating point.
project_simplex <- function(y) {
  u <- sort(y, decreasing = TRUE)
  css <- cumsum(u)
  rho <- max(which(u + (1 - css) / seq_along(u) > 0))
  tau <- (1 - css[rho]) / rho
  pmax(y + tau, 0)
}

#' Softmax reparametrisation of the simplex
#'
#' Maps an unconstrained vector `v` of length `d-1` to a point `alpha` on the
#' standard `(d-1)`-simplex via `softmax(c(0, v))`. The leading zero pins the
#' first coordinate, giving a bijection from `R^{d-1}` to the interior of
#' `Delta^{d-1}`. Suitable as an unconstrained parameterisation for BFGS.
#'
#' @param v Numeric vector of length `d-1`.
#' @return Numeric vector of length `d` summing to 1, with all entries > 0.
#' @keywords internal
#' @noRd
alpha_from_v <- function(v) {
  u <- c(0, v)
  e <- exp(u - max(u))
  e / sum(e)
}

#' Inverse softmax reparametrisation
#'
#' Maps a point `alpha` in the interior of `Delta^{d-1}` back to the
#' unconstrained vector `v` such that `alpha_from_v(v) == alpha`. Uses
#' log-ratios relative to the first coordinate, with an `eps` guard against
#' exact zeros.
#'
#' @param alpha Numeric vector summing to 1, with all entries >= 0.
#' @param eps Small constant added before taking logs to avoid `-Inf`.
#'   Default: `1e-12`.
#' @return Numeric vector of length `length(alpha) - 1`.
#' @keywords internal
#' @noRd
v_from_alpha <- function(alpha, eps = 1e-12) {
  log(alpha[-1L] + eps) - log(alpha[1L] + eps)
}

#' Jacobian of the softmax reparametrisation
#'
#' Computes `d(alpha)/d(v)` at the point corresponding to `alpha`, where
#' `v = v_from_alpha(alpha)`. The result is the `d x (d-1)` matrix obtained by
#' dropping the first column of the full `d x d` softmax Jacobian, using the
#' chain rule through `d(u)/d(v) = rbind(0, I_{d-1})`.
#'
#' @param alpha Numeric vector of length `d` summing to 1 (a simplex point).
#' @return Numeric matrix of dimension `d x (d-1)`.
#' @keywords internal
#' @noRd
softmax_jacobian <- function(alpha) {
  (diag(alpha) - outer(alpha, alpha))[, -1L, drop = FALSE]
}

#' Convex polytope face given by a vertex matrix
#'
#' The workhorse face type: a bounded convex polytope described by the columns
#' of `vertices`, parametrised by convex combinations. The vertex matrix is the
#' constant Jacobian, and the SVD-based left pseudo-inverse (computed once at
#' construction unless supplied) drives the cheap coordinate recovery.
#'
#' The oracle is a softmax + BFGS search with Dirichlet random restarts:
#' `n_seeds` uniform-Dirichlet points are scored with `objective$value_batch`,
#' the best `n_restarts` seed BFGS runs in the unconstrained softmax
#' parametrisation, and the best BFGS result wins. Supplying `seed_alpha` skips
#' the global search (single local refinement, as used by the EM M-step).
#'
#' @param vertices `(d, V)` numeric matrix; columns are the polytope vertices.
#' @param face_index Integer label carried through results.
#' @param init_point_fn Optional function of a reference point returning an
#'   initial atom on the face; defaults to [project()].
#' @param pinv Optional precomputed `(V, d)` left pseudo-inverse of `vertices`.
#' @return A `polytope_face`.
#' @export
polytope_face <- new_class(
  "polytope_face",
  parent = face,
  properties = list(
    vertices = class_any,
    face_index = class_numeric,
    n_vertices = class_numeric,
    pinv = class_any,
    init_point_fn = class_any
  ),
  constructor = function(
    vertices,
    face_index,
    init_point_fn = NULL,
    pinv = NULL
  ) {
    if (!is.matrix(vertices)) {
      stop("`vertices` must be a matrix with one vertex per column")
    }
    if (is.null(pinv)) {
      svd_V <- svd(vertices)
      tol <- max(dim(vertices)) * .Machine$double.eps * max(svd_V$d)
      pos <- svd_V$d > tol
      pinv <- svd_V$v[, pos, drop = FALSE] %*%
        (t(svd_V$u[, pos, drop = FALSE]) / svd_V$d[pos])
    }
    new_object(
      S7_object(),
      vertices = vertices,
      face_index = as.numeric(face_index),
      n_vertices = ncol(vertices),
      pinv = pinv,
      init_point_fn = init_point_fn
    )
  }
)

method(parametrise, polytope_face) <- function(face, alpha) {
  as.vector(face@vertices %*% alpha)
}

method(parametrise_batch, polytope_face) <- function(face, alpha_mat) {
  face@vertices %*% t(alpha_mat)
}

method(face_jacobian, polytope_face) <- function(face) {
  face@vertices
}

method(face_coordinates, polytope_face) <- function(face, theta) {
  alpha <- as.vector(face@pinv %*% theta)
  alpha <- pmax(alpha, 0)
  alpha / sum(alpha)
}

method(project, polytope_face) <- function(face, theta) {
  V <- face@vertices
  # Simplex-constrained least squares min_alpha ||V alpha - theta||^2 via FISTA
  # with gradient restart; handles redundant vertex sets that the pseudo-inverse
  # recovery cannot. The problem is tiny, so a generous iteration cap is cheap.
  n_v <- ncol(V)
  lip <- max(svd(V)$d)^2
  vtv <- crossprod(V)
  vtt <- as.vector(crossprod(V, theta))
  alpha <- rep(1 / n_v, n_v)
  y <- alpha
  t_k <- 1
  for (i in seq_len(20000L)) {
    grad_y <- as.vector(vtv %*% y) - vtt
    alpha_new <- project_simplex(y - grad_y / lip)
    if (sum((y - alpha_new) * (alpha_new - alpha)) > 0) {
      # Momentum points uphill: restart it (O'Donoghue & Candes 2015).
      y <- alpha_new
      t_k <- 1
    } else {
      t_new <- (1 + sqrt(1 + 4 * t_k^2)) / 2
      y <- alpha_new + ((t_k - 1) / t_new) * (alpha_new - alpha)
      t_k <- t_new
    }
    if (max(abs(alpha_new - alpha)) < 1e-14) {
      alpha <- alpha_new
      break
    }
    alpha <- alpha_new
  }
  as.vector(V %*% alpha)
}

method(contains, polytope_face) <- function(face, theta, tol = 1e-8) {
  max(abs(project(face, theta) - theta)) <= tol
}

method(init_point, polytope_face) <- function(face, ref) {
  if (!is.null(face@init_point_fn)) {
    return(face@init_point_fn(ref))
  }
  project(face, ref)
}

# Chain-rule wrapper: (value, grad_theta) -> BFGS objective in softmax coords.
# Caches the last (input, output) pair since `optim` calls `fn` and `gr`
# separately at the same point.
make_face_objective <- function(face, value_fn, grad_theta_fn) {
  J <- face@vertices
  last_v <- NULL
  last_result <- NULL
  function(v) {
    if (!is.null(last_v) && identical(v, last_v)) {
      return(last_result)
    }
    alpha <- alpha_from_v(v)
    theta <- parametrise(face, alpha)
    val <- value_fn(theta)
    grad_theta <- grad_theta_fn(theta)
    grad_v <- as.vector(grad_theta %*% J %*% softmax_jacobian(alpha))
    last_v <<- v
    last_result <<- list(value = -val, gradient = -grad_v)
    last_result
  }
}

method(oracle, polytope_face) <- function(
  face,
  objective,
  ...,
  n_seeds = NULL,
  n_restarts = 25L,
  # A polytope face is parametrised by its own vertices, so Dirichlet seeds are
  # already adapted to its geometry; `seed_centres` is accepted for signature
  # compatibility and deliberately unused.
  seed_alpha = NULL,
  seed_centres = NULL
) {
  n_vertices <- face@n_vertices
  obj_and_grad <- make_face_objective(
    face,
    objective$value,
    objective$grad_theta
  )

  run_bfgs <- function(init_alpha, fallback_value) {
    init_alpha <- pmax(init_alpha, 1e-8)
    init_alpha <- init_alpha / sum(init_alpha)
    v0 <- v_from_alpha(init_alpha)
    tryCatch(
      optim(
        v0,
        fn = function(v) obj_and_grad(v)$value,
        gr = function(v) obj_and_grad(v)$gradient,
        method = "BFGS"
      ),
      error = function(e) list(par = v0, value = fallback_value)
    )
  }

  # Local optimisation from a single seed.
  if (!is.null(seed_alpha)) {
    res <- run_bfgs(seed_alpha, fallback_value = Inf)
    return(list(
      theta = parametrise(face, alpha_from_v(res$par)),
      value = -res$value
    ))
  }

  # Global (hopefully) optimisation via random Dirichlet restarts.
  alpha_mat <- matrix(
    rgamma(n_seeds * n_vertices, shape = 1),
    nrow = n_seeds,
    ncol = n_vertices
  )
  alpha_mat <- alpha_mat / rowSums(alpha_mat)

  neg_obj_grid <- objective$value_batch(parametrise_batch(face, alpha_mat))
  top_idx <- order(neg_obj_grid, decreasing = TRUE)[seq_len(min(
    n_restarts,
    length(neg_obj_grid)
  ))]

  best <- list(par = NULL, value = Inf)
  for (idx in top_idx) {
    res <- run_bfgs(alpha_mat[idx, ], fallback_value = -neg_obj_grid[idx])
    if (res$value < best$value) best <- res
  }
  list(
    theta = parametrise(face, alpha_from_v(best$par)),
    value = -best$value
  )
}

#' Normalise a faces argument into a list of `face` objects
#'
#' Accepts a `null_region`, a list of `face` objects, or a list of plurality
#' face-descriptor lists (whose precomputed jacobian/pinv/init closures are
#' wrapped without recomputation).
#' @param x A `null_region`, list of faces, or list of descriptor lists.
#' @return A list of `face` objects.
#' @keywords internal
as_faces <- function(x) {
  if (S7_inherits(x, null_region)) {
    return(x@faces)
  }
  if (!is.list(x) || length(x) == 0L) {
    stop("`faces` must be a null_region or a non-empty list")
  }
  lapply(x, function(f) {
    if (S7_inherits(f, face)) {
      return(f)
    }
    if (is.list(f) && all(c("jacobian", "pinv", "n_vertices") %in% names(f))) {
      return(polytope_face(
        vertices = f$jacobian,
        face_index = if (is.null(f$face_index)) NA_real_ else f$face_index,
        init_point_fn = f$init_point,
        pinv = f$pinv
      ))
    }
    stop("cannot interpret an element of `faces` as a face")
  })
}

# =============================================================================
# Half-space faces (Gaussian nulls)
# =============================================================================

softplus <- function(s) ifelse(s > 30, s, log1p(exp(s)))
sigmoid_scalar <- function(s) 1 / (1 + exp(-s))
softplus_inv <- function(t) {
  t <- pmax(t, 1e-8)
  ifelse(t > 30, t, log(expm1(t)))
}

#' Half-space face `{theta : <v, theta> <= c}`
#'
#' The Gaussian null building block. Projection and membership are closed form.
#' The face is parametrised over its full interior by
#' `theta(z, s) = anchor + B z - softplus(s) * v / ||v||`, where `anchor` is the
#' boundary projection of the origin, `B` an orthonormal basis of the boundary
#' hyperplane, and `softplus(s) >= 0` the distance into the half-space; the
#' coordinate vector is `c(z, s)`. The oracle runs BFGS over these coordinates
#' from Gaussian random seeds biased toward the boundary.
#'
#' @param v Normal vector of the constraint.
#' @param c Offset: the face is `{theta : <v, theta> <= c}`.
#' @param face_index Integer label carried through results.
#' @return A `halfspace_face`.
#' @export
halfspace_face <- new_class(
  "halfspace_face",
  parent = face,
  properties = list(
    v = class_numeric,
    c = class_numeric,
    v_unit = class_numeric,
    anchor = class_numeric,
    basis = class_any,
    face_index = class_numeric
  ),
  constructor = function(v, c, face_index = NA_real_) {
    v <- as.numeric(v)
    if (all(v == 0)) {
      stop("`v` must be a non-zero vector")
    }
    nv <- sqrt(sum(v^2))
    v_unit <- v / nv
    anchor <- c * v / nv^2
    d <- length(v)
    basis <- if (d >= 2L) {
      full_q <- qr.Q(qr(cbind(v_unit, diag(d))))
      full_q[, 2:d, drop = FALSE]
    } else {
      NULL
    }
    new_object(
      S7_object(),
      v = v,
      c = as.numeric(c),
      v_unit = v_unit,
      anchor = anchor,
      basis = basis,
      face_index = as.numeric(face_index)
    )
  }
)

method(parametrise, halfspace_face) <- function(face, alpha) {
  d <- length(face@v)
  s <- alpha[d]
  boundary <- if (d >= 2L) {
    face@anchor + as.vector(face@basis %*% alpha[-d])
  } else {
    face@anchor
  }
  boundary - softplus(s) * face@v_unit
}

method(parametrise_batch, halfspace_face) <- function(face, alpha_mat) {
  d <- length(face@v)
  matrix(
    vapply(
      seq_len(nrow(alpha_mat)),
      function(i) parametrise(face, alpha_mat[i, ]),
      numeric(d)
    ),
    nrow = d
  )
}

method(face_coordinates, halfspace_face) <- function(face, theta) {
  t_in <- (face@c - sum(face@v * theta)) / sqrt(sum(face@v^2))
  z <- if (length(face@v) >= 2L) {
    as.vector(crossprod(face@basis, theta - face@anchor))
  } else {
    numeric(0L)
  }
  c(z, softplus_inv(t_in))
}

method(project, halfspace_face) <- function(face, theta) {
  slack <- sum(face@v * theta) - face@c
  if (slack <= 0) {
    return(theta)
  }
  theta - (slack / sum(face@v^2)) * face@v
}

method(contains, halfspace_face) <- function(face, theta, tol = 1e-8) {
  sum(face@v * theta) <= face@c + tol * sqrt(sum(face@v^2))
}

method(init_point, halfspace_face) <- function(face, ref) {
  project(face, ref)
}

method(oracle, halfspace_face) <- function(
  face,
  objective,
  ...,
  n_seeds = NULL,
  n_restarts = 25L,
  seed_alpha = NULL,
  seed_centres = NULL
) {
  d <- length(face@v)
  n_par <- d # (d - 1) boundary coordinates plus one slack coordinate

  last_u <- NULL
  last_result <- NULL
  obj_and_grad <- function(u) {
    if (!is.null(last_u) && identical(u, last_u)) {
      return(last_result)
    }
    theta <- parametrise(face, u)
    val <- objective$value(theta)
    grad_theta <- objective$grad_theta(theta)
    jac <- cbind(face@basis, -sigmoid_scalar(u[n_par]) * face@v_unit)
    grad_u <- as.vector(grad_theta %*% jac)
    last_u <<- u
    last_result <<- list(value = -val, gradient = -grad_u)
    last_result
  }

  run_bfgs <- function(u0, fallback_value) {
    tryCatch(
      optim(
        u0,
        fn = function(u) obj_and_grad(u)$value,
        gr = function(u) obj_and_grad(u)$gradient,
        method = "BFGS"
      ),
      error = function(e) list(par = u0, value = fallback_value)
    )
  }

  if (!is.null(seed_alpha)) {
    res <- run_bfgs(seed_alpha, fallback_value = Inf)
    return(list(theta = parametrise(face, res$par), value = -res$value))
  }

  # Boundary-hyperplane coordinates of the seed centres. A half-space is
  # unbounded and so has no intrinsic scale: centring on `anchor` (the boundary
  # point nearest the origin) with unit spread is arbitrary, and misses the
  # optimum badly whenever it sits more than a couple of units away. Centring on
  # the mixture's own atoms removes that failure mode.
  z_centres <- if (is.null(seed_centres)) {
    matrix(0, nrow = n_par - 1L, ncol = 1L)
  } else {
    seed_centres <- as.matrix(seed_centres)
    # Points on other faces are admissible: project first, then take
    # coordinates (the projection of a point already in the face is itself).
    zc <- vapply(
      seq_len(ncol(seed_centres)),
      function(k) face_coordinates(face, project(face, seed_centres[, k]))[
        seq_len(n_par - 1L)
      ],
      numeric(max(n_par - 1L, 1L))
    )
    zc <- matrix(zc, nrow = n_par - 1L)
    zc <- zc[, is.finite(colSums(zc)), drop = FALSE]
    # The anchor is always a centre. Two reasons: with a single atom there is no
    # centre cloud to calibrate a spread from, and without this the search would
    # be *narrower* than the anchor-only default and blind to the region between
    # the support and the anchor. Including it makes the seeded region a strict
    # superset of the unhinted one, and makes the spread scale with how far the
    # support sits from the anchor.
    cbind(zc, matrix(0, nrow = n_par - 1L, ncol = 1L))
  }
  if (ncol(z_centres) == 0L) {
    z_centres <- matrix(0, nrow = n_par - 1L, ncol = 1L)
  }
  m <- ncol(z_centres)

  # Spread: a heavy-tailed radius calibrated so `p_out` of the seeds around each
  # centre land beyond the radius of the centre cloud. Cauchy rather than
  # Gaussian because the screen is cheap and the polish budget capped, which
  # makes far-out seeds nearly free insurance against a mis-estimated scale.
  p_out <- 0.30
  centroid <- rowMeans(z_centres)
  radius <- if (m > 1L) {
    max(sqrt(colSums((z_centres - centroid)^2)))
  } else {
    0
  }
  scale_c <- max(radius, 1) / tan((pi / 2) * (1 - p_out))

  # `n_seeds` is a floor: every centre gets a minimum share, so a growing
  # support is searched proportionally harder.
  seeds_per_centre_min <- 8L
  n_total <- max(as.integer(n_seeds), m * seeds_per_centre_min)
  stratum <- rep(seq_len(m), length.out = n_total)

  # Isotropic direction x half-Cauchy radius: exact control of the proportion
  # outside the cloud in any dimension, unlike per-coordinate draws.
  n_bnd <- n_par - 1L
  offset <- if (n_bnd >= 1L) {
    g <- matrix(rnorm(n_total * n_bnd), nrow = n_total, ncol = n_bnd)
    norms <- sqrt(rowSums(g^2))
    norms[norms == 0] <- 1
    (g / norms) * abs(rcauchy(n_total, location = 0, scale = scale_c))
  } else {
    matrix(0, nrow = n_total, ncol = 0L)
  }
  u_mat <- cbind(
    t(z_centres)[stratum, , drop = FALSE] + offset,
    rnorm(n_total) - 2 # slack: softplus(N(-2, 1)) hugs the boundary
  )

  theta_mat <- parametrise_batch(face, u_mat)
  vals <- objective$value_batch(theta_mat)
  vals[!is.finite(vals)] <- -Inf

  # Stratified polish: the best seed from each centre first, so no atom's
  # neighbourhood is dropped by a noisy screen, then fill by global value.
  n_pol <- min(n_restarts, n_total)
  per_centre <- vapply(
    seq_len(m),
    function(k) {
      idx <- which(stratum == k)
      idx[which.max(vals[idx])]
    },
    integer(1L)
  )
  per_centre <- per_centre[order(vals[per_centre], decreasing = TRUE)]
  top_idx <- utils::head(per_centre, n_pol)
  if (length(top_idx) < n_pol) {
    rest <- setdiff(order(vals, decreasing = TRUE), top_idx)
    top_idx <- c(top_idx, utils::head(rest, n_pol - length(top_idx)))
  }

  best <- list(par = NULL, value = Inf)
  for (idx in top_idx) {
    res <- run_bfgs(u_mat[idx, ], fallback_value = -vals[idx])
    if (res$value < best$value) best <- res
  }
  list(theta = parametrise(face, best$par), value = -best$value)
}
