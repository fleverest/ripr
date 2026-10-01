#' @include region.R
#' @include polyhedra.R
NULL

# --- Polyhedron region --------------------------------------------------------

#' Assemble and check a generator triple
#'
#' Shapes input into `list(v, r, l)`: `NULL` blocks become empty and a cone is
#' anchored at the origin. Content checks are the validator's job.
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
      "generators must be matrices, one generator per row.",
      call. = FALSE
    )
  }
  d <- ncol(given[[1L]])
  as_block <- function(x) {
    if (is.null(x)) no_generators(d) else x
  }
  with_origin_vertex(
    list(v = as_block(vertices), r = as_block(rays), l = as_block(lines))
  )
}


#' A convex polyhedron given by its generators
#'
#' The Minkowski--Weyl form
#' \eqn{\mathrm{conv}(V) + \mathrm{cone}(R) + \mathrm{span}(L)}{conv(V) + cone(R) + span(L)}
#' (Ziegler 1995, Theorem 1.2), base of every convex region in the package.
#' [polytope_region()], [simplex_region()], [halfspace_region()],
#' [point_region()] and [real_region()] are friendlier special cases; use
#' [h_region()] to start from a half-space description.
#'
#' The [chart()] is \eqn{\theta(u) = a V + z L + c R}{theta(u) = a V + z L + c R}
#' with `u = (a, z, c)` (vertex weights, lineality, rays), subject to
#' `a >= 0`, `sum(a) = 1` and `c >= 0`. A lone vertex contributes no
#' coordinate, so `n_par = nv + nl + nr` when `nv > 1` and `nl + nr`
#' otherwise.
#'
#' @param vertices `(nv, d)` numeric matrix, one vertex per row, or `NULL`
#'   for a cone anchored at the origin.
#' @param rays `(nr, d)` numeric matrix of recession directions, one per row,
#'   or `NULL`.
#' @param lines `(nl, d)` numeric matrix whose rows span the lineality space,
#'   or `NULL`.
#' @param .hv Internal use only.
#' @return A `polyhedron_region`.
#' @section Properties:
#' \describe{
#'   \item{`generators`}{`list(v, r, l)` of numeric matrices, one generator
#'   per row. Read-only: derived from the static underlying record.}
#'   \item{`facets`}{The half-space description `(A, B, eq)`, read as
#'   `a %*% theta <= b` for each row, with `eq` flagging when equality
#'   should hold rather than inequality for each row. Read-only.}
#'   \item{`hv`}{Internal record of both descriptions in double and exact
#'   rational form, computed once at construction.}
#' }
#' @references
#'   \insertRef{Ziegler1995}{ripr}
#' @examples
#' # The halfspace `{theta_1 <= 0}` in R^2, by hand:
#' polyhedron_region(
#'   vertices = matrix(c(0, 0), nrow = 1),
#'   rays = matrix(c(-1, 0), nrow = 1),
#'   lines = matrix(c(0, 1), nrow = 1)
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
      class_list,
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
    g <- make_generators(vertices, rays, lines)
    new_object(S7_object(), hv = hv_fill(v = g))
  },
  validator = function(self) {
    g <- self@generators
    if (!identical(names(g), c("v", "r", "l"))) {
      return("`generators` must be a list with elements `v`, `r` and `l`")
    }
    if (!all(vapply(g, \(x) is.matrix(x) && is.numeric(x), logical(1)))) {
      return("every generator block must be a numeric matrix")
    }
    dims <- vapply(g, ncol, integer(1))
    if (length(unique(dims)) > 1L) {
      return(paste0(
        "generator blocks disagree on the ambient dimension: ",
        paste(dims, collapse = ", ")
      ))
    }
    if (!all(vapply(g, \(x) all(is.finite(x)), logical(1)))) {
      return("every generator coordinate must be finite")
    }
    if (nrow(g$v) == 0L) {
      return("`generators$v` must hold at least one point")
    }
    if (any(rowSums(rbind(g$r, g$l) != 0) == 0L)) {
      return("every ray and line must be nonzero")
    }
    f <- self@facets
    if (!all(c("a", "b", "eq") %in% names(f))) {
      return("`facets` must be a list with `a`, `b` and `eq`")
    }
    if (!is.matrix(f$a) || ncol(f$a) != ncol(g$v)) {
      return("`facets$a` must have one column per ambient dimension")
    }
    if (length(f$b) != nrow(f$a) || length(f$eq) != nrow(f$a)) {
      return("`facets` must have one `b` and one `eq` entry per row of `a`")
    }
    NULL
  }
)


#' A convex polyhedron given by its half-space description
#'
#' The set `{theta : a %*% theta <= b}`, with `eq` flagging equality rows
#' (e.g. `sum(theta) == 1`). The rows are kept as the region's facets; the
#' generators are derived exactly.
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
  polyhedron_region(.hv = hv_fill(h = h))
}


method(space_dim, polyhedron_region) <- function(space) {
  ncol(space@generators$v)
}


# --- Representations ----------------------------------------------------------

#' Refuse anything but a `polyhedron_region` a representation
#' @keywords internal
#' @noRd
check_polyhedron <- function(space, what) {
  if (S7_inherits(space, polyhedron_region)) {
    return(invisible(space))
  }
  got <- if (S7_inherits(space)) class_name(space) else class(space)
  stop(
    "`",
    what,
    "()` is defined only for a `polyhedron_region`, not a `",
    got[[1L]],
    "`. Take the representation of each of `parts()` or `cells()` instead.",
    call. = FALSE
  )
}


#' The half-space description of a region
#'
#' `list(a, b, eq)` read as `a %*% theta <= b`, `eq` marking equality rows.
#' @keywords internal
#' @noRd
h_rep <- function(space) check_polyhedron(space, "h_rep")@facets


#' The generator description of a region
#'
#' `list(v, r, l)`, one generator per row, as declared (redundant vertices
#' stay). Rays are determined only modulo the lineality space.
#' @keywords internal
#' @noRd
v_rep <- function(space) check_polyhedron(space, "v_rep")@generators


#' The rational H- and V-representations of a region
#'
#' Exact rcdd matrices under `h_rep()`/`v_rep()`. For [real_region] `q_hrep()`
#' is the row `0 . x <= 1` while `h_rep()` has no rows.
#' @keywords internal
#' @noRd
q_hrep <- function(space) check_polyhedron(space, "q_hrep")@hv@qh


#' @keywords internal
#' @noRd
q_vrep <- function(space) check_polyhedron(space, "q_vrep")@hv@qv


method(is_empty, polyhedron_region) <- function(space) {
  q_is_empty(q_hrep(space))
}


method(is_bounded, polyhedron_region) <- function(space) {
  nrow(space@generators$r) == 0L && nrow(space@generators$l) == 0L
}


method(contains, polyhedron_region) <- function(space, theta, tol = 1e-8) {
  # Normalised by row norm so `tol` is a distance on every facet.
  h <- h_rep(space)
  slack <- as.numeric(h$a %*% theta) - h$b
  scale <- sqrt(rowSums(h$a^2))
  all(ifelse(h$eq, abs(slack) <= tol * scale, slack <= tol * scale))
}


#' Solve for a point's least-squares weights over a generator triple
#'
#' Minimises `|| a V + z L + c R - theta ||^2` over `a` in the simplex,
#' `c >= 0`; the image `theta_hat` is the Euclidean projection, serving both
#' [project()] and `from_theta()`. Tries the unconstrained least-squares
#' solution first (exact for any point in a simplex cell), else solves the
#' quadratic programme with `quadprog` and polishes its answer by re-solving
#' exactly on the generators it kept.
#' @return `list(a, z, c, theta_hat)`.
#' @keywords internal
#' @noRd
generator_weights <- function(g, theta, tol = 1e-9) {
  n_v <- nrow(g$v)
  n_r <- nrow(g$r)
  exact <- exact_weights(g, theta, seq_len(n_v), seq_len(n_r), tol)
  if (!is.null(exact)) {
    return(exact)
  }

  # The constrained problem: `sum(a) == 1`, then `a >= 0`, `c >= 0`. A
  # redundant generator set makes the quadratic term singular, which
  # `quadprog` refuses, so it gets a ridge far below `tol`.
  m <- t(rbind(g$v, g$l, g$r))
  n <- ncol(m)
  i_v <- seq_len(n_v)
  i_r <- n - n_r + seq_len(n_r)
  mtm <- crossprod(m)
  x <- quadprog::solve.QP(
    Dmat = mtm + diag(1e-12 * max(1, diag(mtm)), n),
    dvec = as.numeric(crossprod(m, theta)),
    Amat = cbind(as.numeric(seq_len(n) %in% i_v), diag(n)[, c(i_v, i_r)]),
    bvec = c(1, rep(0, n_v + n_r)),
    meq = 1L
  )$solution

  # The solver leaves weights of order its tolerance on generators that should
  # be exactly zero; dropping them and re-solving on the rest is exact.
  keep_v <- which(x[i_v] > 1e-7)
  keep_r <- which(x[i_r] > 1e-7)
  polished <- exact_weights(g, theta, keep_v, keep_r, tol)
  if (!is.null(polished)) {
    return(polished)
  }
  a <- pmax(x[i_v], 0)
  a <- a / sum(a)
  z <- x[n_v + seq_len(nrow(g$l))]
  cc <- pmax(x[i_r], 0)
  list(a = a, z = z, c = cc, theta_hat = generator_image(g, a, z, cc))
}


#' Least squares over a subset of the vertices and rays, or `NULL` if infeasible
#'
#' Substitutes `a = (1 - sum(b), b)` so `sum(a) = 1` holds, solves for `b`, and
#' accepts the result only if `a, c >= 0` anyway. Weights come back full length,
#' zero off the subset.
#' @keywords internal
#' @noRd
exact_weights <- function(g, theta, keep_v, keep_r, tol) {
  v <- g$v[keep_v, , drop = FALSE]
  r <- g$r[keep_r, , drop = FALSE]
  n_v <- nrow(v)
  n_l <- nrow(g$l)
  m <- t(rbind(add_by_col(v[-1L, , drop = FALSE], -v[1L, ]), g$l, r))
  x <- min_norm_solve(m, theta - v[1L, ])
  b <- x[seq_len(n_v - 1L)]
  a_kept <- c(1 - sum(b), b)
  c_kept <- x[(n_v - 1L) + n_l + seq_along(keep_r)]
  if (any(a_kept < -tol) || any(c_kept < -tol)) {
    return(NULL)
  }
  a <- numeric(nrow(g$v))
  a[keep_v] <- pmax(a_kept, 0) / sum(pmax(a_kept, 0))
  cc <- numeric(nrow(g$r))
  cc[keep_r] <- pmax(c_kept, 0)
  z <- x[(n_v - 1L) + seq_len(n_l)]
  list(a = a, z = z, c = cc, theta_hat = generator_image(g, a, z, cc))
}


#' `a V + z L + c R`, the point a set of generator weights describes
#' @keywords internal
#' @noRd
generator_image <- function(g, a, z, cc) {
  as.numeric(a %*% g$v + z %*% g$l + cc %*% g$r)
}


#' Minimum-norm least-squares solution of `m x = rhs`
#'
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
  n_v <- nrow(v)
  n_l <- nrow(l)
  n_r <- nrow(r)
  # A lone vertex contributes no coordinate: `theta = v_1 + z L + c R`.
  free_v <- if (n_v > 1L) n_v else 0L
  i_v <- seq_len(free_v)
  base <- if (free_v == 0L) as.numeric(v[1L, ]) else numeric(ncol(v))
  # `theta = base + u gens`.
  gens <- rbind(v[i_v, , drop = FALSE], l, r)
  jac <- t(gens)
  n_par <- free_v + n_l + n_r

  list(
    n_par = n_par,
    to_theta = function(u) base + as.numeric(u %*% gens),
    to_theta_batch = function(u_mat) add_by_col(u_mat %*% gens, base),
    from_theta = function(theta) {
      w <- generator_weights(g, theta)
      c(if (free_v > 0L) w$a, w$z, w$c)
    },
    jacobian = function(u) jac,
    lower = c(rep(0, free_v), rep(-Inf, n_l), rep(0, n_r)),
    heq = if (free_v > 0L) {
      function(u) sum(u[i_v]) - 1
    },
    heqjac = if (free_v > 0L) {
      function(u) matrix(c(rep(1, free_v), rep(0, n_l + n_r)), nrow = 1L)
    },
    # a ~ uniform Dirichlet, z ~ N(0, 1), c ~ Exp(1).
    seed = function(n) {
      u_v <- if (free_v > 0L) {
        dirichlet_draws(rep(1, n_v), n)
      } else {
        matrix(numeric(0), nrow = n, ncol = 0L)
      }
      cbind(
        u_v,
        matrix(stats::rnorm(n * n_l), nrow = n, ncol = n_l, byrow = TRUE),
        matrix(
          stats::rgamma(n * n_r, shape = 1),
          nrow = n,
          ncol = n_r,
          byrow = TRUE
        )
      )
    }
  )
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
#' Used by `format()` and `print()`.
#' @keywords internal
#' @noRd
region_phrase <- new_generic("region_phrase", "space", function(space) {
  S7::S7_dispatch()
})


method(region_phrase, polyhedron_region) <- function(space) {
  g <- space@generators
  counts <- c(
    count_label(nrow(g$v), "vertex", "vertices"),
    if (nrow(g$r) > 0L) count_label(nrow(g$r), "ray"),
    if (nrow(g$l) > 0L) count_label(nrow(g$l), "line")
  )
  sprintf("%s in R^%d", paste(counts, collapse = ", "), space_dim(space))
}


#' @rdname polyhedron_region
#' @usage NULL
#' @export
method(format, polyhedron_region) <- function(x, ...) {
  sprintf("%s: %s", class_name(x), region_phrase(x))
}


#' @description `print()` summarises generator and facet counts, and lists the
#'   vertices of a small polytope.
#' @rdname polyhedron_region
#' @usage NULL
#' @export
method(print, polyhedron_region) <- function(x, ...) {
  cat("<", class_name(x), ">\n", sep = "")
  cat(
    "  ",
    region_phrase(x),
    if (!is_bounded(x) && !S7_inherits(x, real_region)) ", unbounded",
    "\n",
    sep = ""
  )
  f <- x@facets
  if (length(f$eq) == 0L) {
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
      nrow(x@vertices) <= 8L
  ) {
    cat("  vertices:\n")
    print(x@vertices)
  }
  invisible(x)
}
