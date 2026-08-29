#' @include region.R
#' @include polyhedra.R
NULL

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
#' [point_region()] and [real_region()] are all special cases that
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
#'   don't start from a rounding.}
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
    if (!is_bounded(x) && !S7_inherits(x, real_region)) ", unbounded",
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
