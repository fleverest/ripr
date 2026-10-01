# --- The rcdd bridge ----------------------------------------------------------
#
# What `rcdd` is (according to CRAN):
# R interface to (some of) cddlib (<https://github.com/cddlib/cddlib>). Converts
# back and forth between two representations of a convex polytope: as solution
# of a set of linear equalities and inequalities and as convex hull of set of
# points and rays. Also does linear programming and redundant generator
# elimination (for example, convex hull in n dimensions). All functions can use
# exact infinite-precision rational arithmetic.
#
# Every call into `rcdd` (cddlib: H/V double description and LP in exact
# rationals), and every matrix carrying its flag columns, sits here; the rest
# of the package sees plain lists with flags stripped and signs undone.

# --- Elementwise conversions --------------------------------------------------

#' Doubles to rcdd's rational form
#'
#' Exact for the double given: `d2q(1/3)` is the double nearest a third, not a
#' third.
#' @keywords internal
#' @noRd
as_qmatrix <- function(x) {
  storage.mode(x) <- "double"
  if (!all(is.finite(x))) {
    # rcdd's own error here is uninformative.
    stop("`x` must be finite to have a rational form.", call. = FALSE)
  }
  out <- rcdd::d2q(x)
  dim(out) <- dim(x)
  dimnames(out) <- dimnames(x)
  out
}


#' rcdd's rational form back to doubles, shaped like `x`
#' @keywords internal
#' @noRd
from_qmatrix <- function(x) {
  out <- if (is.character(x)) rcdd::q2d(x) else as.numeric(x)
  dim(out) <- dim(x)
  dimnames(out) <- dimnames(x)
  out
}


# --- Flag-column matrices in and out ------------------------------------------
#
# H-row `(l, b, -a)` means `a . x <= b` when `l == 0` and `a . x == b` when
# `l == 1`. `cddlib` itself stores `b - a.x >= 0`.
#
# V-row `(l, t, x)`: `t == 1` is a point, `t == 0` is a ray, and `l == 1` marks
# a line (lineality, a ray whose negation is also in the set).

#' A H-representation list to an rcdd H-representation qmatrix
#'
#' No constraints become the trivially true row `0 . x <= 1`, as cddlib itself
#' emits, since cddlib needs a row.
#' @keywords internal
#' @noRd
as_hmatrix <- function(h) {
  d <- ncol(h$a)
  if (nrow(h$a) == 0L) {
    return(rcdd::makeH(
      a1 = as_qmatrix(matrix(0, 1L, d)),
      b1 = as_qmatrix(1)
    ))
  }
  a <- as_qmatrix(h$a)
  b <- as_qmatrix(h$b)
  eq <- as.logical(h$eq)
  x <- NULL
  if (any(!eq)) {
    x <- rcdd::makeH(
      a1 = a[!eq, , drop = FALSE],
      b1 = b[!eq],
      x = x
    )
  }
  if (any(eq)) {
    x <- rcdd::makeH(
      a2 = a[eq, , drop = FALSE],
      b2 = b[eq],
      x = x
    )
  }
  x
}


#' An rcdd H-representation to a list
#' @keywords internal
#' @noRd
from_hmatrix <- function(m) {
  list(
    # The coordinate block is `-a`, so undo the sign here and nowhere else.
    a = -from_qmatrix(m[, -(1:2), drop = FALSE]),
    b = as.vector(from_qmatrix(m[, 2L, drop = FALSE])),
    eq = as.vector(from_qmatrix(m[, 1L, drop = FALSE])) == 1
  )
}


#' A package-native V list (one generator per row) to an rcdd V-representation
#' @keywords internal
#' @noRd
as_vmatrix <- function(v) {
  x <- NULL
  if (nrow(v$v) > 0L) {
    x <- rcdd::makeV(points = as_qmatrix(v$v), x = x)
  }
  if (nrow(v$r) > 0L) {
    x <- rcdd::makeV(rays = as_qmatrix(v$r), x = x)
  }
  if (nrow(v$l) > 0L) {
    x <- rcdd::makeV(lines = as_qmatrix(v$l), x = x)
  }
  x
}


#' An rcdd V-representation to a V list
#' @keywords internal
#' @noRd
from_vmatrix <- function(m) {
  linear <- as.vector(from_qmatrix(m[, 1L, drop = FALSE])) == 1
  point <- as.vector(from_qmatrix(m[, 2L, drop = FALSE])) == 1
  x <- from_qmatrix(m[, -(1:2), drop = FALSE])
  list(
    v = x[point & !linear, , drop = FALSE],
    r = x[!point & !linear, , drop = FALSE],
    l = x[!point & linear, , drop = FALSE]
  )
}


#' Normalise a V list to hold at least one point
#'
#' cddlib omits the point block for a cone (`cone(r) + span(l)`, which holds
#' the origin); add the origin, since a chart needs at least one vertex.
#' @keywords internal
#' @noRd
with_origin_vertex <- function(v) {
  if (nrow(v$v) > 0L) {
    return(v)
  }
  if (nrow(v$r) == 0L && nrow(v$l) == 0L) {
    stop(
      "a representation with no generators is the empty set, not a cone.",
      call. = FALSE
    )
  }
  v$v <- matrix(0, nrow = 1L, ncol = ncol(v$v))
  v
}


# --- The rational layer -------------------------------------------------------
#
# `q_` values are flag-column matrices in GMP rationals. Chains of algebra
# (intersection, complement, triangulation) stay here and convert to doubles
# once, at the end.

#' Evaluate an rcdd call without disturbing the global RNG stream
#'
#' cddlib draws from R's RNG in `scdd()` and `lpcdd()`, though its results are
#' deterministic, so restore the seed afterwards.
#' @keywords internal
#' @noRd
without_rng <- function(expr) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", seed, envir = globalenv()))
  }
  expr
}


#' Convert a rational representation to the other one
#'
#' The exact double-description step, in rationals.
#' @keywords internal
#' @noRd
q_scdd <- function(m) {
  without_rng(rcdd::scdd(m, representation = q_kind(m))$output)
}


#' Which representation an rcdd matrix carries
#' @keywords internal
#' @noRd
q_kind <- function(m) {
  kind <- attr(m, "representation")
  if (is.null(kind)) {
    stop(
      "not an rcdd representation: the `representation` attribute is missing.",
      call. = FALSE
    )
  }
  kind
}


#' Stack two representations of the same kind
#'
#' For H this is intersection; for V it is the convex hull, *not* the union.
#' @keywords internal
#' @noRd
q_rbind <- function(a, b) {
  if (!identical(q_kind(a), q_kind(b))) {
    stop(
      "cannot stack an ",
      q_kind(a),
      "-representation onto a ",
      q_kind(b),
      "-representation.",
      call. = FALSE
    )
  }
  if (ncol(a) != ncol(b)) {
    stop(
      "representations have different ambient dimensions: ",
      ncol(a) - 2L,
      " and ",
      ncol(b) - 2L,
      ".",
      call. = FALSE
    )
  }
  out <- rbind(a, b)
  attr(out, "representation") <- q_kind(a)
  out
}


#' Drop rows a representation does not need
#'
#' Redundant inequalities (H) or non-extreme generators (V).
#' @keywords internal
#' @noRd
q_nonredundant <- function(m) {
  # `rcdd::redundant()` errors on a single row.
  if (nrow(m) < 2L) {
    return(m)
  }
  without_rng(rcdd::redundant(m, representation = q_kind(m))$output)
}


#' Exact feasibility of a rational H-representation
#'
#' An exact zero-objective LP.
#' @keywords internal
#' @noRd
q_is_empty <- function(m) {
  d <- ncol(m) - 2L
  lp <- without_rng(rcdd::lpcdd(
    m,
    objgrd = as_qmatrix(rep(0, d)),
    objcon = "0",
    minimize = TRUE
  ))
  status <- lp$solution.type
  switch(
    status,
    "Optimal" = FALSE,
    "Inconsistent" = TRUE,
    "StrucInconsistent" = TRUE,
    stop(
      "`lpcdd()` returned the status \"",
      status,
      "\", which a feasibility program with a zero objective should not ",
      "produce. This is a bug in ripr.",
      call. = FALSE
    )
  )
}


#' Does every row of an H-representation hold on the region `qv` generates?
#'
#' Exact: `a . x <= b` holds everywhere when it holds at every vertex, with
#' `a . r <= 0` on every ray and `a . l == 0` on every line; an equality row
#' must be tight at every generator.
#' @keywords internal
#' @noRd
q_holds <- function(qv, rows) {
  if (nrow(qv) == 0L || nrow(rows) == 0L) {
    return(TRUE)
  }
  # Entry (i, j) is `t_i b_j - a_j . x_i`, with `t_i` 1 for a vertex, 0 else.
  s <- rcdd::qsign(rcdd::qmatmult(
    qv[, -1L, drop = FALSE],
    t(rows[, -1L, drop = FALSE])
  ))
  tight <- outer(qv[, 1L] == "1", rows[, 1L] == "1", `|`)
  all(s >= 0L) && all(s[tight] == 0L)
}


#' Subset rows of a representation, keeping its kind
#'
#' Plain `[` drops the `representation` attribute `q_rbind()` and `q_scdd()`
#' dispatch on.
#' @keywords internal
#' @noRd
q_subrows <- function(m, rows) {
  out <- m[rows, , drop = FALSE]
  attr(out, "representation") <- q_kind(m)
  out
}


#' Reverse an inequality row: `a . x <= b` becomes `a . x >= b`
#'
#' Closed (`>=`, not `>`), so cells overlap only on their boundaries.
#' @keywords internal
#' @noRd
q_reverse_ineq <- function(m, rows) {
  m[rows, -1L] <- rcdd::qneg(m[rows, -1L])
  m
}


#' The exact dimension of a rational representation
#'
#' The dimension of the affine hull, which is the ambient dimension minus the
#' number of independent equalities.
#' @keywords internal
#' @noRd
q_dim <- function(m) {
  h <- if (identical(q_kind(m), "H")) q_nonredundant(m) else q_scdd(m)
  (ncol(h) - 2L) - sum(h[, 1L] == "1")
}


#' Feasibility of a package-native H list
#'
#' [is_empty()] for constraints that are not yet a region.
#' @keywords internal
#' @noRd
h_is_empty <- function(h) q_is_empty(as_hmatrix(h))


#' The facets of a generator set, and which vertices lie on each
#'
#' Returns `h`, the H-representation, and `on`, per row of `h` the rows of `v`
#' lying on it (free from cddlib's incidence output).
#' @keywords internal
#' @noRd
q_facets <- function(v) {
  out <- without_rng(rcdd::scdd(v, representation = "V", incidence = TRUE))
  list(h = out$output, on = out$incidence)
}


# --- The dual-representation record -------------------------------------------

#' All internal representations of a convex polyhedron
#'
#' H and V, each in double (`h` = `(a, b, eq)`, `v` = `list(v, r, l)`) and
#' rational (`qh`, `qv`) form. Build with `hv_fill()`.
#' @keywords internal
#' @noRd
hv_region <- new_class(
  "hv_region",
  properties = list(
    h = class_list,
    v = class_list,
    qh = class_any,
    qv = class_any
  )
)


#' An `hv_region` from whichever representation is at hand
#'
#' Give one of `qh`, `qv`, `v` or `h`; the rest follows by one exact dd step.
#' `qh` is made non-redundant first. `v` alongside `h` overrides only the
#' double generators (for display); the rational V is still derived from `h`.
#' @keywords internal
#' @noRd
hv_fill <- function(qh = NULL, qv = NULL, h = NULL, v = NULL) {
  if (!is.null(h)) {
    qh <- q_nonredundant(as_hmatrix(h))
  } else if (!is.null(qh)) {
    qh <- q_nonredundant(qh)
  }
  if (!is.null(qh)) {
    qv <- q_scdd(qh)
  } else {
    if (is.null(qv)) {
      qv <- as_vmatrix(v)
    }
    qh <- q_scdd(qv)
  }
  hv_region(
    h = if (is.null(h)) from_hmatrix(qh) else h,
    v = if (is.null(v)) with_origin_vertex(from_vmatrix(qv)) else v,
    qh = qh,
    qv = qv
  )
}


#' An empty generator block of the right ambient dimension
#' @keywords internal
#' @noRd
no_generators <- function(d) matrix(numeric(0), nrow = 0L, ncol = d)
