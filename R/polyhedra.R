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
#
# Every call into `rcdd`, and every matrix carrying rcdd's flag columns,
# sits here. Regions describe themselves by generators (vertices, or
# a normal and an offset) but complement, intersection, triangulation and
# emptiness pruning are all defined facet-by-facet, so they need the dual
# description. `h_rep()` and `v_rep()` convert between the dual representations;
# everything else in the package sees package-native lists with the flags
# already stripped and the sign convention already undone.
#
# The conversion runs in GMP rationals rather than doubles because both
# consumers need it to. A complement decomposition has to see two cells share a
# facet exactly, or it fails to prune the sliver between them and may produce
# cells of width 1e-17 that survive into the downstream loop.

# --- Elementwise conversions --------------------------------------------------

#' Doubles to rcdd's rational form
#'
#' [rcdd::d2q()] is exact for the double it is given: `d2q(0.5)` is `"1/2"`, and
#' `d2q(1/3)` is `"6004799503160661/18014398509481984"`, which is not one third
#' but is precisely the double that `1/3` evaluates to.
#'
#' No rounding here, so a vertex that goes in comes out describing the same
#' point, and two cells given the same double always get the same rational and
#' so genuinely share a facet.
#'
#' @param x A numeric vector or matrix.
#' @return A character vector or matrix of rationals, shaped like `x`.
#' @keywords internal
#' @noRd
as_qmatrix <- function(x) {
  storage.mode(x) <- "double"
  if (!all(is.finite(x))) {
    # rcdd's own message for this names no argument and no value, and every
    # path here is deep inside a conversion the caller did not ask for.
    stop("`x` must be finite to have a rational form.", call. = FALSE)
  }
  out <- rcdd::d2q(x)
  dim(out) <- dim(x)
  dimnames(out) <- dimnames(x)
  out
}


#' rcdd's rational form back to doubles
#' @param x A character vector or matrix of rationals.
#' @return A numeric vector or matrix, shaped like `x`.
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
# `l == 1`. Note the sign: cddlib stores `b - a.x >= 0`, so the coordinate block
# is `-a`. Getting this wrong reflects the region through the origin, which
# passes most smoke tests.
#
# V-row `(l, t, x)`: `t == 1` is a point, `t == 0` is a ray, and `l == 1` marks
# a line (lineality, a ray whose negation is also in the set).

#' A H-representation list to an rcdd H-representation qmatrix
#'
#' A region with no constraints has no rows to hand cddlib, so it is written as
#' the trivially true row `0 . x <= 1`. That is what cddlib itself emits for the
#' same set, and it keeps every caller off a zero-row special case.
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
  d <- ncol(m) - 2L
  list(
    # The coordinate block is `-a`, so undo the sign here and nowhere else.
    a = -from_qmatrix(m[, -(1:2), drop = FALSE]),
    b = as.vector(from_qmatrix(m[, 2L, drop = FALSE])),
    eq = as.vector(from_qmatrix(m[, 1L, drop = FALSE])) == 1
  )
}


#' A package-native V list to an rcdd V-representation
#'
#' Generators arrive one per column, as vertices do everywhere else in the
#' package; cddlib wants one per row.
#' @keywords internal
#' @noRd
as_vmatrix <- function(v) {
  x <- NULL
  if (ncol(v$v) > 0L) {
    x <- rcdd::makeV(points = t(as_qmatrix(v$v)), x = x)
  }
  if (ncol(v$r) > 0L) {
    x <- rcdd::makeV(rays = t(as_qmatrix(v$r)), x = x)
  }
  if (ncol(v$l) > 0L) {
    x <- rcdd::makeV(lines = t(as_qmatrix(v$l)), x = x)
  }
  x
}


#' An rcdd V-representation to a V list
#' @keywords internal
#' @noRd
from_vmatrix <- function(m) {
  d <- ncol(m) - 2L
  linear <- as.vector(from_qmatrix(m[, 1L, drop = FALSE])) == 1
  point <- as.vector(from_qmatrix(m[, 2L, drop = FALSE])) == 1
  x <- t(from_qmatrix(m[, -(1:2), drop = FALSE]))
  dim(x) <- c(d, nrow(m))
  list(
    v = x[, point & !linear, drop = FALSE],
    r = x[, !point & !linear, drop = FALSE],
    l = x[, !point & linear, drop = FALSE]
  )
}


#' Normalise a V list to hold at least one point
#'
#' cddlib omits the point block entirely when the polyhedron is a cone: the
#' halfspace `{t1 <= 0}` in `R^3` comes back as `nv = 0, nr = 1, nl = 2`. A
#' V-representation with no points denotes `cone(r) + span(l)`, which always
#' contains the origin. We just add the point, since a chart's barycentric
#' coordinates on zero vertices is undefined, and `v_rep(halfspace_region)`
#' already yields the anchored form, so this way everything is consistent.
#'
#' @keywords internal
#' @noRd
with_origin_vertex <- function(v) {
  if (ncol(v$v) > 0L) {
    return(v)
  }
  if (ncol(v$r) == 0L && ncol(v$l) == 0L) {
    stop(
      "a representation with no generators is the empty set, not a cone.",
      call. = FALSE
    )
  }
  v$v <- matrix(0, nrow = nrow(v$v), ncol = 1L)
  v
}


# --- The rational layer -------------------------------------------------------
#
# `h_rep()` and `v_rep()` below return doubles, and doubles are where the
# boundary between this file and the rest of the package sits: charts,
# projection, optimisation and the Bernstein enclosure are all floating point
# and always will be. But a facet derived from vertices in general position is
# an exact rational with a numerator far too long for a double, so a *chain* of
# operations that goes back through doubles at every step rounds away the
# exactness each individual step had. Measured: ten points in general position
# in `R^3` have a hull with 9 vertices, and a chain that round-trips through the
# double form returns 24.
#
# So the rule is to cross the boundary once per chain, not once per operation.
# These functions are the inside of that chain. They speak rcdd's flag-column
# matrices, they never round, and the algebra of a later phase -- intersection,
# complement, triangulation -- should compose here and convert to doubles only
# when it hands a finished cell back.
#
# `q_` marks that a value is a flag-column matrix in GMP rationals. Nothing
# carrying that prefix should leave this file except into the geometry algebra.

#' Evaluate an rcdd call without disturbing the global RNG stream
#'
#' cddlib randomises internally through R's own RNG for `scdd()` and
#' `lpcdd()`. The randomness is a runtime hedge inherited from cddlib, rerouted
#' through R's RNG by rcdd only to abide by CRAN policy (no direct `rand()` or
#' `srand()`), but the result itself is deterministic, so we just reset the
#' PRNG state after such calls.
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
#' The exact double-description step, staying in rationals throughout. Reads
#' which representation it was handed from the matrix rather than being told.
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
#' For an H-representation this is intersection: a point satisfies the stack
#' exactly when it satisfies both. For a V-representation it is *not* union:
#' it is the convex hull of the two generator sets, so region algebra needs
#' to do unions structurally rather than by stacking.
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
#' Redundant inequalities for an H-representation, non-extreme generators for a
#' V-representation. Intersection produces both routinely, so we will trim them
#' with this.
#' @keywords internal
#' @noRd
q_nonredundant <- function(m) {
  without_rng(rcdd::redundant(m, representation = q_kind(m))$output)
}


#' Exact feasibility of a rational H-representation
#'
#' A linear program with a zero objective: the region is empty exactly when the
#' constraints have no common solution. Runs in GMP rationals, so the answer is
#' a decision rather than an estimate.
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
  # cddlib reports feasibility through the status string, and the two
  # inconsistent statuses are the answer we came for rather than a failure.
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


#' Feasibility of a package-native H list
#'
#' [is_empty()] for constraints that are not yet a region. Intersection produces
#' H lists before it produces regions, and they still need pruning.
#' @param h A list with `a`, `b` and `eq`, as `h_rep()` returns.
#' @return `TRUE` or `FALSE`.
#' @keywords internal
#' @noRd
h_is_empty <- function(h) q_is_empty(as_hmatrix(h))


#' The facets of a generator set, and which vertices lie on each
#'
#' Triangulation cones from an interior point to each facet, so it needs each
#' facet's vertex set.
#'
#' cddlib solves this during the dd run and will hand it over exactly at no
#' extra cost.
#'
#' @param v An rcdd V-representation.
#' @return A list with `h`, the H-representation, and `on`, a list with one
#'   integer vector per row of `h` giving the rows of `v` lying on it.
#' @keywords internal
#' @noRd
q_facets <- function(v) {
  out <- without_rng(rcdd::scdd(v, representation = "V", incidence = TRUE))
  list(h = out$output, on = out$incidence)
}


#' The exact centroid of a representation's points
#'
#' The point that triangulation cones from. The centroid of a vertex set lies
#' in the relative interior of its hull, which is what a cone decomposition
#' uses.
#'
#' @param v An rcdd V-representation.
#' @return A length-`d` character vector of rationals.
#' @keywords internal
#' @noRd
q_centroid <- function(v) {
  point <- v[, 2L] == "1"
  if (!any(point)) {
    stop("a representation with no points has no centroid.", call. = FALSE)
  }
  x <- v[point, -(1:2), drop = FALSE]
  n <- as.character(nrow(x))
  vapply(
    seq_len(ncol(x)),
    function(j) rcdd::qdq(rcdd::qsum(x[, j]), n),
    character(1)
  )
}


# --- The double description itself --------------------------------------------

#' Convert a V list to a package-native H list
#' @keywords internal
#' @noRd
v_to_h <- function(v) {
  from_hmatrix(q_scdd(as_vmatrix(v)))
}


#' Convert a H list to a package-native V list
#' @keywords internal
#' @noRd
h_to_v <- function(h) {
  from_vmatrix(q_scdd(as_hmatrix(h)))
}


#' An empty generator block of the right ambient dimension
#' @keywords internal
#' @noRd
no_generators <- function(d) matrix(numeric(0), nrow = d, ncol = 0L)


# --- Redundancy ---------------------------------------------------------------

#' Drop vertices that are not extreme
#'
#' A thin wrapper on [rcdd::redundant()] for the vertex case: given one vertex
#' per column, returns the subset that are actually vertices of the hull. An
#' intersection produces redundant generators routinely, and carrying them
#' costs every later conversion.
#'
#' @param v A `d` by `n` numeric matrix, one vertex per column.
#' @return A `d` by `k` numeric matrix, a column subset of `v`.
#' @keywords internal
#' @noRd
redundant_vertices <- function(v) {
  d <- nrow(v)
  from_vmatrix(q_nonredundant(as_vmatrix(
    list(v = v, r = no_generators(d), l = no_generators(d))
  )))$v
}
