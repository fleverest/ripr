# Tests for the rcdd bridge in R/hrep.R.

# Columns are only determined up to order, and vertex sets up to permutation.
sorted_cols <- function(m) {
  m[, order(apply(m, 2L, paste, collapse = ",")), drop = FALSE]
}

# Two orthonormal bases span the same subspace iff each projects the other onto
# itself, so compare projectors rather than the arbitrary basis choice.
projector <- function(b) {
  if (ncol(b) == 0L) matrix(0, nrow(b), nrow(b)) else tcrossprod(qr.Q(qr(b)))
}

plurality_cell <- function() {
  simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1)))
}

h_to_v <- function(h) {
  from_vmatrix(q_scdd(as_hmatrix(h)))
}


# --- Round trips --------------------------------------------------------------

test_that("H and V representations round trip", {
  for (s in list(
    plurality_cell(),
    simplex_region(vertices = diag(3)),
    polytope_region(vertices = cbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)))
  )) {
    back <- h_to_v(h_rep(s))
    expect_equal(ncol(back$r), 0L)
    expect_equal(ncol(back$l), 0L)
    expect_identical(
      sorted_cols(back$v),
      sorted_cols(s@vertices)
    )
  }
})


# --- The sign convention ------------------------------------------------------

test_that("h_rep of a simplex holds its own vertices and excludes outside points", {
  s <- simplex_region(vertices = diag(3))
  h <- h_rep(s)

  for (j in seq_len(ncol(s@vertices))) {
    expect_true(all(h$a %*% s@vertices[, j] <= h$b + rounding_tol(1)))
  }
  # Reflecting the region through the origin would admit this one.
  expect_false(all(h$a %*% c(2, -1, 0) <= h$b + rounding_tol(1)))
})


# --- The equality row ---------------------------------------------------------

test_that("a cell in the standard simplex carries exactly one equality row", {
  h <- h_rep(simplex_region(vertices = diag(3)))

  expect_equal(sum(h$eq), 1L)
  row <- which(h$eq)
  # cddlib is free to return the row scaled by any non-zero constant, and does
  # return `-sum(theta) == -1` here, so normalise by a coefficient.
  a <- h$a[row, ]
  expect_equal(a / a[1L], c(1, 1, 1))
  expect_equal(h$b[row] / a[1L], 1)
})


# --- Cross-check against the hand-written generators --------------------------

test_that("halfspace_region's own generators agree with cddlib's", {
  # `anchor`, `-unit` and `basis` are computed by the constructor without any
  # reference to rcdd, so this compares two independently derived answers and
  # is the cheapest available check that the sign convention is right.
  s <- halfspace_region(normal = c(1, -1, 0), offset = 2)
  ours <- v_rep(s)
  theirs <- h_to_v(h_rep(s))

  # The point may be any point of the halfspace, so check membership, not
  # equality.
  expect_equal(ncol(theirs$v), 1L)
  expect_true(contains(s, ours$v[, 1L]))
  expect_true(contains(s, theirs$v[, 1L]))

  expect_equal(projector(theirs$l), projector(ours$l), tolerance = rounding_tol(1))

  # A ray is only determined modulo the lineality space, so compare the two
  # after projecting the lineality directions out.
  off <- diag(3) - projector(ours$l)
  ray_ours <- off %*% ours$r[, 1L]
  ray_theirs <- off %*% theirs$r[, 1L]
  expect_equal(
    ray_theirs / sqrt(sum(ray_theirs^2)),
    ray_ours / sqrt(sum(ray_ours^2)),
    tolerance = rounding_tol(1)
  )
})


# --- Exactness ----------------------------------------------------------------

test_that("conversion to rationals is exact for the double", {
  # Nothing rounds, so this is an identity, not an approximation. `1/3` becomes
  # the rational the double actually is -- long, and not one third.
  x <- c(0.5, 1 / 3, 0, -0.25, 2, pi, 55 / 10000024)
  expect_identical(rcdd::q2d(ripr:::as_qmatrix(x)), x)
  expect_equal(ripr:::as_qmatrix(0.5), "1/2")
})


test_that("as_qmatrix keeps the shape of its input", {
  m <- matrix(c(0.5, 1 / 3, 0.25, 0.75), 2L, 2L)
  expect_identical(rcdd::q2d(ripr:::as_qmatrix(m)), m)
})


test_that("the same double always gives the same rational", {
  # This is the facet-sharing property the package relies on: two cells handed
  # the same vertex are seen to share it. It holds because the conversion is a
  # function of the double and nothing rounds.
  expect_identical(ripr:::as_qmatrix(1 / 3), ripr:::as_qmatrix((1 / 6) * 2))

  # And the honest limit of it: one ulp apart is a different point, and is
  # treated as one. A consumer needing more than this needs a tolerance-based
  # facet matcher, not a different encoding here.
  expect_false(identical(1 / 3, 1 - 2 / 3))
  expect_false(ripr:::as_qmatrix(1 / 3) == ripr:::as_qmatrix(1 - 2 / 3))
})


test_that("a double H-representation is a lossy intermediate for derived facets", {
  # The limit of the double-valued API
  set.seed(1)
  v <- matrix(stats::runif(30), 3L, 10L)
  h <- ripr:::v_to_h(list(
    v = v,
    r = ripr:::no_generators(3),
    l = ripr:::no_generators(3)
  ))
  back <- h_to_v(h)

  # The true answer, which staying in rationals throughout would have given.
  expect_equal(nrow(unique(round(t(back$v), 9L))), 9L)
  # What the double round trip actually returns.
  expect_gt(ncol(back$v), 9L)

  # It is a perturbation, not a wrong answer: every reported vertex is on the
  # region, and the region still holds the points that generated it.
  expect_lt(max(h$a %*% back$v - h$b), rounding_tol(1))
  expect_lt(max(h$a %*% v - h$b), rounding_tol(1))
})


test_that("a value with no rational form is refused by name", {
  expect_error(ripr:::as_qmatrix(NA_real_), "finite")
})


# --- Predicates ---------------------------------------------------------------

test_that("is_empty is FALSE for regions that hold something", {
  expect_false(is_empty(plurality_cell()))
  expect_false(is_empty(simplex_region(vertices = diag(3))))
  expect_false(is_empty(point_region(theta = c(0.5, 0.3, 0.2))))
  expect_false(is_empty(unconstrained_region(3L)))
  expect_false(is_empty(halfspace_region(normal = c(1, -1, 0))))
})


test_that("is_empty is TRUE for contradictory constraints", {
  # `{theta_1 <= 0}` and `{theta_1 >= 0.5}` inside the standard simplex.
  # Built by concatenating H-rows by hand: there is no `intersect()` yet
  # at the time of writing.
  h <- list(
    a = rbind(
      c(1, 0, 0),
      c(-1, 0, 0),
      c(-1, 0, 0),
      c(0, -1, 0),
      c(0, 0, -1),
      c(1, 1, 1)
    ),
    b = c(0, -0.5, 0, 0, 0, 1),
    eq = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
  )
  expect_true(ripr:::h_is_empty(h))

  # Drop the contradictory row and it is the sub-simplex `{theta_1 <= 0}`.
  keep <- -2L
  expect_false(ripr:::h_is_empty(
    list(a = h$a[keep, ], b = h$b[keep], eq = h$eq[keep])
  ))
})


test_that("is_bounded distinguishes the bounded regions from the unbounded", {
  expect_true(is_bounded(plurality_cell()))

  expect_false(is_bounded(halfspace_region(normal = c(1, -1, 0))))
  expect_false(is_bounded(unconstrained_region(3L)))
})


test_that("both predicates reduce over a union with all()", {
  s <- plurality_cell()
  h <- halfspace_region(normal = c(1, -1, 0))

  expect_true(is_bounded(union_region(s, s)))
  expect_false(is_bounded(union_region(s, h)))
  expect_false(is_empty(union_region(s, h)))
})


# --- Refusals -----------------------------------------------------------------

test_that("a union has neither representation, and says so", {
  u <- union_region(
    plurality_cell(),
    simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
  )
  expect_error(h_rep(u), "union_region")
  expect_error(v_rep(u), "union_region")
})


# --- Containment ---------------------------------------------------------------

test_that("no function outside R/hrep.R touches rcdd", {
  root <- testthat::test_path("..", "..")
  skip_if(!dir.exists(file.path(root, "R")), "not running from the source tree")

  files <- list.files(
    file.path(root, "R"),
    pattern = "[.]R$",
    full.names = TRUE
  )
  offenders <- Filter(
    function(f) {
      basename(f) != "polyhedra.R" &&
        any(grepl("rcdd::", readLines(f, warn = FALSE), fixed = TRUE))
    },
    files
  )
  expect_equal(basename(offenders), character(0))
})


# --- RNG neutrality -----------------------------------------------------------

test_that("conversions leave the global RNG stream alone", {
  # cddlib randomises internally through R's own RNG, so a bare `scdd()` call
  # advances `.Random.seed`. The bridge saves and restores it: geometry is
  # exact and deterministic, and "same seed, same fit" must not depend on how
  # many conversions a region's construction happened to run.
  set.seed(42)
  before <- .Random.seed
  s <- simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1)))
  h_rep(s)
  v_rep(s)
  is_empty(s)
  invisible(q_scdd(q_vrep(s)))
  invisible(q_nonredundant(q_vrep(s)))
  expect_identical(.Random.seed, before)
})
