# Tests for the polyhedron base class in R/region.R.
#
# Every implemented convex region is a polyhedron (at the time of writing), and
# the base chart is written once against the generator triple. These tests
# restate the old per-class formulas inline and check there is zero difference.

# Draw random local coords for comparison.
random_u <- function(n_par, n = 50L) {
  matrix(stats::rnorm(n * n_par, sd = 2), nrow = n_par, ncol = n)
}

# --- The chart identity -------------------------------------------------------

test_that("the base chart reproduces polytope_region's chart exactly", {
  square <- polytope_region(
    vertices = cbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  )
  vertices <- square@vertices
  old_to_theta <- function(u) as.vector(vertices %*% softmax0(u))
  old_jacobian <- function(u) vertices %*% softmax0_jacobian(softmax0(u))

  ch <- chart(square)
  expect_identical(ch$n_par, 3L)
  set.seed(1)
  u_mat <- random_u(ch$n_par)
  for (i in seq_len(ncol(u_mat))) {
    expect_identical(ch$to_theta(u_mat[, i]), old_to_theta(u_mat[, i]))
    expect_identical(ch$jacobian(u_mat[, i]), old_jacobian(u_mat[, i]))
  }
  expect_identical(
    ch$to_theta_batch(u_mat),
    vertices %*%
      vapply(seq_len(ncol(u_mat)), \(i) softmax0(u_mat[, i]), numeric(4))
  )
})


test_that("the base chart reproduces simplex_region's chart exactly", {
  plurality <- simplex_region(
    vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))
  )
  vertices <- plurality@vertices

  ch <- chart(plurality)
  expect_identical(ch$n_par, 2L)
  set.seed(2)
  u_mat <- random_u(ch$n_par)
  for (i in seq_len(ncol(u_mat))) {
    u <- u_mat[, i]
    expect_identical(ch$to_theta(u), as.vector(vertices %*% softmax0(u)))
    expect_identical(
      ch$jacobian(u),
      vertices %*% softmax0_jacobian(softmax0(u))
    )
  }
})


test_that("the base chart reproduces halfspace_region's chart exactly", {
  for (space in list(
    halfspace_region(normal = c(1, -1), offset = 0),
    halfspace_region(normal = c(1, -1, 0), offset = 2),
    halfspace_region(normal = c(2, 1, -1), offset = -1)
  )) {
    d <- length(space@normal)
    anchor <- space@anchor
    basis <- space@basis
    unit <- space@unit
    old_to_theta <- function(u) {
      anchor + as.vector(basis %*% u[-d]) - softplus(u[d]) * unit
    }
    old_jacobian <- function(u) cbind(basis, -sigmoid(u[d]) * unit)

    ch <- chart(space)
    expect_identical(ch$n_par, d)
    set.seed(3)
    u_mat <- random_u(ch$n_par)
    for (i in seq_len(ncol(u_mat))) {
      expect_identical(ch$to_theta(u_mat[, i]), old_to_theta(u_mat[, i]))
      expect_identical(ch$jacobian(u_mat[, i]), old_jacobian(u_mat[, i]))
    }
  }
})


test_that("the base chart reproduces point_region's chart exactly", {
  theta <- c(0.5, 0.3, 0.2)
  space <- point_region(theta = theta)

  ch <- chart(space)
  expect_identical(ch$n_par, 0L)
  expect_identical(ch$to_theta(numeric(0)), theta)
  expect_identical(ch$jacobian(numeric(0)), matrix(0, 3L, 0L))
  expect_identical(ch$from_theta(theta), numeric(0))
  expect_identical(dim(ch$seed(5L)), c(0L, 5L))
})


test_that("the base chart reproduces unconstrained_region's chart exactly", {
  space <- unconstrained_region(3L)

  ch <- chart(space)
  expect_identical(ch$n_par, 3L)
  set.seed(4)
  u_mat <- random_u(ch$n_par)
  for (i in seq_len(ncol(u_mat))) {
    expect_identical(ch$to_theta(u_mat[, i]), u_mat[, i])
    expect_identical(ch$jacobian(u_mat[, i]), diag(3L))
  }
  expect_identical(ch$to_theta_batch(u_mat), u_mat)
})


# --- Generators and facets ----------------------------------------------------

test_that("every class carries the generator triple it used to hand-write", {
  # We store both the original construction vertices and the polyhedron
  # representations
  square <- cbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  g <- polytope_region(vertices = square)@generators
  expect_identical(g$v, square)
  expect_identical(dim(g$r), c(2L, 0L))
  expect_identical(dim(g$l), c(2L, 0L))

  h <- halfspace_region(normal = c(1, -1, 0), offset = 2)
  g <- h@generators
  expect_identical(g$v, matrix(h@anchor, ncol = 1L))
  expect_identical(g$r, matrix(-h@unit, ncol = 1L))
  expect_identical(g$l, h@basis)

  g <- point_region(theta = c(0.5, 0.5))@generators
  expect_identical(g$v, matrix(c(0.5, 0.5), ncol = 1L))

  g <- unconstrained_region(2L)@generators
  expect_identical(g$v, matrix(0, 2L, 1L))
  expect_identical(g$l, diag(2L) + 0)
})


test_that("declared facets are stored exactly; derived facets describe the hull", {
  h <- halfspace_region(normal = c(1, -1, 0), offset = 2)
  expect_identical(
    h@facets,
    list(a = matrix(c(1, -1, 0), nrow = 1L), b = 2, eq = FALSE)
  )

  p <- point_region(theta = c(0.5, 0.3, 0.2))
  expect_identical(p@facets$eq, rep(TRUE, 3L))
  expect_identical(p@facets$b, c(0.5, 0.3, 0.2))

  expect_identical(nrow(unconstrained_region(3L)@facets$a), 0L)

  # A polytope derives its facets once at construction: the unit square has
  # four, and they cut out exactly the square.
  square <- polytope_region(
    vertices = cbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  )
  expect_identical(nrow(square@facets$a), 4L)
  expect_true(contains(square, c(0.5, 0.5)))
  expect_true(contains(square, c(0, 0)))
  expect_false(contains(square, c(1.1, 0.5)))
  expect_false(contains(square, c(0.5, -0.1)))
})


# --- The new constructors -----------------------------------------------------

test_that("polyhedron_region() builds from generators", {
  # The halfspace `{theta_1 <= 0}` in R^2, by hand.
  space <- polyhedron_region(
    vertices = matrix(c(0, 0), ncol = 1L),
    rays = matrix(c(-1, 0), ncol = 1L),
    lines = matrix(c(0, 1), ncol = 1L)
  )
  expect_true(S7_inherits(space, polyhedron_region))
  expect_identical(space_dim(space), 2L)
  expect_false(is_bounded(space))
  expect_true(contains(space, c(-3, 5)))
  expect_false(contains(space, c(1, 0)))

  # A cone with no vertices is anchored at the origin.
  cone <- polyhedron_region(rays = cbind(c(1, 0), c(0, 1)))
  expect_identical(cone@generators$v, matrix(0, 2L, 1L))
  expect_true(contains(cone, c(2, 3)))
  expect_false(contains(cone, c(-1, 1)))

  expect_error(polyhedron_region(), "cannot all be NULL")
  expect_error(
    polyhedron_region(vertices = cbind(c(0, Inf))),
    "finite"
  )
})


test_that("h_region() builds from the half-space form and keeps the rows", {
  square <- h_region(a = rbind(diag(2), -diag(2)), b = c(1, 1, 0, 0))
  expect_true(S7_inherits(square, polyhedron_region))
  expect_true(is_bounded(square))
  v <- v_rep(square)$v
  expect_identical(ncol(v), 4L)
  expect_true(contains(square, c(0.5, 0.5)))
  expect_false(contains(square, c(1.5, 0.5)))
  # The declared rows are the facets, untouched.
  expect_identical(h_rep(square)$b, c(1, 1, 0, 0))

  half <- h_region(a = matrix(c(1, -1), nrow = 1L), b = 0)
  expect_false(is_bounded(half))
  expect_true(contains(half, c(-1, 1)))

  # Infeasible constraints are refused rather than fabricating an anchor.
  expect_error(
    h_region(a = rbind(c(1, 0), c(-1, 0)), b = c(0, -1)),
    "no common solution"
  )
})


test_that("the parent validator runs through every subclass constructor", {
  # S7 skips recursive validation at construction when the parent class is
  # concrete, assuming the parent's constructor already validated its part.
  # That holds only because every subclass constructor delegates to
  # `polyhedron_region()`; built directly on `new_object(S7_object(), ...)`,
  # a logical vertex matrix used to slip through unvalidated.
  expect_error(polytope_region(vertices = diag(2) > 0), "numeric matrix")
  expect_error(simplex_region(vertices = diag(2) > 0), "numeric matrix")
  expect_error(polyhedron_region(vertices = diag(2) > 0), "numeric matrix")

  # Integer matrices are numeric already and flow through as they are.
  square <- polytope_region(
    vertices = matrix(c(0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L), nrow = 2L)
  )
  expect_identical(nrow(h_rep(square)$a), 4L)
  expect_true(contains(square, c(0.5, 0.5)))
  expect_false(contains(square, c(1.5, 0.5)))
})


# --- from_theta: the constrained recovery -------------------------------------

test_that("from_theta round-trips on polytopes outside the probability simplex", {
  # The old pseudo-inverse recovery never constrained `sum(alpha) = 1`, which
  # is invisible inside the probability simplex and wrong by ~0.27 on a unit
  # square. These are the shapes where the constrained solve earns its keep.
  square <- polytope_region(
    vertices = cbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  )
  hexagon <- polytope_region(
    vertices = vapply(
      0:5,
      \(k) c(cos(k * pi / 3), sin(k * pi / 3)),
      numeric(2)
    )
  )
  set.seed(5)
  for (space in list(square, hexagon)) {
    ch <- chart(space)
    for (i in seq_len(20L)) {
      u <- stats::rnorm(ch$n_par)
      theta <- ch$to_theta(u)
      round_trip <- ch$to_theta(ch$from_theta(theta))
      expect_true(all(is.finite(round_trip)))
      expect_equal(round_trip, theta, tolerance = 1e-9)
    }
  }
})


test_that("from_theta round-trips on an offset halfspace", {
  # Exact for an offset halfspace, where a plain pseudo-inverse over
  # `cbind(V, L, R)` is rank-deficient and lands far away.
  space <- halfspace_region(normal = c(0, 1), offset = 2)
  ch <- chart(space)
  set.seed(6)
  for (i in seq_len(20L)) {
    theta <- ch$to_theta(stats::rnorm(ch$n_par))
    expect_equal(ch$to_theta(ch$from_theta(theta)), theta, tolerance = 1e-9)
  }
})


test_that("from_theta of an outside point recovers the projection", {
  # `generator_weights()` minimises the distance to the region, so the
  # recovered coordinates map to the Euclidean projection: from_theta is now
  # safe to call without composing with project() first.
  square <- polytope_region(
    vertices = cbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  )
  ch <- chart(square)
  outside <- c(2, 0.5)
  expect_equal(
    ch$to_theta(ch$from_theta(outside)),
    c(1, 0.5),
    tolerance = 1e-6
  )
  expect_equal(project(square, outside), c(1, 0.5), tolerance = 1e-6)
})


# --- Structure ----------------------------------------------------------------

test_that("the hierarchy nests as documented", {
  s <- simplex_region(vertices = diag(3))
  p <- point_region(theta = c(1, 0))
  for (space in list(
    s,
    polytope_region(vertices = diag(3)),
    p,
    halfspace_region(normal = c(1, -1)),
    unconstrained_region(2L)
  )) {
    expect_true(S7_inherits(space, polyhedron_region))
    expect_true(S7_inherits(space, convex_region))
  }
  expect_true(S7_inherits(s, polytope_region))
  # A point is the one-vertex polytope.
  expect_true(S7_inherits(p, polytope_region))
})
