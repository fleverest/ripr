# Tests for the polyhedron base class in R/region.R.
#
# Every implemented convex region is a polyhedron (at the time of writing), and
# the base chart is written once against the generator triple. Coordinates are
# the generator weights themselves (barycentric on the vertices, non-negative
# ray coefficients and free lineality coords) so the map is linear and the
# constraints are handed to the optimiser explicitly through `lower` and `heq`
# rather than being reparametrised away (like we did in an earlier version via
# softmax / softplus for stats::optim).

# Check a coordinate matrix against the chart's own constraint declaration.
expect_feasible <- function(ch, u_mat) {
  expect_true(all(
    u_mat >=
      matrix(
        ch$lower,
        nrow(u_mat),
        ncol(u_mat)
      )
  ))
  if (!is.null(ch$heq)) {
    for (i in seq_len(ncol(u_mat))) {
      expect_equal(ch$heq(u_mat[, i]), 0, tolerance = 1e-12)
    }
  }
}

# --- The chart identity -------------------------------------------------------

test_that("the chart is the linear generator map on a polytope", {
  square <- polytope_region(
    vertices = cbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  )
  vertices <- square@vertices

  ch <- chart(square)
  # One coordinate per vertex: the simplex constraint is declared, not
  # substituted out.
  expect_identical(ch$n_par, 4L)
  expect_identical(ch$lower, rep(0, 4L))
  set.seed(1)
  u_mat <- ch$seed(50L)
  expect_feasible(ch, u_mat)
  for (i in seq_len(ncol(u_mat))) {
    u <- u_mat[, i]
    expect_equal(ch$to_theta(u), as.vector(vertices %*% u))
    # The map is linear, so the jacobian is the generator matrix itself.
    expect_equal(ch$jacobian(u), vertices)
  }
  expect_equal(ch$to_theta_batch(u_mat), vertices %*% u_mat)
  # Barycentric coordinates hit the vertices exactly.
  for (j in 1:4) {
    expect_identical(ch$to_theta(diag(4)[, j]), vertices[, j])
  }
})


test_that("the chart is the linear generator map on a simplex", {
  plurality <- simplex_region(
    vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))
  )
  vertices <- plurality@vertices

  ch <- chart(plurality)
  expect_identical(ch$n_par, 3L)
  set.seed(2)
  u_mat <- ch$seed(50L)
  expect_feasible(ch, u_mat)
  for (i in seq_len(ncol(u_mat))) {
    u <- u_mat[, i]
    expect_equal(ch$to_theta(u), as.vector(vertices %*% u))
    expect_equal(ch$jacobian(u), vertices)
  }
})


test_that("the chart anchors a halfspace and frees its hyperplane", {
  for (space in list(
    halfspace_region(normal = c(1, -1), offset = 0),
    halfspace_region(normal = c(1, -1, 0), offset = 2),
    halfspace_region(normal = c(2, 1, -1), offset = -1)
  )) {
    d <- length(space@normal)
    anchor <- space@anchor
    basis <- space@basis
    unit <- space@unit
    # Coordinates are (z, c): the position on the bounding hyperplane and the
    # distance inward along the one ray. A lone vertex contributes no
    # coordinate of its own.
    to_theta <- function(u) {
      anchor + as.vector(basis %*% u[-d]) - u[d] * unit
    }

    ch <- chart(space)
    expect_identical(ch$n_par, d)
    expect_identical(ch$lower, c(rep(-Inf, d - 1L), 0))
    expect_null(ch$heq)
    set.seed(3)
    u_mat <- ch$seed(50L)
    expect_feasible(ch, u_mat)
    for (i in seq_len(ncol(u_mat))) {
      expect_equal(ch$to_theta(u_mat[, i]), to_theta(u_mat[, i]))
      expect_equal(ch$jacobian(u_mat[, i]), cbind(basis, -unit))
    }
  }
})


test_that("the chart of a point region is empty", {
  theta <- c(0.5, 0.3, 0.2)
  space <- point_region(theta = theta)

  ch <- chart(space)
  expect_identical(ch$n_par, 0L)
  expect_identical(ch$to_theta(numeric(0)), theta)
  expect_identical(ch$jacobian(numeric(0)), matrix(0, 3L, 0L))
  expect_identical(ch$from_theta(theta), numeric(0))
  expect_identical(dim(ch$seed(5L)), c(0L, 5L))
})


test_that("the chart of an unconstrained region is the identity", {
  space <- unconstrained_region(3L)

  ch <- chart(space)
  expect_identical(ch$n_par, 3L)
  expect_identical(ch$lower, rep(-Inf, 3L))
  expect_null(ch$heq)
  set.seed(4)
  u_mat <- ch$seed(50L)
  for (i in seq_len(ncol(u_mat))) {
    expect_equal(ch$to_theta(u_mat[, i]), u_mat[, i])
    expect_equal(ch$jacobian(u_mat[, i]), diag(3L))
  }
  expect_equal(ch$to_theta_batch(u_mat), u_mat)
})


test_that("from_theta returns feasible coordinates, exactly at a vertex", {
  # The optimiser is seeded from `from_theta` of the current atoms, so the
  # coordinates must satisfy the declared constraints; SLSQP starts from them.
  s <- simplex_region(
    vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))
  )
  ch <- chart(s)
  interior <- as.vector(s@vertices %*% c(0.2, 0.5, 0.3))
  expect_feasible(ch, matrix(ch$from_theta(interior), ncol = 1L))

  # A vertex is a boundary point the old softmax chart could only approach:
  # its recovery clamped an exact zero to an eps. The direct chart represents
  # it, so the recovery is exact up to the floating point of the
  # least-squares solve rather than short of the boundary by design.
  vertex <- s@vertices[, 2L]
  u <- ch$from_theta(vertex)
  expect_equal(u, c(0, 1, 0), tolerance = 1e-12)
  expect_equal(ch$to_theta(u), vertex, tolerance = 1e-12)
  expect_feasible(ch, matrix(u, ncol = 1L))
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

test_that("the chart round-trips both ways on a full polyhedron", {
  # One region exercising every generator block at once: four vertices, two
  # rays and two lines in R^7, in general position so the seven free
  # directions are independent and the coordinates are identifiable. That
  # makes the round-trip two-sided: coordinates survive to_theta then
  # from_theta, and points survive from_theta then to_theta.
  set.seed(5)
  space <- polyhedron_region(
    vertices = matrix(round(stats::rnorm(7 * 4), 1), nrow = 7),
    rays = matrix(round(stats::rnorm(7 * 2), 1), nrow = 7),
    lines = matrix(round(stats::rnorm(7 * 2), 1), nrow = 7)
  )
  ch <- chart(space)
  expect_identical(ch$n_par, 8L)

  u_mat <- ch$seed(20L)
  expect_feasible(ch, u_mat)
  for (i in seq_len(ncol(u_mat))) {
    u <- u_mat[, i]
    theta <- ch$to_theta(u)
    u_rt <- ch$from_theta(theta)
    expect_equal(u_rt, u, tolerance = 1e-9)
    expect_equal(ch$to_theta(u_rt), theta, tolerance = 1e-9)
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


# --- Printing -----------------------------------------------------------------

test_that("regions print a geometric summary, not a property dump", {
  square <- polytope_region(
    vertices = cbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  )
  expect_output(print(square), "4 vertices in R\\^2")
  expect_output(print(square), "facets: 4 inequalities")
  expect_output(print(square), "vertices:")

  h <- halfspace_region(normal = c(1, -1, 0), offset = 2)
  expect_output(print(h), "unbounded")
  expect_output(print(h), "facets: 1 inequality")

  p <- point_region(theta = c(0.5, 0.3, 0.2))
  expect_output(print(p), "the point \\(0.5, 0.3, 0.2\\)")
  expect_output(print(p), "facets: 3 equalities")

  expect_output(print(unconstrained_region(3L)), "all of R\\^3")
  expect_output(print(unconstrained_region(3L)), "facets: none")

  expect_match(
    format(simplex_region(vertices = diag(3))),
    "3 vertices in R\\^3"
  )

  # A union lists each part by its own format line while short.
  u <- union_region(
    simplex_region(vertices = diag(3)),
    halfspace_region(normal = c(1, -1, 0))
  )
  expect_output(print(u), "simplex_region: 3 vertices")
  expect_output(print(u), "halfspace_region: the halfspace")
})


test_that("facets cannot be declared, only derived", {
  tri_v <- cbind(c(0, 0), c(1, 0), c(0, 1))
  expect_error(
    polytope_region(vertices = tri_v, facets = list(a = 1)),
    "unused argument"
  )
  # The derived views are read-only, so the halves cannot desync after
  # construction.
  tri <- polytope_region(vertices = tri_v)
  expect_error(tri@facets <- list(a = 1), "read-only")
  expect_error(tri@generators <- list(), "read-only")
  expect_error(tri@hv <- list(), "must be <ripr::hv_region>")
  # And the internal door takes only the record type.
  expect_error(
    polyhedron_region(.hv = list(h = NULL, v = list())),
    "must be <ripr::hv_region>"
  )
  expect_error(
    polyhedron_region(vertices = tri_v, .hv = tri@hv),
    "cannot both be given"
  )
})
