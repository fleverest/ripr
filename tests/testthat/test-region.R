# Tests for spaces and regions in R/region.R.
#
# The split between `polytope_region` and `simplex_region` exists so that
# `certify()` can dispatch on a validated class.

# --- The validator ------------------------------------------------------------

test_that("simplex vertices are affinely independent", {
  # So all of these are simplices:
  # the 2-simplex in Delta,
  expect_no_error(simplex_region(vertices = diag(3)))
  # a segment inside it,
  expect_no_error(simplex_region(vertices = cbind(c(1, 0, 0), c(0, 1, 0))))
  # a tetrahedron in R^3
  expect_no_error(simplex_region(
    vertices = cbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
  ))

  # And these are not:
  # affinely dependent,
  expect_error(
    simplex_region(vertices = cbind(c(0, 0, 0), c(1, 0, 0), c(2, 0, 0))),
    "affinely dependent"
  )
  # too many points to be independent in this ambient,
  expect_error(
    simplex_region(vertices = cbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1))),
    "at most 3 points can be affinely independent in 2 dimensions"
  )
  # non-finite
  expect_error(
    simplex_region(vertices = cbind(c(Inf, 0, 0), c(0, 1, 0))),
    "finite"
  )
  expect_error(
    simplex_region(vertices = cbind(c(NaN, 0, 0), c(0, 1, 0))),
    "finite"
  )
})

test_that("simplex degeneracy test works in high dimensions", {
  for (k in c(3L, 8L, 30L)) {
    vertices <- diag(k)
    vertices[, 1L] <- replace(numeric(k), c(1L, 2L), 0.5)
    expect_no_error(simplex_region(vertices = vertices))
  }

  tiny <- 1e-6 * cbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
  expect_lt(abs(det(tiny)), 1e-12)
  expect_no_error(simplex_region(vertices = tiny))
})

test_that("polytope_region accepts simplices and non-simplices", {
  expect_no_error(
    polytope_region(vertices = diag(3))
  )

  square <- cbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  expect_no_error(polytope_region(vertices = square))
  expect_false(S7_inherits(polytope_region(vertices = square), simplex_region))

  # Redundant and degenerate vertex sets are polytopes too.
  expect_no_error(polytope_region(vertices = cbind(diag(3), c(0.5, 0.5, 0))))
  expect_no_error(
    polytope_region(vertices = cbind(c(1, 0, 0), c(0.5, 0.5, 0), c(0, 1, 0)))
  )
})

test_that("polytope_region still rejects a non-matrix or empty vertex set", {
  expect_error(polytope_region(vertices = c(0.5, 0.5)), "must be a matrix")
  expect_error(
    polytope_region(vertices = matrix(numeric(0), nrow = 2L, ncol = 0L)),
    "must be a matrix"
  )
})

test_that("a non-finite vertex is named rather than crashing the SVD", {
  # Left to `svd()` this is an unattributed LAPACK error, and it fires before
  # `simplex_region`'s validator gets a chance to say anything useful.
  bad <- cbind(c(NA, 1, 0), c(0, 1, 0), c(0, 0, 1))
  expect_error(polytope_region(vertices = bad), "must all be finite")
  expect_error(simplex_region(vertices = bad), "must all be finite")
  expect_error(
    polytope_region(vertices = cbind(c(Inf, 0), c(0, 1))),
    "must all be finite"
  )
})

# --- The geometry is unchanged ------------------------------------------------

test_that("the subclass inherits every method unchanged", {
  # Same vertex matrix through both constructors: the split is a refactor, so
  # nothing about the geometry may depend on which class holds it.
  p <- polytope_region(vertices = diag(3))
  s <- simplex_region(vertices = diag(3))

  expect_equal(space_dim(p), space_dim(s))

  cp <- chart(p)
  cs <- chart(s)
  expect_equal(cp$n_par, cs$n_par)
  u <- c(0.3, -0.7)
  expect_equal(cp$to_theta(u), cs$to_theta(u))
  expect_equal(cp$jacobian(u), cs$jacobian(u))
  theta <- c(0.2, 0.5, 0.3)
  expect_equal(cp$from_theta(theta), cs$from_theta(theta))

  for (point in list(c(2, -1, 0), c(1 / 3, 1 / 3, 1 / 3), c(0.7, 0.1, 0.2))) {
    expect_equal(project(p, point), project(s, point))
    expect_equal(contains(p, point), contains(s, point))
  }

  expect_equal(init_point(p, c(1, 0, 0)), init_point(s, c(1, 0, 0)))
})

# --- parts() and cells() ------------------------------------------------------

test_that("every convex region is its own only part and its own only cell", {
  spaces <- list(
    simplex_region(vertices = diag(3)),
    polytope_region(vertices = cbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))),
    halfspace_region(normal = c(1, -1), offset = 0),
    point_region(theta = c(0.5, 0.3, 0.2)),
    unconstrained_region(2L)
  )
  for (space in spaces) {
    for (decomposition in list(parts(space), cells(space))) {
      expect_type(decomposition, "list")
      expect_length(decomposition, 1L)
      expect_identical(decomposition[[1L]], space)
    }
  }
})
