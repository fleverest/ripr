# Properties of R/region.R and R/null.R.
#
# The geometry is checked by round-trips and idempotence rather than against
# stored coordinates, since a chart is only required to be *a* parametrisation,
# not a particular one.

# The plurality part {theta : theta_1 <= theta_j} on the K-simplex, in both
# representations, so the two can be checked against each other.
plurality_simplex <- function(k, j) {
  basis <- lapply(setdiff(seq_len(k), 1L), function(i) {
    v <- numeric(k)
    v[i] <- 1
    v
  })
  tie <- numeric(k)
  tie[c(1L, j)] <- 0.5
  simplex_region(vertices = do.call(cbind, c(basis, list(tie))))
}

plurality_halfspace <- function(k, j) {
  a <- numeric(k)
  a[1L] <- 1
  a[j] <- -1
  halfspace_region(normal = a, offset = 0)
}

# --- Charts -------------------------------------------------------------------

test_that("a chart round-trips points in the part", {
  for (s in list(plurality_simplex(4, 2), plurality_halfspace(4, 2))) {
    ch <- chart(s)
    set.seed(1)
    for (i in 1:5) {
      u <- ch$seed(1L)[, 1L]
      theta <- ch$to_theta(u)
      expect_true(contains(s, theta))
      expect_equal(
        ch$to_theta(ch$from_theta(theta)),
        theta,
        tolerance = rounding_tol(1)
      )
    }
  }
})

test_that("to_theta_batch agrees with to_theta column by column", {
  for (s in list(plurality_simplex(4, 3), plurality_halfspace(4, 3))) {
    ch <- chart(s)
    set.seed(2)
    u_mat <- ch$seed(6L)
    batched <- ch$to_theta_batch(u_mat)
    expect_equal(ncol(batched), 6L)
    for (i in 1:6) {
      expect_equal(batched[, i], ch$to_theta(u_mat[, i]))
    }
  }
})

test_that("the chart Jacobian matches a finite difference", {
  for (s in list(plurality_simplex(4, 2), plurality_halfspace(3, 2))) {
    ch <- chart(s)
    set.seed(3)
    u <- ch$seed(1L)[, 1L]
    jac <- ch$jacobian(u)
    eps <- 1e-6
    for (j in seq_len(ch$n_par)) {
      e <- numeric(ch$n_par)
      e[j] <- eps
      fd <- (ch$to_theta(u + e) - ch$to_theta(u - e)) / (2 * eps)
      expect_equal(jac[, j], fd, tolerance = 1e-5)
    }
  }
})

test_that("chart seeds land inside the part", {
  for (s in list(plurality_simplex(5, 4), plurality_halfspace(5, 4))) {
    ch <- chart(s)
    set.seed(4)
    u_mat <- ch$seed(40L)
    expect_equal(nrow(u_mat), ch$n_par)
    thetas <- ch$to_theta_batch(u_mat)
    expect_true(all(apply(thetas, 2L, \(t) contains(s, t))))
  }
})

# --- Membership and projection ------------------------------------------------

test_that("projection is idempotent and lands in the part", {
  set.seed(5)
  for (s in list(plurality_simplex(4, 2), plurality_halfspace(4, 2))) {
    for (i in 1:5) {
      theta <- stats::runif(4)
      theta <- theta / sum(theta)
      p <- project(s, theta)
      expect_true(contains(s, p))
      expect_equal(project(s, p), p, tolerance = rounding_tol(1))
    }
  }
})

test_that("projection leaves points already inside untouched", {
  s <- plurality_halfspace(3, 2)
  theta <- c(0.2, 0.5, 0.3) # theta_1 <= theta_2
  expect_true(contains(s, theta))
  expect_equal(project(s, theta), theta)
})

test_that("the two representations agree on membership within the simplex", {
  # The vertex hull and the halfspace describe the same set once intersected
  # with the probability simplex, so they must agree on simplex points.
  set.seed(6)
  hull <- plurality_simplex(4, 2)
  half <- plurality_halfspace(4, 2)
  for (i in 1:40) {
    theta <- stats::rgamma(4, 1)
    theta <- theta / sum(theta)
    expect_equal(
      contains(hull, theta, tol = 1e-6),
      contains(half, theta, tol = 1e-6)
    )
  }
})

test_that("a simplex part contains its own vertices", {
  s <- plurality_simplex(4, 3)
  for (i in seq_len(ncol(s@vertices))) {
    expect_true(contains(s, s@vertices[, i], tol = 1e-6))
  }
})

test_that("a singleton contains only its own point", {
  s <- point_region(theta = c(0.5, 0.5))
  expect_true(contains(s, c(0.5, 0.5)))
  expect_false(contains(s, c(0.6, 0.4)))
  expect_equal(project(s, c(0.9, 0.1)), c(0.5, 0.5))
  expect_equal(chart(s)$n_par, 0L)
})

# --- The oracle ---------------------------------------------------------------

test_that("maximise_over finds a maximum interior to the chart", {
  # The peak must be interior in *vertex-weight* coordinates, not merely inside
  # the set: the softmax chart covers only the relative interior, so a target
  # with a zero vertex weight is approached and never reached. Building the
  # target as an interior convex combination guarantees the chart can attain it.
  s <- plurality_simplex(3, 2)
  target <- as.vector(s@vertices %*% rep(1 / 3, 3))
  obj <- objective(
    value = \(theta) -sum((theta - target)^2),
    grad = \(theta) -2 * (theta - target)
  )
  set.seed(7)
  res <- maximise_over(s, obj, n_seeds = 50L, n_restarts = 5L)
  expect_equal(res$theta, target, tolerance = 1e-5)
  expect_equal(res$value, 0, tolerance = 1e-9)
  expect_true(contains(s, res$theta, tol = 1e-6))
})

test_that("a maximum at a vertex is attained exactly", {
  # The old softmax chart covered only the relative interior, so a vertex
  # maximum was approached and never reached. The direct chart represents the
  # vertex and SLSQP's active constraints pin it exactly. The objective here
  # is concave, so the multistart is guaranteed the right basin and the test
  # is deterministic; for a non-convex oracle the exactness holds only once
  # the search finds the right face, and maximise_over remains a lower bound
  # on the supremum.
  s <- plurality_simplex(3, 2)
  target <- s@vertices[, 1L] # vertex weight (1, 0, 0)
  obj <- objective(
    value = \(theta) -sum((theta - target)^2),
    grad = \(theta) -2 * (theta - target)
  )
  set.seed(12)
  res <- maximise_over(s, obj, n_seeds = 50L, n_restarts = 5L)

  expect_equal(res$value, 0, tolerance = rounding_tol(0))
  expect_equal(res$theta, target, tolerance = 1e-9)
  expect_true(contains(s, res$theta, tol = 1e-6))
})

test_that("maximise_over is at least as good as the chart image of its seeds", {
  # The guarantee the duality gap leans on. Stated against the chart image of
  # the seed, not the seed itself: the direct chart makes the two agree to
  # floating point, but the guarantee is about what the search was actually
  # given, so the statement stays in this form.
  s <- plurality_simplex(4, 2)
  peak <- as.vector(s@vertices %*% c(0.7, 0.1, 0.1, 0.1))
  obj <- objective(
    value = \(theta) -sum((theta - peak)^2),
    grad = \(theta) -2 * (theta - peak)
  )
  ch <- chart(s)
  seed_pt <- project(s, peak)
  seed_image <- ch$to_theta(ch$from_theta(seed_pt))

  set.seed(8)
  res <- maximise_over(
    s,
    obj,
    seeds = matrix(seed_pt, ncol = 1L),
    n_seeds = 10L,
    n_restarts = 3L
  )
  expect_gte(res$value, obj$value(seed_image) - rounding_tol(1))
})

test_that("the chart round-trip is lossless in the interior and at a vertex", {
  # The direct chart represents boundary points exactly -- there is no softmax
  # guard losing O(eps) at a vertex, so a correctly seeded gap can no longer
  # come out negative through the round-trip.
  s <- plurality_simplex(4, 2)
  ch <- chart(s)

  interior <- as.vector(s@vertices %*% rep(0.25, 4))
  expect_equal(
    ch$to_theta(ch$from_theta(interior)),
    interior,
    tolerance = rounding_tol(1)
  )

  # Exact up to the floating point of the least-squares recovery: bit-identical
  # on some vertex matrices, an ulp or two off on others, and BLAS-dependent
  # either way -- so tested at 1e-12, not identical().
  vertex <- s@vertices[, 1L]
  expect_equal(
    ch$to_theta(ch$from_theta(vertex)),
    vertex,
    tolerance = rounding_tol(1)
  )
})

test_that("maximise_over accepts seeds lying outside the part", {
  # Atoms live on other parts, so seeds are projected before use.
  s <- plurality_halfspace(3, 2)
  target <- c(0.2, 0.6, 0.2)
  obj <- objective(
    value = \(theta) -sum((theta - target)^2),
    grad = \(theta) -2 * (theta - target)
  )
  outside <- matrix(c(0.9, 0.05, 0.05), ncol = 1L)
  expect_false(contains(s, outside[, 1L]))
  set.seed(9)
  res <- maximise_over(s, obj, seeds = outside, n_seeds = 20L, n_restarts = 3L)
  expect_true(contains(s, res$theta, tol = 1e-6))
})

test_that("maximise_over on a singleton evaluates the point", {
  s <- point_region(theta = c(0.5, 0.5))
  obj <- objective(value = \(theta) sum(theta^2), grad = \(theta) 2 * theta)
  res <- maximise_over(s, obj)
  expect_equal(res$theta, c(0.5, 0.5))
  expect_equal(res$value, 0.5)
})

test_that("objective supplies a working default batch evaluator", {
  obj <- objective(value = \(theta) sum(theta), grad = \(theta) {
    rep(1, length(theta))
  })
  m <- cbind(c(1, 2), c(3, 4), c(5, 6))
  expect_equal(obj$value_batch(m), c(3, 7, 11))
})

test_that("a supplied batch evaluator is used", {
  called <- 0L
  obj <- objective(
    value = \(theta) sum(theta),
    grad = \(theta) rep(1, length(theta)),
    value_batch = function(m) {
      called <<- called + 1L
      colSums(m)
    }
  )
  expect_equal(obj$value_batch(cbind(c(1, 2))), 3)
  expect_equal(called, 1L)
})

# --- null_model ---------------------------------------------------------------

test_that("null_model bundles a family with its parts", {
  fam <- multinomial_family(n_trials = 10, k = 4)
  null <- null_model(fam, lapply(2:4, \(j) plurality_simplex(4, j)))
  expect_equal(n_parts(null@region), 3L)
  expect_identical(null@family, fam)
})

test_that("a null takes its geometry as a part, a list or a union alike", {
  # `null_model()` coerces, so none of the three needs the caller to know
  # which form the others take.
  fam <- multinomial_family(n_trials = 10, k = 3)
  s1 <- plurality_simplex(3, 2)
  s2 <- plurality_simplex(3, 3)
  theta <- c(0.3, 0.5, 0.2)

  from_list <- null_model(fam, list(s1, s2))
  from_union <- null_model(fam, union_region(s1, s2))
  expect_equal(n_parts(from_list@region), 2L)
  expect_equal(n_parts(from_union@region), 2L)
  expect_identical(from_list@region, from_union@region)
  expect_equal(in_null(from_list, theta), in_null(from_union, theta))

  # A lone convex region is already a region, so it is stored as it came --
  # not wrapped in a one-element union.
  bare <- null_model(fam, s1)
  wrapped <- null_model(fam, list(s1))
  expect_identical(bare@region, s1)
  expect_identical(wrapped@region, s1)
  expect_equal(n_parts(bare@region), 1L)
  expect_equal(in_null(bare, theta), in_null(wrapped, theta))
})

test_that("a null takes its decomposition once, at construction", {
  fam <- multinomial_family(n_trials = 10, k = 3)
  # A part that is already a simplex is its own cell, and the same object:
  # nothing has been rebuilt behind the caller's back.
  simplices <- null_model(fam, lapply(2:3, \(j) plurality_simplex(3, j)))
  expect_identical(simplices@cells, parts(simplices@region))
  expect_identical(simplices@cell_part, 1:2)

  # A part that is a convex hull is several cells, all filed under it.
  square <- polytope_region(
    vertices = cbind(
      c(0.5, 0.5, 0),
      c(0, 0.5, 0.5),
      c(0, 0, 1),
      c(0.5, 0, 0.5)
    )
  )
  mixed <- null_model(fam, list(plurality_simplex(3, 2), square))
  expect_length(mixed@cells, 3L)
  expect_identical(mixed@cell_part, c(1L, 2L, 2L))
  expect_true(all(vapply(
    mixed@cells,
    \(c) S7_inherits(c, simplex_region),
    logical(1)
  )))

  # The cells are in part order and are that part's own `cells()`.
  expect_identical(mixed@cells[[1L]], parts(mixed@region)[[1L]])
  expect_identical(mixed@cells[2:3], cells(square))
})

test_that("null_model()'s max_cells caps the decomposition across parts", {
  fam <- multinomial_family(n_trials = 10, k = 3)
  square <- polytope_region(
    vertices = cbind(
      c(0.5, 0.5, 0),
      c(0, 0.5, 0.5),
      c(0, 0, 1),
      c(0.5, 0, 0.5)
    )
  )
  # Two triangles per square, so two squares need four in total.
  null <- null_model(fam, list(square, square), max_cells = 4L)
  expect_length(null@cells, 4L)
  # One budget spans the parts, and the refusal still names the one that
  # exhausted it.
  expect_error(
    null_model(fam, list(square, square), max_cells = 3L),
    "could not decompose part 2.*max_cells = 3"
  )
})

test_that("in_null is the union of the pieces", {
  fam <- multinomial_family(n_trials = 10, k = 3)
  null <- null_model(fam, lapply(2:3, \(j) plurality_halfspace(3, j)))
  # theta_1 largest: outside every piece, so outside the null.
  expect_false(in_null(null, c(0.5, 0.3, 0.2)))
  # theta_1 not the strict plurality winner: inside at least one piece.
  expect_true(in_null(null, c(0.2, 0.5, 0.3)))
  expect_true(in_null(null, c(1 / 3, 1 / 3, 1 / 3)))
})

test_that("null_model rejects an empty or malformed part list", {
  fam <- multinomial_family(n_trials = 4, k = 2)
  expect_error(null_model(fam, list()), "non-empty")
  expect_error(
    null_model(fam, list("not a region")),
    "must be a `convex_region`"
  )
})

test_that("simplex_region rejects a malformed vertex matrix", {
  expect_error(simplex_region(vertices = c(0.5, 0.5)), "must be a matrix")
  expect_error(
    simplex_region(vertices = matrix(numeric(0), nrow = 2L, ncol = 0L)),
    "must be a matrix"
  )
})

test_that("halfspace_region rejects a zero normal", {
  expect_error(halfspace_region(normal = c(0, 0)), "non-zero")
})

# --- Dimension ----------------------------------------------------------------

test_that("every region reports the dimension of the parameters it holds", {
  expect_equal(space_dim(simplex_region(vertices = diag(3))), 3L)
  expect_equal(space_dim(halfspace_region(normal = c(1, -1, 0))), 3L)
  expect_equal(space_dim(point_region(theta = c(0.5, 0.3, 0.2))), 3L)
  expect_equal(space_dim(unconstrained_region(3L)), 3L)
})

test_that("region_dim is the affine dimension, at most the ambient", {
  # `space_dim` is how many coordinates a point carries; `region_dim` is what
  # the geometry spans. Constraints separate the two.
  expect_identical(region_dim(point_region(theta = c(0.5, 0.3, 0.2))), 0L)
  segment <- polytope_region(vertices = cbind(c(0, 0, 0), c(1, 1, 1)))
  expect_identical(region_dim(segment), 1L)
  expect_identical(region_dim(simplex_region(vertices = diag(3))), 2L)
  expect_identical(region_dim(halfspace_region(normal = c(1, -1, 0))), 3L)
  expect_identical(region_dim(unconstrained_region(2L)), 2L)

  # A union spans what its largest part spans.
  tri <- simplex_region(
    vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))
  )
  expect_identical(region_dim(union_region(segment, tri)), 2L)

  # The empty set sits strictly below a point.
  expect_identical(region_dim(empty_region()), -1L)
})

test_that("a region of the wrong dimension is refused at construction", {
  # Without this the comparison inside `contains()` recycles instead of
  # complaining, and `in_null()` returns TRUE for a parameter it never checked:
  #   null_model(multinomial_family(4, 3), list(halfspace_region(c(1, -1))))
  #   in_null(null, c(0.2, 0.5, 0.3))  # TRUE, silently wrong
  fam <- multinomial_family(n_trials = 4, k = 3)
  expect_error(
    null_model(fam, list(halfspace_region(normal = c(1, -1)))),
    "dimension 3"
  )
  expect_error(
    null_model(fam, list(point_region(theta = c(0.5, 0.5)))),
    "dimension 3"
  )
  # Parts that disagree with each other fail earlier and for a better reason:
  # `region`'s validator refuses them before the family is consulted at all,
  # since a union of a plane and a line has no ambient space to live in.
  expect_error(
    null_model(
      fam,
      list(simplex_region(vertices = diag(3)), unconstrained_region(2L))
    ),
    "same ambient dimension"
  )
  expect_silent(null_model(fam, list(unconstrained_region(3L))))
})

test_that("a family's parameter space has the family's own dimension", {
  # `param_dim()` used to answer this; `space_dim()` on the parameter space
  # subsumes it, and unlike a per-family method it cannot disagree with the
  # geometry the null is checked against.
  fam <- multinomial_family(n_trials = 7, k = 5)
  expect_equal(space_dim(fam@parameter_space), 5L)
  expect_true(contains(fam@parameter_space, rep(1 / 5, 5)))
  expect_false(contains(fam@parameter_space, c(0.5, 0.5, 0.5, 0.5, 0.5)))
})

# --- The unconstrained region -------------------------------------------------

test_that("a real region contains every finite point and moves none", {
  r <- unconstrained_region(2L)
  expect_true(contains(r, c(1e6, -3)))
  expect_false(contains(r, c(Inf, 0)))
  expect_equal(project(r, c(3, -1)), c(3, -1))
})

test_that("the real region's chart is the identity", {
  r <- unconstrained_region(3L)
  ch <- chart(r)
  expect_equal(ch$n_par, 3L)
  theta <- c(0.4, -2, 7)
  expect_equal(ch$to_theta(ch$from_theta(theta)), theta)
  expect_equal(ch$jacobian(theta), diag(3))
  expect_equal(dim(ch$seed(5L)), c(3L, 5L))
  u <- cbind(c(1, 2, 3), c(-1, 0, 1))
  expect_equal(ch$to_theta_batch(u), u)
})

test_that("maximise_over finds an interior optimum on a real region", {
  # An interior optimum on an unbounded region: no constraint is active, so
  # this exercises the plain quasi-Newton behaviour of the refinement.
  set.seed(1)
  target <- c(1.5, -0.5)
  obj <- objective(
    value = function(theta) -sum((theta - target)^2),
    grad = function(theta) -2 * (theta - target)
  )
  found <- maximise_over(
    unconstrained_region(2L),
    obj,
    n_seeds = 50L,
    n_restarts = 5L
  )
  expect_equal(found$theta, target, tolerance = 1e-6)
  expect_equal(found$value, 0, tolerance = 1e-10)
})
