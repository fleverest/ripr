# Properties of R/subnull.R.
#
# The geometry is checked by round-trips and idempotence rather than against
# stored coordinates, since a chart is only required to be *a* parametrisation,
# not a particular one.

# The plurality subnull {theta : theta_1 <= theta_j} on the K-simplex, in both
# representations, so the two can be checked against each other.
plurality_simplex <- function(k, j) {
  basis <- lapply(setdiff(seq_len(k), 1L), function(i) {
    v <- numeric(k)
    v[i] <- 1
    v
  })
  tie <- numeric(k)
  tie[c(1L, j)] <- 0.5
  simplex_null(vertices = do.call(cbind, c(basis, list(tie))))
}

plurality_halfspace <- function(k, j) {
  a <- numeric(k)
  a[1L] <- 1
  a[j] <- -1
  halfspace_null(normal = a, offset = 0)
}

# --- Charts -------------------------------------------------------------------

test_that("a chart round-trips points in the subnull", {
  for (s in list(plurality_simplex(4, 2), plurality_halfspace(4, 2))) {
    ch <- chart(s)
    set.seed(1)
    for (i in 1:5) {
      u <- ch$seed(1L)[, 1L]
      theta <- ch$to_theta(u)
      expect_true(contains(s, theta))
      expect_equal(ch$to_theta(ch$from_theta(theta)), theta, tolerance = 1e-6)
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

test_that("chart seeds land inside the subnull", {
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

test_that("projection is idempotent and lands in the subnull", {
  set.seed(5)
  for (s in list(plurality_simplex(4, 2), plurality_halfspace(4, 2))) {
    for (i in 1:5) {
      theta <- stats::runif(4)
      theta <- theta / sum(theta)
      p <- project(s, theta)
      expect_true(contains(s, p))
      expect_equal(project(s, p), p, tolerance = 1e-7)
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

test_that("a simplex subnull contains its own vertices", {
  s <- plurality_simplex(4, 3)
  for (i in seq_len(ncol(s@vertices))) {
    expect_true(contains(s, s@vertices[, i], tol = 1e-6))
  }
})

test_that("a singleton contains only its own point", {
  s <- singleton_null(theta = c(0.5, 0.5))
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

test_that("a maximum at a vertex is approached but not attained", {
  # Documents the limitation rather than working around it: this is why
  # maximise_over returns a lower bound on the supremum, and why a certified
  # upper bound cannot come from this search.
  s <- plurality_simplex(3, 2)
  target <- s@vertices[, 1L] # vertex weight (1, 0, 0)
  obj <- objective(
    value = \(theta) -sum((theta - target)^2),
    grad = \(theta) -2 * (theta - target)
  )
  set.seed(12)
  res <- maximise_over(s, obj, n_seeds = 50L, n_restarts = 5L)

  expect_lt(res$value, 0) # never reaches the supremum of 0
  expect_gt(res$value, -1e-3) # but gets close
  expect_true(contains(s, res$theta, tol = 1e-6))
})

test_that("maximise_over is at least as good as the chart image of its seeds", {
  # The guarantee the duality gap leans on. Stated against the chart image of
  # the seed, not the seed itself: from_theta then to_theta loses O(eps) through
  # the softmax guard, so the search cannot promise to match a seed sitting
  # exactly on a vertex. That residual is why a correctly seeded gap can still
  # come out very slightly negative.
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
  expect_gte(res$value, obj$value(seed_image) - 1e-10)
})

test_that("the chart round-trip is lossless in the interior and lossy at a vertex", {
  # Quantifies the O(eps) above, so a change to the softmax guard shows up here
  # rather than as a mysterious negative gap much later.
  s <- plurality_simplex(4, 2)
  ch <- chart(s)

  interior <- as.vector(s@vertices %*% rep(0.25, 4))
  expect_equal(ch$to_theta(ch$from_theta(interior)), interior, tolerance = 1e-9)

  vertex <- s@vertices[, 1L]
  round_tripped <- ch$to_theta(ch$from_theta(vertex))
  expect_false(isTRUE(all.equal(round_tripped, vertex, tolerance = 0)))
  expect_equal(round_tripped, vertex, tolerance = 1e-6)
})

test_that("maximise_over accepts seeds lying outside the subnull", {
  # Atoms live on other subnulls, so seeds are projected before use.
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
  s <- singleton_null(theta = c(0.5, 0.5))
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

test_that("null_model bundles a family with its subnulls", {
  fam <- multinomial_family(n_trials = 10, k = 4)
  null <- null_model(fam, lapply(2:4, \(j) plurality_simplex(4, j)))
  expect_equal(n_subnulls(null), 3L)
  expect_identical(null@family, fam)
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

test_that("null_model rejects an empty or malformed subnull list", {
  fam <- multinomial_family(n_trials = 4, k = 2)
  expect_error(null_model(fam, list()), "non-empty")
  expect_error(null_model(fam, list("not a subnull")), "must be a `subnull`")
})

test_that("simplex_null rejects a malformed vertex matrix", {
  expect_error(simplex_null(vertices = c(0.5, 0.5)), "must be a matrix")
  expect_error(
    simplex_null(vertices = matrix(numeric(0), nrow = 2L, ncol = 0L)),
    "must be a matrix"
  )
})

test_that("halfspace_null rejects a zero normal", {
  expect_error(halfspace_null(normal = c(0, 0)), "non-zero")
})
