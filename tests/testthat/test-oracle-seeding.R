# Half-space oracle seeding. A half-space is unbounded and has no intrinsic
# scale, so the search must be centred on something problem-derived; these tests
# pin the contract that makes that safe.

quad_objective <- function(f, peak) {
  # A smooth objective with a known unique maximum at `peak`, so the oracle's
  # answer can be checked exactly -- no Monte Carlo noise involved.
  list(
    value = function(theta) -sum((theta - peak)^2),
    grad_theta = function(theta) -2 * (theta - peak),
    value_batch = function(theta_mat) {
      -colSums((theta_mat - peak)^2)
    }
  )
}

test_that("seed_centres steers the search to a distant optimum", {
  f <- halfspace_face(v = c(1, -1), c = 0, face_index = 1)
  u1 <- f@basis[, 1]
  peak <- as.vector(40 * u1) # far from the anchor at the origin
  obj <- quad_objective(f, peak)

  set.seed(3)
  hinted <- replicate(10, {
    r <- oracle(f, obj, n_seeds = 40L, seed_centres = matrix(peak, ncol = 1L))
    sqrt(sum((r$theta - peak)^2))
  })
  # Centred on the right region, the oracle lands on the optimum every time.
  expect_lt(max(hinted), 1e-3)
})

test_that("the anchor is always among the seed centres", {
  # A single centre far from the optimum must not blind the search: the anchor
  # is included, so the region between the support and the anchor is covered.
  f <- halfspace_face(v = c(1, -1), c = 0, face_index = 1)
  u1 <- f@basis[, 1]
  peak <- as.vector(0.5 * u1) # near the anchor
  decoy <- matrix(as.vector(40 * u1), ncol = 1L) # the only atom, far away
  obj <- quad_objective(f, peak)

  set.seed(4)
  got <- replicate(10, {
    r <- oracle(f, obj, n_seeds = 40L, seed_centres = decoy)
    sqrt(sum((r$theta - peak)^2))
  })
  expect_lt(max(got), 1e-3)
})

test_that("centres on other faces are projected rather than rejected", {
  f2 <- halfspace_face(v = c(1, 1), c = 0, face_index = 2)
  # This point lies on face 1 and is NOT in face 2.
  foreign <- matrix(c(1.548, 1.548), ncol = 1L)
  expect_false(contains(f2, as.vector(foreign)))

  obj <- quad_objective(f2, as.vector(project(f2, c(-2, -1))))
  set.seed(5)
  r <- oracle(f2, obj, n_seeds = 30L, seed_centres = foreign)
  expect_true(all(is.finite(r$theta)))
  expect_true(contains(f2, r$theta))
  expect_true(is.finite(r$value))
})

test_that("many centres are handled under a capped polish budget", {
  f <- halfspace_face(v = c(1, -1), c = 0, face_index = 1)
  u1 <- f@basis[, 1]
  # 40 centres but only 25 BFGS polishes available.
  centres <- do.call(cbind, lapply(seq(-20, 20, length.out = 40),
                                   function(z) as.vector(z * u1)))
  peak <- as.vector(7 * u1)
  obj <- quad_objective(f, peak)

  set.seed(6)
  r <- oracle(f, obj, n_seeds = 10L, seed_centres = centres, n_restarts = 25L)
  expect_true(contains(f, r$theta))
  expect_lt(sqrt(sum((r$theta - peak)^2)), 1e-3)
})

test_that("n_seeds is a floor: more centres draw proportionally more seeds", {
  f <- halfspace_face(v = c(1, -1), c = 0, face_index = 1)
  u1 <- f@basis[, 1]
  peak <- as.vector(2 * u1)
  n_evaluated <- 0L
  obj <- quad_objective(f, peak)
  counting <- list(
    value = obj$value,
    grad_theta = obj$grad_theta,
    value_batch = function(theta_mat) {
      n_evaluated <<- n_evaluated + ncol(theta_mat)
      obj$value_batch(theta_mat)
    }
  )
  centres <- do.call(cbind, lapply(seq(-5, 5, length.out = 20),
                                   function(z) as.vector(z * u1)))
  set.seed(7)
  oracle(f, counting, n_seeds = 10L, seed_centres = centres)
  # 21 centres (20 atoms + anchor) x 8 minimum each = 168 >> the n_seeds floor.
  expect_gte(n_evaluated, 21L * 8L)
})

test_that("polytope faces ignore seed_centres", {
  face <- polytope_face(vertices = cbind(c(0, 1), c(0.5, 0.5)), face_index = 1)
  obj <- quad_objective(face, c(0.5, 0.5))
  set.seed(8)
  a <- oracle(face, obj, n_seeds = 40L)
  set.seed(8)
  b <- oracle(face, obj, n_seeds = 40L, seed_centres = matrix(c(0, 1), ncol = 1L))
  expect_equal(a$theta, b$theta)
  expect_equal(a$value, b$value)
})

test_that("the seeded search is unaffected when a single seed_alpha is given", {
  # The EM M-step path must stay a pure local refinement.
  f <- halfspace_face(v = c(1, -1), c = 0, face_index = 1)
  peak <- as.vector(3 * f@basis[, 1])
  obj <- quad_objective(f, peak)
  start <- ripr:::face_coordinates(f, as.vector(2.5 * f@basis[, 1]))
  set.seed(9)
  a <- oracle(f, obj, seed_alpha = start)
  set.seed(9)
  b <- oracle(f, obj, seed_alpha = start,
              seed_centres = matrix(as.vector(-30 * f@basis[, 1]), ncol = 1L))
  expect_equal(a$theta, b$theta)
})
