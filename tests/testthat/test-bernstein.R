# The Bernstein machinery, tested as numerics rather than as statistics. Two
# things matter here and nothing else does: the subdivision is *exact* (its
# child coefficients are the true Bernstein form of the same polynomial, so the
# bound stays valid however deep it goes), and the coefficient indexing is not
# quietly assumed to match the engine's outcome enumeration.

# G evaluated the slow, honest way: sum over the support of coefficient times
# multinomial pmf. Uses the package's -Inf-safe multiply so a vertex with a zero
# component (0 * log 0) does not become NaN.
bern_eval <- function(coef, lat, theta) {
  log_p <- ripr:::log_multinom_coef(lat$tally, lat$n) +
    as.vector(ripr:::matmul_0_ninf(lat$tally, log(theta)))
  sum(coef * exp(log_p))
}

test_that("split_1d subdivides a 1-D Bernstein array at the midpoint", {
  lr <- ripr:::split_1d(c(0.4, 1.8, 0.8))
  expect_equal(lr$left, c(0.4, 1.1, 1.2))
  expect_equal(lr$right, c(1.2, 1.3, 0.8))

  # Degenerate degrees are pass-throughs, not errors: a lattice line of length
  # one appears whenever the spectator indices already use up the whole degree.
  expect_equal(ripr:::split_1d(c(2.5)), list(left = 2.5, right = 2.5))
})

test_that("bisection is exact: child vertex coefficients are G at the vertices", {
  set.seed(7)
  for (K in 2:4) {
    n <- 5L
    lat <- ripr:::bernstein_lattice(n, K)
    coef <- runif(lat$n_coef, 0.05, 4)
    root <- list(V = diag(K), coef = coef)

    kids <- ripr:::bisect(root, 1L, K, lat)
    for (kid in kids) {
      direct <- apply(kid$V, 2L, function(th) bern_eval(coef, lat, th))
      expect_equal(kid$coef[lat$vertex], direct, tolerance = 1e-12)
      # And the child still bounds the polynomial it describes.
      expect_lte(max(direct), max(kid$coef))
    }
  }
})

test_that("the bound is order-agnostic in the tally enumeration", {
  set.seed(11)
  n <- 5L
  K <- 3L
  canonical <- ripr:::compositions(n, K)
  perm <- sample(nrow(canonical))
  permuted <- canonical[perm, , drop = FALSE]

  # One objective, expressed against each enumeration.
  coef_canonical <- runif(nrow(canonical), 0.05, 4)
  coef_permuted <- coef_canonical[perm]

  V <- cbind(c(0, 1, 0), c(0, 0, 1), c(0.5, 0.5, 0))

  run <- function(tally, coef) {
    lat <- ripr:::bernstein_lattice(n, K, tally = tally)
    box <- list(V = V, coef = ripr:::reparametrise_to(coef, lat, V))
    ripr:::certify_sup(list(box), lat, tol = 1e-9, max_iter = 300L)
  }

  a <- run(canonical, coef_canonical)
  b <- run(permuted, coef_permuted)
  d <- run(ripr:::build_counts_matrix(n, K), coef_canonical)

  expect_equal(a$bound, b$bound, tolerance = 1e-9)
  expect_equal(a$bound, d$bound, tolerance = 1e-9)
  expect_equal(a$incumbent, b$incumbent, tolerance = 1e-9)
})

test_that("a malformed tally is rejected rather than silently mis-indexed", {
  n <- 4L
  K <- 3L
  good <- ripr:::compositions(n, K)

  expect_error(ripr:::bernstein_lattice(n, K, tally = good[, 1:2]))
  expect_error(ripr:::bernstein_lattice(n, K, tally = good[-1L, ]))
  expect_error(
    ripr:::bernstein_lattice(n, K, tally = rbind(good[-1L, ], good[1L, ] + 1L))
  )
  # Duplicated rows: right shape, right row sums, wrong set.
  dup <- good
  dup[1L, ] <- dup[2L, ]
  expect_error(ripr:::bernstein_lattice(n, K, tally = dup))
})

test_that("build_counts_matrix() and compositions() agree row for row", {
  # The two enumerations are independently maintained -- stars-and-bars in
  # multinomial.R, a prefix recursion here -- and their agreement is load
  # bearing: the Bernstein coefficients are indexed by `engine@outcomes`, which
  # comes from the former. This pins the coincidence so a change to either
  # fails loudly here instead of producing a plausible wrong bound.
  for (n in 1:12) {
    for (K in 2:6) {
      if (choose(n + K - 1L, K - 1L) > 3000) next
      expect_identical(
        ripr:::build_counts_matrix(n, K),
        ripr:::compositions(n, K),
        info = sprintf("n = %d, K = %d", n, K)
      )
    }
  }
})

test_that("reparametrise_to() agrees with bisect() on a plurality face", {
  set.seed(3)
  n <- 6L
  K <- 3L
  lat <- ripr:::bernstein_lattice(n, K)
  coef <- runif(lat$n_coef, 0.05, 4)
  root <- list(V = diag(K), coef = coef)

  # The plurality face {theta_1 <= theta_2} is one bisection of edge (1, 2),
  # keeping vertex 2 -- so the cheap path and the general path must coincide.
  for (i in 2:K) {
    kid <- ripr:::bisect(root, 1L, i, lat)[[2L]]
    general <- ripr:::reparametrise_to(coef, lat, kid$V)
    expect_equal(general, kid$coef, tolerance = 1e-12)
  }
})

test_that("reparametrise_to() reproduces the polynomial on a general simplex", {
  set.seed(5)
  n <- 4L
  K <- 3L
  lat <- ripr:::bernstein_lattice(n, K)
  coef <- runif(lat$n_coef, 0.05, 4)

  # The identity map must be the identity.
  expect_equal(ripr:::reparametrise_to(coef, lat, diag(K)), coef)

  # A simplex reachable by no sequence of edge bisections.
  V <- cbind(c(0.5, 0.3, 0.2), c(0.1, 0.8, 0.1), c(0.25, 0.25, 0.5))
  rp <- ripr:::reparametrise_to(coef, lat, V)

  # Vertex interpolation: the corner coefficients are exact values of G.
  expect_equal(
    rp[lat$vertex],
    apply(V, 2L, function(th) bern_eval(coef, lat, th)),
    tolerance = 1e-12
  )

  # And an interior point, evaluated through the new form, matches G.
  a <- c(0.2, 0.5, 0.3)
  lat_a <- ripr:::bernstein_lattice(n, K)
  expect_equal(
    bern_eval(rp, lat_a, a),
    bern_eval(coef, lat, as.vector(V %*% a)),
    tolerance = 1e-12
  )
})

test_that("certify_sup() brackets the true supremum and stays valid early", {
  set.seed(21)
  n <- 6L
  K <- 3L
  lat <- ripr:::bernstein_lattice(n, K)
  coef <- runif(lat$n_coef, 0.05, 4)
  box <- list(V = diag(K), coef = coef)

  grid <- as.matrix(expand.grid(a = seq(0, 1, by = 0.005), b = seq(0, 1, by = 0.005)))
  grid <- grid[rowSums(grid) <= 1, , drop = FALSE]
  theta <- cbind(grid, 1 - rowSums(grid))
  gmax <- max(apply(theta, 1L, function(th) bern_eval(coef, lat, th)))

  # Validity does not depend on convergence: the bound holds at iteration 0.
  for (mi in c(0L, 1L, 5L, 50L, 500L)) {
    res <- ripr:::certify_sup(list(box), lat, tol = 1e-10, max_iter = mi)
    expect_lte(res$incumbent, res$bound)
    expect_gte(res$bound, gmax)
    expect_lte(res$incumbent, gmax + 1e-12)
    expect_lte(res$iterations, mi)
  }

  # And it does converge: at a tight tolerance the bracket closes on gmax.
  tight <- ripr:::certify_sup(list(box), lat, tol = 1e-8, max_iter = 5000L)
  expect_equal(tight$bound, gmax, tolerance = 1e-5)
  expect_true(ripr:::contains(
    polytope_face(vertices = diag(K), face_index = 1),
    tight$theta
  ))
})

test_that("keep_argmax = TRUE forbids a positive slack", {
  lat <- ripr:::bernstein_lattice(3L, 3L)
  box <- list(V = diag(3L), coef = rep(1, lat$n_coef))
  expect_error(
    ripr:::certify_sup(list(box), lat, keep_argmax = TRUE, slack = 0.1),
    "must be 0"
  )
})

test_that("the cost model runs and is the pessimistic one", {
  set.seed(9)
  lat <- ripr:::bernstein_lattice(5L, 3L)
  coef <- runif(lat$n_coef, 0.05, 4)
  box <- list(V = diag(3L), coef = coef)

  sdn <- ripr:::second_difference_norm(box, lat)
  expect_true(is.finite(sdn))
  expect_gte(sdn, 0)

  cost <- ripr:::leroy_cost(lat, sdn, tol = 1e-3)
  expect_true(all(is.finite(unlist(cost))))
  # It ignores the cut-off test entirely, so it must sit far above what the
  # branch and bound actually spends. This is why it is not user-facing.
  measured <- ripr:::certify_sup(list(box), lat, tol = 1e-3, max_iter = 5000L)
  expect_gt(cost$subsimplices, measured$iterations)
})
