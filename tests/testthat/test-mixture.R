# Properties of R/mixture.R.
#
# A mixture is the pushforward of a mixing measure through a family, so these
# check the integral and its sampler agree with each other and with the family
# they are built from.

# --- Induced density ----------------------------------------------------------

test_that("a point mixing induces the family at its own parameter", {
  fam <- multinomial_family(n_trials = 10, k = 2)
  theta <- c(0.75, 0.25)
  p <- mixture(point_mixing(theta_star = theta), fam)
  x <- enumerate_space(fam@sample_space)
  expect_equal(dist_log_density(p, x), log_density(fam, theta, x))
})

test_that("an induced mixture is itself a probability distribution", {
  fam <- multinomial_family(n_trials = 9, k = 3)
  comp <- cbind(c(0.6, 0.3, 0.1), c(0.2, 0.2, 0.6), c(1 / 3, 1 / 3, 1 / 3))
  p <- mixture(
    finite_mixing(components = comp, weights = c(0.2, 0.5, 0.3)),
    fam
  )
  expect_equal(
    sum(exp(dist_log_density(p, enumerate_space(fam@sample_space)))),
    1,
    tolerance = 1e-12
  )
})

test_that("a finite mixture is the weighted sum of its components", {
  fam <- multinomial_family(n_trials = 6, k = 3)
  w <- c(0.4, 0.6)
  comp <- cbind(c(0.5, 0.3, 0.2), c(0.1, 0.1, 0.8))
  x <- enumerate_space(fam@sample_space)
  p <- mixture(finite_mixing(components = comp, weights = w), fam)

  manual <- log(
    w[1] *
      exp(log_density(fam, comp[, 1], x)) +
      w[2] * exp(log_density(fam, comp[, 2], x))
  )
  expect_equal(dist_log_density(p, x), manual, tolerance = 1e-12)
})

test_that("a degenerate finite mixing agrees with the point mixing", {
  fam <- multinomial_family(n_trials = 5, k = 3)
  theta <- c(0.5, 0.3, 0.2)
  x <- enumerate_space(fam@sample_space)
  a <- mixture(point_mixing(theta_star = theta), fam)
  b <- mixture(
    finite_mixing(components = matrix(theta, ncol = 1L), weights = 1),
    fam
  )
  expect_equal(dist_log_density(a, x), dist_log_density(b, x))
})

test_that("a zero-weight atom contributes nothing", {
  fam <- multinomial_family(n_trials = 5, k = 3)
  x <- enumerate_space(fam@sample_space)
  live <- c(0.5, 0.3, 0.2)
  with_dead <- mixture(
    finite_mixing(
      components = cbind(live, c(0.1, 0.1, 0.8)),
      weights = c(1, 0)
    ),
    fam
  )
  expect_equal(dist_log_density(with_dead, x), log_density(fam, live, x))
})

test_that("a mixture carrying a boundary atom stays finite where it should", {
  # P_W must be strictly positive wherever Q is, or the ratio is undefined; a
  # boundary atom in the mixture is the usual way that gets violated.
  fam <- multinomial_family(n_trials = 5, k = 3)
  p <- mixture(
    finite_mixing(
      components = cbind(c(0.5, 0.5, 0), c(1 / 3, 1 / 3, 1 / 3)),
      weights = c(0.5, 0.5)
    ),
    fam
  )
  ld <- dist_log_density(p, enumerate_space(fam@sample_space))
  expect_true(all(is.finite(ld)))
  expect_equal(sum(exp(ld)), 1, tolerance = 1e-12)
})

test_that("dist_log_density accepts a bare vector as one outcome", {
  fam <- multinomial_family(n_trials = 10, k = 2)
  p <- mixture(point_mixing(theta_star = c(0.75, 0.25)), fam)
  expect_length(dist_log_density(p, c(8, 2)), 1L)
})

# --- Sampling -----------------------------------------------------------------

test_that("dist_draw returns the right shape and respects the trial total", {
  fam <- multinomial_family(n_trials = 7, k = 3)
  p <- mixture(
    finite_mixing(
      components = cbind(c(0.6, 0.3, 0.1), c(0.1, 0.1, 0.8)),
      weights = c(0.5, 0.5)
    ),
    fam
  )
  set.seed(2)
  d <- dist_draw(p, 200L)
  expect_equal(dim(d), c(200L, 3L))
  expect_true(all(rowSums(d) == 7))
})

test_that("dist_draw handles a mixture where one component is never selected", {
  # The per-component loop must not assume every component draws at least once.
  fam <- multinomial_family(n_trials = 4, k = 2)
  p <- mixture(
    finite_mixing(
      components = cbind(c(0.5, 0.5), c(0.9, 0.1)),
      weights = c(1, 0)
    ),
    fam
  )
  set.seed(5)
  d <- dist_draw(p, 30L)
  expect_equal(dim(d), c(30L, 2L))
  expect_false(anyNA(d))
})

test_that("draws from a mixture are consistent with its density", {
  # A slow-but-honest check that induced_draw samples the law
  # induced_log_density describes: empirical frequencies converge to the pmf.
  skip_on_cran()
  fam <- multinomial_family(n_trials = 4, k = 2)
  p <- mixture(
    finite_mixing(
      components = cbind(c(0.75, 0.25), c(0.25, 0.75)),
      weights = c(0.5, 0.5)
    ),
    fam
  )
  set.seed(42)
  d <- dist_draw(p, 2e5)
  x <- enumerate_space(fam@sample_space)
  empirical <- vapply(
    seq_len(nrow(x)),
    \(i) mean(d[, 1] == x[i, 1]),
    numeric(1)
  )
  expect_equal(empirical, exp(dist_log_density(p, x)), tolerance = 0.01)
})

# --- Family recovery ----------------------------------------------------------

test_that("dist_family reports the family for a mixture and NULL otherwise", {
  fam <- multinomial_family(n_trials = 3, k = 2)
  p <- mixture(point_mixing(theta_star = c(0.5, 0.5)), fam)
  expect_identical(dist_family(p), fam)

  direct <- new_class("direct_law", parent = outcome_distribution)
  expect_null(dist_family(direct()))
})

test_that("mixture rejects components of the wrong type", {
  fam <- multinomial_family(n_trials = 3, k = 2)
  expect_error(mixture(fam, fam))
  expect_error(mixture(point_mixing(theta_star = c(0.5, 0.5)), "not a family"))
})
