# Properties of R/family.R.
#
# Nothing here compares against the old package: these are checks against the
# definition of a multinomial, plus the invariants that make compiling a
# log-likelihood safe.

# --- Sample space -------------------------------------------------------------

test_that("a family carries the sample space its outcomes live in", {
  # Enumeration, membership and dimension belong to the space; a family only
  # holds a reference to one. Their properties live in `test-sample_space.R`.
  fam <- multinomial_family(n_trials = 4, k = 3)
  expect_identical(fam@sample_space, count_space(n = 4L, k = 3L))
  expect_identical(gaussian_family(dim = 2)@sample_space, real_space(2L))
})

test_that("param_dim for multinomial is the number of categories", {
  expect_equal(param_dim(multinomial_family(n_trials = 7, k = 5)), 5L)
})

# --- Log density --------------------------------------------------------------

test_that("the multinomial pmf sums to 1 over its support", {
  fam <- multinomial_family(n_trials = 10, k = 3)
  x <- enumerate_space(fam@sample_space)
  for (theta in list(
    c(1 / 3, 1 / 3, 1 / 3),
    c(0.7, 0.2, 0.1),
    c(0.5, 0.5, 0) # a boundary parameter: log(0) = -Inf in the kernel
  )) {
    expect_equal(sum(exp(log_density(fam, theta, x))), 1, tolerance = 1e-12)
  }
})

test_that("a boundary parameter gives -Inf, not NaN, off its face", {
  fam <- multinomial_family(n_trials = 4, k = 3)
  ld <- log_density(fam, c(0.5, 0.5, 0), enumerate_space(fam@sample_space))
  expect_false(anyNA(ld))
  dead <- enumerate_space(fam@sample_space)[, 3] > 0 # positive count in the zero-probability category
  expect_true(all(ld[dead] == -Inf))
  expect_true(all(is.finite(ld[!dead])))
})

test_that("log_density_batch column c equals log_density of column c", {
  fam <- multinomial_family(n_trials = 6, k = 4)
  x <- enumerate_space(fam@sample_space)
  theta_mat <- cbind(
    c(0.25, 0.25, 0.25, 0.25),
    c(0.7, 0.1, 0.1, 0.1),
    c(0.4, 0.3, 0.2, 0.1)
  )
  batched <- log_density_batch(fam, theta_mat, x)
  for (i in seq_len(ncol(theta_mat))) {
    expect_equal(batched[, i], log_density(fam, theta_mat[, i], x))
  }
})

test_that("log_density accepts a bare vector as one outcome", {
  fam <- multinomial_family(n_trials = 10, k = 2)
  expect_length(log_density(fam, c(0.5, 0.5), c(8, 2)), 1L)
})

# --- Compiled log-likelihood --------------------------------------------------
# Compiling is a pure optimisation: it must never change an answer. These are
# what allow an engine to compile once and call the result blindly thereafter.

test_that("a compiled evaluator matches the one-off wrapper", {
  fam <- multinomial_family(n_trials = 12, k = 4)
  x <- enumerate_space(fam@sample_space)
  ld <- compile_loglik(fam, x)
  theta_mat <- cbind(
    c(0.25, 0.25, 0.25, 0.25),
    c(0.7, 0.1, 0.1, 0.1),
    c(0.5, 0.5, 0, 0) # boundary column: -Inf enters the kernel
  )
  expect_equal(ld(theta_mat), log_density_batch(fam, theta_mat, x))
})

test_that("a compiled evaluator is reusable across different theta", {
  fam <- multinomial_family(n_trials = 8, k = 3)
  x <- enumerate_space(fam@sample_space)
  ld <- compile_loglik(fam, x)
  for (theta in list(
    c(1 / 3, 1 / 3, 1 / 3),
    c(0.8, 0.1, 0.1),
    c(0.5, 0.5, 0)
  )) {
    expect_equal(
      as.vector(ld(matrix(theta, ncol = 1L))),
      log_density(fam, theta, x)
    )
  }
})

test_that("a compiled evaluator does not alias its outcome matrix", {
  # It closes over `x`; mutating the caller's copy afterwards must not change
  # what it computes. Copy-on-modify ensures this, just including as a sanity
  # check.
  fam <- multinomial_family(n_trials = 6, k = 3)
  x <- enumerate_space(fam@sample_space)
  ld <- compile_loglik(fam, x)
  before <- ld(matrix(c(0.5, 0.3, 0.2), ncol = 1L))
  x[1L, 1L] <- 99L
  expect_equal(ld(matrix(c(0.5, 0.3, 0.2), ncol = 1L)), before)
})

test_that("compiling against a subset of outcomes evaluates on that subset", {
  # Monte Carlo and quadrature engines compile against their own nodes, not the
  # full support, so this is the ordinary case rather than an edge one.
  fam <- multinomial_family(n_trials = 10, k = 3)
  nodes <- enumerate_space(fam@sample_space)[c(2L, 5L, 9L), , drop = FALSE]
  ld <- compile_loglik(fam, nodes)
  theta <- c(0.5, 0.3, 0.2)
  expect_equal(
    as.vector(ld(matrix(theta, ncol = 1L))),
    log_density(fam, theta, nodes)
  )
})

test_that("compiling works for repeated outcomes, as Monte Carlo draws give", {
  fam <- multinomial_family(n_trials = 4, k = 2)
  set.seed(7)
  nodes <- draw(fam, c(0.6, 0.4), 50L)
  ld <- compile_loglik(fam, nodes)
  expect_equal(
    as.vector(ld(matrix(c(0.5, 0.5), ncol = 1L))),
    log_density(fam, c(0.5, 0.5), nodes)
  )
})

test_that("compiling accepts a bare vector as a single outcome", {
  fam <- multinomial_family(n_trials = 10, k = 2)
  ld <- compile_loglik(fam, c(8, 2))
  expect_equal(dim(ld(matrix(c(0.5, 0.5), ncol = 1L))), c(1L, 1L))
})

test_that("a family with no compile_loglik method errors", {
  # compile_loglik is the one density method a family must supply; there is no
  # default, so an incomplete family fails loudly rather than silently.
  toy <- new_class("toy_family", parent = parametric_family)
  toy_fam <- toy(sample_space = real_space(1L))
  expect_error(compile_loglik(toy_fam, matrix(1:4, nrow = 2)))
})

# --- Score --------------------------------------------------------------------

test_that("Mnom score matches a central finite difference of log_density", {
  # Differentiated along simplex-tangent directions, since theta must stay on
  # the simplex. The score is returned in raw coordinates, so the directional
  # derivative is what lines up.
  fam <- multinomial_family(n_trials = 8, k = 3)
  x <- enumerate_space(fam@sample_space)
  theta <- c(0.45, 0.35, 0.20)
  s <- score(fam, theta, x)
  eps <- 1e-6

  for (dir in list(c(1, -1, 0), c(0, 1, -1), c(1, 0, -1))) {
    fd <- (log_density(fam, theta + eps * dir, x) -
      log_density(fam, theta - eps * dir, x)) /
      (2 * eps)
    expect_equal(as.vector(s %*% dir), fd, tolerance = 1e-5)
  }
})

test_that("Mnom score is x_j / theta_j", {
  fam <- multinomial_family(n_trials = 8, k = 3)
  theta <- c(0.45, 0.35, 0.20)
  x <- enumerate_space(fam@sample_space)
  expect_equal(score(fam, theta, x), sweep(x, 2L, theta, "/"))
})

test_that("Mnom score is 0 rather than NaN where a zero count meets zero mass", {
  fam <- multinomial_family(n_trials = 3, k = 3)
  s <- score(fam, c(0.5, 0.5, 0), enumerate_space(fam@sample_space))
  zero_count <- enumerate_space(fam@sample_space)[, 3] == 0
  expect_false(any(is.nan(s[zero_count, 3])))
  expect_true(all(s[zero_count, 3] == 0))
})

# --- Sampling -----------------------------------------------------------------

test_that("draw returns the right shape and respects the trial total", {
  n <- 12
  for (n_trials in c(1, 25, 50)) {
    for (k in 1:5) {
      fam <- multinomial_family(n_trials = n_trials, k = k)
      d <- draw(fam, seq_len(k) / sum(seq_len(k)), n)
      expect_equal(dim(d), c(n, k))
      expect_true(all(rowSums(d) == n_trials))
    }
  }
})

test_that("multinomial_family rejects malformed parameters", {
  expect_error(multinomial_family(n_trials = -1, k = 2), "non-negative")
  expect_error(multinomial_family(n_trials = 5, k = 0), "k` must be")
  expect_error(multinomial_family(n_trials = c(1, 2), k = 2), "single")
})
