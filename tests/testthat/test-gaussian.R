# Properties of R/gaussian.R.
#
# The Gaussian is the continuous distribution that we include in this repo.
# Infinite support, and no certified gap bound available. Correctness is checked
# against closed forms rather than against enumeration.

test_that("the univariate Gaussian density matches dnorm", {
  fam <- gaussian_family(dim = 1, sigma = matrix(4))
  x <- matrix(seq(-5, 5, by = 0.5), ncol = 1L)
  expect_equal(
    kernel_loglik(fam, 1.5, x),
    stats::dnorm(x[, 1], mean = 1.5, sd = 2, log = TRUE)
  )
})

test_that("the multivariate density matches an explicit quadratic form", {
  sigma <- matrix(c(2, 0.6, 0.6, 1), 2, 2)
  fam <- gaussian_family(dim = 2, sigma = sigma)
  theta <- c(0.3, -0.7)
  x <- rbind(c(0, 0), c(1, 2), c(-1.5, 0.4))

  manual <- apply(x, 1L, function(xi) {
    dx <- xi - theta
    -0.5 *
      (2 * log(2 * pi) + log(det(sigma)) + drop(dx %*% solve(sigma) %*% dx))
  })
  expect_equal(kernel_loglik(fam, theta, x), manual)
})

test_that("compile_loglik matches the wrapper across parameter columns", {
  # The Gaussian compiles by splitting the quadratic form, not by caching a
  # constant, so this checks a different algebraic path from the multinomial.
  sigma <- matrix(c(1.5, -0.4, -0.4, 0.8), 2, 2)
  fam <- gaussian_family(dim = 2, sigma = sigma)
  set.seed(1)
  x <- matrix(rnorm(40), ncol = 2L)
  theta_mat <- cbind(c(0, 0), c(1, -1), c(2.5, 0.25))

  ld <- compile_loglik(fam, x)
  for (i in seq_len(ncol(theta_mat))) {
    expect_equal(ld(theta_mat)[, i], kernel_loglik(fam, theta_mat[, i], x))
  }
})

test_that("score matches a central finite difference", {
  sigma <- matrix(c(1.2, 0.3, 0.3, 2), 2, 2)
  fam <- gaussian_family(dim = 2, sigma = sigma)
  set.seed(2)
  x <- matrix(rnorm(20), ncol = 2L)
  theta <- c(0.4, -0.2)
  s <- score(fam, theta, x)
  eps <- 1e-6

  for (j in 1:2) {
    e <- numeric(2)
    e[j] <- eps
    fd <- (kernel_loglik(fam, theta + e, x) - kernel_loglik(fam, theta - e, x)) /
      (2 * eps)
    expect_equal(s[, j], fd, tolerance = 1e-6)
  }
})

test_that("draws have the right mean and covariance", {
  skip_on_cran()
  sigma <- matrix(c(2, 0.5, 0.5, 1), 2, 2)
  fam <- gaussian_family(dim = 2, sigma = sigma)
  theta <- c(1, -2)
  set.seed(3)
  d <- kernel_draw(fam, matrix(theta, nrow = length(theta), ncol = 2e5))
  expect_equal(dim(d), c(2e5L, 2L))
  expect_equal(colMeans(d), theta, tolerance = 0.02)
  expect_equal(stats::cov(d), sigma, tolerance = 0.05)
})

test_that("the covariance must be symmetric positive definite", {
  expect_error(
    gaussian_family(dim = 2, sigma = matrix(c(1, 2, 3, 4), 2)),
    "symmetric"
  )
  expect_error(
    gaussian_family(dim = 2, sigma = matrix(c(1, 2, 2, 1), 2)),
    "positive definite"
  )
  expect_error(gaussian_family(dim = 2, sigma = diag(3)), "2 by 2")
})

# --- Conjugate mixing measure -------------------------------------------------

test_that("a Gaussian prior induces a mixture with inflated covariance", {
  sigma <- matrix(c(1, 0.2, 0.2, 1.5), 2, 2)
  prior_cov <- matrix(c(0.5, 0, 0, 0.25), 2, 2)
  fam <- gaussian_family(dim = 2, sigma = sigma)
  m <- c(0.5, -1)
  p <- induced_distribution(fam, gaussian_mixing(prior_mean = m, prior_cov = prior_cov))

  x <- rbind(c(0, 0), c(1, 1), c(-2, 0.5))
  direct <- gaussian_family(dim = 2, sigma = sigma + prior_cov)
  expect_equal(log_density(p, x), kernel_loglik(direct, m, x))
})

test_that("draws from the induced mixture match the inflated covariance", {
  skip_on_cran()
  sigma <- diag(c(1, 2))
  prior_cov <- diag(c(0.5, 0.5))
  fam <- gaussian_family(dim = 2, sigma = sigma)
  p <- induced_distribution(fam, gaussian_mixing(prior_mean = c(1, 0), prior_cov = prior_cov))
  set.seed(4)
  d <- draw(p, 2e5)
  expect_equal(colMeans(d), c(1, 0), tolerance = 0.02)
  expect_equal(stats::cov(d), sigma + prior_cov, tolerance = 0.05)
})

test_that("a Gaussian mixing measure has no finite atom count", {
  # It is a continuous measure over the parameter space, so `n_atoms` has no
  # answer rather than the answer 1.
  g <- gaussian_mixing(prior_mean = 0, prior_cov = matrix(1))
  expect_true(is.na(n_atoms(g)))
})

test_that("the prior covariance must match the prior mean", {
  expect_error(
    gaussian_mixing(prior_mean = c(0, 0), prior_cov = matrix(1)),
    "2 by 2"
  )
})
