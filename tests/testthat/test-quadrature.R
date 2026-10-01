# Properties of R/quadrature.R.
#
# An engine is a rule for E_Q[.]: nodes, log weights, and log q at those nodes.
# Every quantity the optimiser touches is such an expectation, so these
# invariants are the contract the whole algorithm rests on. In particular
# `sum_i w_i = 1` is what makes `sum_j w_j G(theta_j) = 1` hold downstream.
#
# A caller meets an engine as a spec handed to `ripr_init()`, which resolves it
# against the alternative and carries the resolved rule on the state it
# returns. That is where these tests read it from.

q_binomial <- function(n = 10, p = 0.75) {
  fam <- multinomial_family(n_trials = n, k = 2)
  list(
    family = fam,
    Q = mixture(fam, dirac(theta = c(p, 1 - p))),
    at = c(0.5, 0.5)
  )
}

q_gaussian <- function(d = 1, mean = NULL, sigma = NULL) {
  fam <- gaussian_family(d = d, sigma = sigma)
  if (is.null(mean)) {
    mean <- rep(0, d)
  }
  list(
    family = fam,
    Q = mixture(fam, dirac(theta = mean)),
    at = rep(0, d)
  )
}

# The rule `ripr_init()` resolves `spec` to. The null is incidental here, so it
# is a single point of the parameter space.
resolved <- function(spec, s) {
  null <- null_model(s$family, point_region(theta = s$at))
  ripr_init(s$Q, null, engine = spec)@engine
}

# E_Q[v] under a resolved rule: the weighted sum over its nodes.
e_q <- function(eng, v) sum(exp(eng@log_w) * v)

# A likelihood ratio between two parameters of the same family, which is the
# shape the oracle evaluates. Deliberately not `p_theta / q`: that integrates to
# 1 for every theta, so it would pass on a broken engine.
ratio_between <- function(eng, num, den) {
  ld <- compile_loglik(eng@family, eng@nodes)
  log_v <- as.vector(ld(matrix(num, nrow = 1L))) -
    as.vector(ld(matrix(den, nrow = 1L)))
  e_q(eng, exp(log_v))
}

# --- Exact engine -------------------------------------------------------------

test_that("the exact engine integrates against the true Q", {
  s <- q_binomial()
  eng <- resolved(exact_engine(), s)

  # A deterministic rule over the whole support, with weights summing to one.
  expect_true(eng@deterministic)
  expect_identical(
    nrow(eng@nodes),
    nrow(enumerate_space(s$family@sample_space))
  )
  expect_lt(abs(sum(exp(eng@log_w)) - 1), rounding_tol(1))

  # E_Q[X_1] for Binomial(10, 0.75) is 7.5.
  expect_lt(abs(e_q(eng, eng@nodes[, 1]) - 7.5), rounding_tol(7.5))
})

test_that("the exact engine drops nodes carrying no Q mass", {
  # A degenerate Q puts zero mass on most of the support. Those nodes must be
  # screened out, or `log_q - log_p_W` gives 0 * -Inf = NaN downstream.
  s <- q_binomial(n = 4, p = 1)
  state <- ripr_init(
    s$Q,
    null_model(s$family, point_region(theta = s$at)),
    engine = exact_engine()
  )
  eng <- state@engine

  expect_equal(nrow(eng@nodes), 1L)
  expect_equal(eng@nodes[1, ], c(4, 0))
  expect_true(all(is.finite(eng@log_q)))
  # And the fit starts from a finite KL rather than a NaN.
  expect_true(is.finite(state@trace$kl[[1L]]))
})

test_that("the exact engine refuses a family with no enumerable support", {
  expect_error(resolved(exact_engine(), q_gaussian()), "cannot be enumerated")
})

# --- Monte Carlo engine -------------------------------------------------------

test_that("a Monte Carlo engine freezes its draws at resolution", {
  s <- q_binomial()
  set.seed(4)
  eng <- resolved(mc_engine(100L), s)
  expect_false(eng@deterministic)
  expect_identical(nrow(eng@nodes), 100L)
  first <- e_q(eng, eng@nodes[, 1])
  # Consuming random numbers in between must not change what the engine holds.
  invisible(runif(10))
  expect_equal(e_q(eng, eng@nodes[, 1]), first)
})

test_that("resolving the same spec twice gives independent draws", {
  # Certification resolves a fresh engine rather than reusing the fit's nodes,
  # so a spec must not be sticky.
  s <- q_binomial()
  spec <- mc_engine(200L)
  set.seed(5)
  a <- resolved(spec, s)
  b <- resolved(spec, s)
  expect_false(identical(a@nodes, b@nodes))
})

# --- Agreement between engines ------------------------------------------------

test_that("the exact engine reproduces a closed-form likelihood ratio", {
  # For Q = Bin(n, q), E_Q[p_a / p_b] = (q a / b + (1 - q)(1 - a)/(1 - b))^n by
  # the binomial theorem.
  n <- 10
  q <- 0.75
  a <- 0.6
  b <- 0.5
  s <- q_binomial(n = n, p = q)
  eng <- resolved(exact_engine(), s)
  closed <- (q * a / b + (1 - q) * (1 - a) / (1 - b))^n

  expect_lt(
    abs(ratio_between(eng, c(a, 1 - a), c(b, 1 - b)) - closed),
    rounding_tol(closed)
  )
})

test_that("exact and Monte Carlo agree to Monte Carlo error", {
  # The justification for collapsing the engines into one interface: the same
  # expression must be computable either way.
  s <- q_binomial(n = 10, p = 0.75)
  exact <- resolved(exact_engine(), s)
  set.seed(7)
  mc <- resolved(mc_engine(2e5), s)

  expect_equal(
    e_q(mc, mc@nodes[, 1]),
    e_q(exact, exact@nodes[, 1]),
    tolerance = 0.01
  )
  expect_equal(
    ratio_between(mc, c(0.6, 0.4), c(0.5, 0.5)),
    ratio_between(exact, c(0.6, 0.4), c(0.5, 0.5)),
    tolerance = 0.02
  )
})

# --- Validation ---------------------------------------------------------------

test_that("mc_engine rejects a non-positive draw count", {
  expect_error(mc_engine(0L), "whole number larger than or equal to 1")
})

# --- Gauss-Hermite ------------------------------------------------------------

test_that("a one-node rule is the mean, not an error", {
  # With `n = 1` the Jacobi matrix is the scalar 0: the single node carries all
  # the weight.
  s <- q_gaussian(d = 1, mean = 1.5)
  eng <- resolved(gh_engine(1L), s)
  expect_identical(nrow(eng@nodes), 1L)
  expect_equal(eng@nodes[1, 1], 1.5)
  expect_lt(abs(sum(exp(eng@log_w)) - 1), rounding_tol(1))
})

test_that("a three-node rule has the known nodes and weights", {
  # The probabilists' three-point rule sits at `0, +-sqrt(3)` with weights
  # `2/3, 1/6, 1/6`, shifted and scaled by the alternative's mean and sd.
  s <- q_gaussian(d = 1, mean = 1.5, sigma = matrix(4))
  eng <- resolved(gh_engine(3L), s)
  o <- order(eng@nodes[, 1])
  expect_lt(
    max(abs(eng@nodes[o, 1] - (1.5 + 2 * c(-sqrt(3), 0, sqrt(3))))),
    rounding_tol(5)
  )
  expect_lt(
    max(abs(exp(eng@log_w[o]) - c(1 / 6, 2 / 3, 1 / 6))),
    rounding_tol(1)
  )
  expect_true(eng@deterministic)

  # An n-point rule is exact to degree 2n - 1, so three points reach the fifth
  # moment of N(1.5, 4).
  x <- eng@nodes[, 1]
  m <- 1.5
  v <- 4
  moments <- c(
    m,
    m^2 + v,
    m^3 + 3 * m * v,
    m^4 + 6 * m^2 * v + 3 * v^2,
    m^5 + 10 * m^3 * v + 15 * m * v^2
  )
  for (p in 1:5) {
    expect_lt(abs(e_q(eng, x^p) - moments[[p]]), 1e-12 * moments[[p]] * 10)
  }
})

test_that("Gauss-Hermite integrates low-order polynomials exactly", {
  # An n-point rule is exact for degree <= 2n - 1, so the first two moments of
  # the Gaussian should come back to machine precision, not to quadrature error.
  s <- q_gaussian(d = 1, mean = 1.5, sigma = matrix(4))
  eng <- resolved(gh_engine(10L), s)
  x <- eng@nodes[, 1]
  expect_lt(abs(sum(exp(eng@log_w)) - 1), rounding_tol(1))
  expect_lt(abs(e_q(eng, x) - 1.5), rounding_tol(1.5))
  expect_lt(abs(e_q(eng, x^2) - (1.5^2 + 4)), rounding_tol(1.5^2 + 4))
  expect_lt(
    abs(e_q(eng, x^3) - (1.5^3 + 3 * 1.5 * 4)),
    rounding_tol(1.5^3 + 3 * 1.5 * 4)
  )
})

test_that("Gauss-Hermite handles a correlated bivariate alternative", {
  sigma <- matrix(c(2, 0.7, 0.7, 1), 2, 2)
  s <- q_gaussian(d = 2, mean = c(0.5, -1), sigma = sigma)
  eng <- resolved(gh_engine(15L), s)
  n_nodes <- nrow(eng@nodes)

  expect_equal(n_nodes, 15L^2)
  expect_lt(abs(sum(exp(eng@log_w)) - 1), rounding_tol(1))
  w <- exp(eng@log_w)
  expect_lt(max(abs(colSums(w * eng@nodes) - c(0.5, -1))), rounding_tol(1))
  centred <- eng@nodes - rep(c(0.5, -1), each = n_nodes)
  expect_lt(
    max(abs(crossprod(centred * sqrt(w), centred * sqrt(w)) - sigma)),
    rounding_tol(2)
  )
})

test_that("Gauss-Hermite reproduces a closed-form likelihood ratio", {
  # With unit variance, p_0(x) / p_2(x) = exp(2 - 2x), so under Q = N(1, 1) the
  # expectation is exp(2) exactly.
  s <- q_gaussian(d = 1, mean = 1, sigma = matrix(1))
  gh <- resolved(gh_engine(40L), s)
  expect_equal(ratio_between(gh, 0, 2), exp(2), tolerance = 1e-6)
})

test_that("Gauss-Hermite and Monte Carlo agree on a likelihood ratio", {
  s <- q_gaussian(d = 1, mean = 1, sigma = matrix(1))
  gh <- resolved(gh_engine(40L), s)
  set.seed(11)
  mc <- resolved(mc_engine(2e5), s)
  expect_equal(
    ratio_between(mc, 0, 2),
    ratio_between(gh, 0, 2),
    tolerance = 0.05
  )
})

test_that("Gauss-Hermite refuses a non-Gaussian alternative", {
  expect_error(
    resolved(gh_engine(10L), q_binomial()),
    "needs a Gaussian alternative"
  )
})

test_that("Gauss-Hermite refuses a grid larger than max_nodes", {
  expect_error(
    resolved(gh_engine(40L, max_nodes = 100), q_gaussian(d = 2)),
    "above `max_nodes`"
  )
})

test_that("Gauss-Hermite works for a Gaussian-prior alternative", {
  fam <- gaussian_family(d = 1, sigma = matrix(1))
  s <- list(
    family = fam,
    Q = mixture(fam, gaussian_dist(mean = 0.5, cov = matrix(2))),
    at = 0
  )
  eng <- resolved(gh_engine(25L), s)
  # The induced mixture is N(0.5, 1 + 2), so its variance is 3.
  x <- eng@nodes[, 1]
  expect_lt(abs(e_q(eng, x) - 0.5), rounding_tol(0.5))
  expect_lt(abs(e_q(eng, (x - 0.5)^2) - 3), rounding_tol(3))
})

test_that("only a single-atom finite mixing measure is Gaussian", {
  # One atom is a Gaussian alternative with the family's own covariance; two
  # atoms make a Gaussian mixture, which the rule cannot take.
  fam <- gaussian_family(d = 1L, sigma = matrix(2))
  one <- resolved(
    gh_engine(10L),
    list(family = fam, Q = mixture(fam, dirac(1.5)), at = 0)
  )
  x <- one@nodes[, 1]
  expect_lt(abs(e_q(one, x) - 1.5), rounding_tol(1.5))
  expect_lt(abs(e_q(one, (x - 1.5)^2) - 2), rounding_tol(2))

  two <- finite_dist(atoms = rbind(0, 1), weights = c(0.5, 0.5))
  expect_error(
    resolved(
      gh_engine(10L),
      list(family = fam, Q = mixture(fam, two), at = 0)
    ),
    "needs a Gaussian alternative"
  )
})

test_that("gh_engine rejects a non-positive node count", {
  expect_error(gh_engine(0L), "positive")
})

# --- Printing -----------------------------------------------------------------

test_that("engines and specs print summaries, not dumps", {
  s <- q_binomial(n = 3L, p = 0.5)
  out <- paste(
    capture.output(print(resolved(exact_engine(), s))),
    collapse = "\n"
  )
  expect_match(out, "<quadrature>", fixed = TRUE)
  expect_match(out, "4 nodes, deterministic", fixed = TRUE)
  expect_no_match(out, "@")

  set.seed(1)
  mc <- resolved(mc_engine(200L), s)
  expect_match(
    format(mc),
    "quadrature: 200 nodes over multinomial_family, stochastic",
    fixed = TRUE
  )

  expect_output(
    print(exact_engine()),
    "<engine_spec> exact_engine()",
    fixed = TRUE
  )
  expect_output(
    print(mc_engine(200L)),
    "<engine_spec> mc_engine(200)",
    fixed = TRUE
  )
  expect_output(
    print(gh_engine(5L)),
    "<engine_spec> gh_engine(5)",
    fixed = TRUE
  )
  expect_identical(format(gh_engine(5L)), "gh_engine(5)")
})

# --- Internal kernels ---------------------------------------------------------

test_that("log_expect_q beats exponentiating before integrating", {
  # The log-space expectation behind the oracle's `G(theta)`. Its stability is
  # only visible publicly at audit scale, where likelihood ratios overflow a
  # double -- far too slow to fit in a unit test -- so it is pinned directly.
  # The naive route from logs -- exponentiate, integrate, take the log -- is
  # what overflows, and is what log_expect_q exists to avoid.
  eng <- resolved(exact_engine(), q_binomial())
  m <- nrow(eng@nodes)

  log_v <- rep(1000, m)
  expect_true(is.infinite(log(e_q(eng, exp(log_v)))))
  expect_lt(abs(log_expect_q(eng, log_v) - 1000), rounding_tol(1000))

  # And the same in the underflow direction, where the naive route gives -Inf.
  log_v <- rep(-1000, m)
  expect_equal(log(e_q(eng, exp(log_v))), -Inf)
  expect_lt(abs(log_expect_q(eng, log_v) + 1000), rounding_tol(1000))
})
