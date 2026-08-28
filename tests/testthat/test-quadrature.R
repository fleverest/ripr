# Properties of R/quadrature.R.
#
# An engine is a rule for E_Q[.]: nodes, log weights, and log q at those nodes.
# Every quantity the optimiser touches is such an expectation, so these
# invariants are the contract the whole algorithm rests on. In particular
# `sum_i w_i = 1` is what makes `sum_j w_j G(theta_j) = 1` hold downstream.

q_binomial <- function(n = 10, p = 0.75) {
  fam <- multinomial_family(n_trials = n, k = 2)
  list(
    family = fam,
    Q = induced_distribution(fam, point_mixing(theta_star = c(p, 1 - p)))
  )
}

# --- The normalisation invariant ----------------------------------------------

test_that("quadrature weights sum to one for every engine", {
  s <- q_binomial()
  exact <- resolve_engine(exact_engine(), s$Q, s$family)
  set.seed(1)
  mc <- resolve_engine(mc_engine(500L), s$Q, s$family)

  expect_lt(abs(sum(exp(exact@log_w)) - 1), rounding_tol(1))
  expect_lt(abs(sum(exp(mc@log_w)) - 1), rounding_tol(1))
})

test_that("E_Q[1] = 1 under every engine", {
  # The invariant restated as an expectation: the identity everything else
  # inherits.
  s <- q_binomial()
  for (eng in list(
    resolve_engine(exact_engine(), s$Q, s$family),
    {
      set.seed(2)
      resolve_engine(mc_engine(500L), s$Q, s$family)
    }
  )) {
    expect_lt(abs(expect_q(eng, rep(1, n_nodes(eng))) - 1), rounding_tol(1))
    expect_lt(abs(log_expect_q(eng, rep(0, n_nodes(eng)))), rounding_tol(0))
  }
})

test_that("resolve_engine errors when the weights do not sum to one", {
  # A malformed rule invalidates every downstream quantity, so this must be an
  # error and not a warning.
  s <- q_binomial()
  bad <- structure(
    function(alternative, family) {
      quadrature(
        nodes = enumerate_space(family@sample_space),
        log_w = log(rep(1 / 3, nrow(enumerate_space(family@sample_space)))),
        log_q = log_density(alternative, enumerate_space(family@sample_space)),
        family = family,
        deterministic = TRUE
      )
    },
    class = "ripr_engine_spec"
  )
  expect_error(resolve_engine(bad, s$Q, s$family), "sum to 1")
})

# --- Exact engine -------------------------------------------------------------

test_that("the exact engine integrates against the true Q", {
  s <- q_binomial()
  eng <- resolve_engine(exact_engine(), s$Q, s$family)
  x <- enumerate_space(s$family@sample_space)

  # E_Q[X_1] for Binomial(10, 0.75) is 7.5.
  expect_lt(abs(expect_q(eng, x[, 1]) - 7.5), rounding_tol(7.5))
  expect_true(deterministic(eng))
  expect_equal(expect_se(eng, x[, 1]), 0)
})

test_that("the exact engine drops nodes carrying no Q mass", {
  # A degenerate Q puts zero mass on most of the support. Those nodes must be
  # screened out, or `log_q - log_p_W` gives 0 * -Inf = NaN downstream.
  fam <- multinomial_family(n_trials = 4, k = 2)
  Q <- induced_distribution(fam, point_mixing(theta_star = c(1, 0)))
  eng <- resolve_engine(exact_engine(), Q, fam)

  expect_equal(n_nodes(eng), 1L)
  expect_equal(eng@nodes[1, ], c(4, 0))
  expect_lt(abs(sum(exp(eng@log_w)) - 1), rounding_tol(1))
  expect_true(all(is.finite(eng@log_q)))
})

test_that("the exact engine refuses a family with no enumerable support", {
  toy <- new_class("toy_family", parent = parametric_family)
  direct <- new_class("direct_law", parent = distribution)
  expect_error(
    resolve_engine(
      exact_engine(),
      direct(sample_space = real_space(1L)),
      toy(
        sample_space = real_space(1L),
        parameter_space = unconstrained_region(1L)
      )
    ),
    "cannot be enumerated"
  )
})

# --- Monte Carlo engine -------------------------------------------------------

test_that("the Monte Carlo engine has uniform weights over its draws", {
  s <- q_binomial()
  set.seed(3)
  eng <- resolve_engine(mc_engine(250L), s$Q, s$family)

  expect_equal(n_nodes(eng), 250L)
  expect_equal(eng@log_w, rep(-log(250), 250L))
  expect_false(deterministic(eng))
})

test_that("a Monte Carlo engine freezes its draws at resolution", {
  s <- q_binomial()
  set.seed(4)
  eng <- resolve_engine(mc_engine(100L), s$Q, s$family)
  first <- expect_q(eng, eng@nodes[, 1])
  # Consuming random numbers in between must not change what the engine holds.
  invisible(runif(10))
  expect_equal(expect_q(eng, eng@nodes[, 1]), first)
})

test_that("resolving the same spec twice gives independent draws", {
  # Certification resolves a fresh engine rather than reusing the fit's nodes,
  # so a spec must not be sticky.
  s <- q_binomial()
  spec <- mc_engine(200L)
  set.seed(5)
  a <- resolve_engine(spec, s$Q, s$family)
  b <- resolve_engine(spec, s$Q, s$family)
  expect_false(identical(a@nodes, b@nodes))
})

test_that("expect_se is the standard error of the mean under Monte Carlo", {
  s <- q_binomial()
  set.seed(6)
  eng <- resolve_engine(mc_engine(400L), s$Q, s$family)
  v <- eng@nodes[, 1]
  expect_lt(
    abs(expect_se(eng, v) - sd(v) / sqrt(400)),
    rounding_tol(sd(v) / sqrt(400))
  )
})

# --- Agreement between engines ------------------------------------------------

# A likelihood ratio between two parameters of the same family, which is the
# shape the oracle evaluates. Deliberately not `p_theta / q`: that integrates to
# 1 for every theta, so it would pass on a broken engine.
ratio_between <- function(eng, num, den) {
  ld <- compile_engine(eng)
  exp(log_expect_q(
    eng,
    as.vector(ld(matrix(num, ncol = 1L))) -
      as.vector(ld(matrix(den, ncol = 1L)))
  ))
}

test_that("the exact engine reproduces a closed-form likelihood ratio", {
  # For Q = Bin(n, q), E_Q[p_a / p_b] = (q a / b + (1 - q)(1 - a)/(1 - b))^n by
  # the binomial theorem.
  n <- 10
  q <- 0.75
  a <- 0.6
  b <- 0.5
  s <- q_binomial(n = n, p = q)
  eng <- resolve_engine(exact_engine(), s$Q, s$family)
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
  exact <- resolve_engine(exact_engine(), s$Q, s$family)
  set.seed(7)
  mc <- resolve_engine(mc_engine(2e5), s$Q, s$family)

  expect_equal(
    expect_q(mc, mc@nodes[, 1]),
    expect_q(exact, exact@nodes[, 1]),
    tolerance = 0.01
  )
  expect_equal(
    ratio_between(mc, c(0.6, 0.4), c(0.5, 0.5)),
    ratio_between(exact, c(0.6, 0.4), c(0.5, 0.5)),
    tolerance = 0.02
  )
})

test_that("log_expect_q beats exponentiating before integrating", {
  # The two take different arguments: expect_q takes values, log_expect_q takes
  # their logs. The naive route from logs -- exponentiate, integrate, take the
  # log -- is what overflows, and is what log_expect_q exists to avoid.
  s <- q_binomial()
  eng <- resolve_engine(exact_engine(), s$Q, s$family)

  log_v <- rep(1000, n_nodes(eng))
  expect_true(is.infinite(log(expect_q(eng, exp(log_v)))))
  expect_lt(abs(log_expect_q(eng, log_v) - 1000), rounding_tol(1000))

  # And the same in the underflow direction, where the naive route gives -Inf.
  log_v <- rep(-1000, n_nodes(eng))
  expect_equal(log(expect_q(eng, exp(log_v))), -Inf)
  expect_lt(abs(log_expect_q(eng, log_v) + 1000), rounding_tol(1000))
})

test_that("log_expect_q is -Inf when the integrand is everywhere zero", {
  s <- q_binomial()
  eng <- resolve_engine(exact_engine(), s$Q, s$family)
  expect_equal(log_expect_q(eng, rep(-Inf, n_nodes(eng))), -Inf)
})

# --- Compilation --------------------------------------------------------------

test_that("compile_engine evaluates the family at the engine's own nodes", {
  s <- q_binomial()
  eng <- resolve_engine(exact_engine(), s$Q, s$family)
  ld <- compile_engine(eng)
  theta_mat <- cbind(c(0.5, 0.5), c(0.25, 0.75))
  expect_equal(
    ld(theta_mat),
    kernel_loglik_batch(s$family, theta_mat, eng@nodes)
  )
})

test_that("compile_engine is consistent across repeated compilations", {
  s <- q_binomial()
  eng <- resolve_engine(exact_engine(), s$Q, s$family)
  theta <- matrix(c(0.5, 0.5), ncol = 1L)
  expect_equal(compile_engine(eng)(theta), compile_engine(eng)(theta))
})

# --- Validation ---------------------------------------------------------------

test_that("quadrature rejects mismatched node, weight and q lengths", {
  fam <- multinomial_family(n_trials = 3, k = 2)
  x <- enumerate_space(fam@sample_space)
  expect_error(
    quadrature(
      nodes = x,
      log_w = log(rep(1 / nrow(x), nrow(x) - 1L)),
      log_q = rep(0, nrow(x)),
      family = fam,
      deterministic = TRUE
    ),
    "one weight per node"
  )
  expect_error(
    quadrature(
      nodes = x,
      log_w = log(rep(1 / nrow(x), nrow(x))),
      log_q = rep(0, nrow(x) - 1L),
      family = fam,
      deterministic = TRUE
    ),
    "one log_q per node"
  )
})

test_that("mc_engine rejects a non-positive draw count", {
  expect_error(mc_engine(0L), "positive")
  expect_error(mc_engine(-5L), "positive")
})

# --- Gauss-Hermite ------------------------------------------------------------

q_gaussian <- function(d = 1, mean = NULL, sigma = NULL) {
  fam <- gaussian_family(dim = d, sigma = sigma)
  if (is.null(mean)) {
    mean <- rep(0, d)
  }
  list(
    family = fam,
    Q = induced_distribution(fam, point_mixing(theta_star = mean))
  )
}

test_that("Gauss-Hermite weights sum to one", {
  s <- q_gaussian(d = 1)
  eng <- resolve_engine(gh_engine(20L), s$Q, s$family)
  expect_lt(abs(sum(exp(eng@log_w)) - 1), rounding_tol(1))
  expect_true(deterministic(eng))
  expect_equal(n_nodes(eng), 20L)
})


test_that("a one-node rule is the mean, not an error", {
  # `n = 1` skips the eigendecomposition entirely: the single node carries all
  # the weight.
  s <- q_gaussian(d = 1, mean = 1.5)
  eng <- resolve_engine(gh_engine(1L), s$Q, s$family)
  expect_equal(n_nodes(eng), 1L)
  expect_lt(abs(sum(exp(eng@log_w)) - 1), rounding_tol(1))
})

test_that("Gauss-Hermite integrates low-order polynomials exactly", {
  # An n-point rule is exact for degree <= 2n - 1, so the first two moments of
  # the Gaussian should come back to machine precision, not to quadrature error.
  s <- q_gaussian(d = 1, mean = 1.5, sigma = matrix(4))
  eng <- resolve_engine(gh_engine(10L), s$Q, s$family)
  x <- eng@nodes[, 1]
  expect_lt(abs(expect_q(eng, x) - 1.5), rounding_tol(1.5))
  expect_lt(abs(expect_q(eng, x^2) - (1.5^2 + 4)), rounding_tol(1.5^2 + 4))
  expect_lt(
    abs(expect_q(eng, x^3) - (1.5^3 + 3 * 1.5 * 4)),
    rounding_tol(1.5^3 + 3 * 1.5 * 4)
  )
})

test_that("Gauss-Hermite handles a correlated bivariate alternative", {
  sigma <- matrix(c(2, 0.7, 0.7, 1), 2, 2)
  s <- q_gaussian(d = 2, mean = c(0.5, -1), sigma = sigma)
  eng <- resolve_engine(gh_engine(15L), s$Q, s$family)

  expect_equal(n_nodes(eng), 15L^2)
  expect_lt(abs(sum(exp(eng@log_w)) - 1), rounding_tol(1))
  w <- exp(eng@log_w)
  expect_lt(max(abs(colSums(w * eng@nodes) - c(0.5, -1))), rounding_tol(1))
  centred <- eng@nodes - rep(c(0.5, -1), each = n_nodes(eng))
  expect_lt(
    max(abs(crossprod(centred * sqrt(w), centred * sqrt(w)) - sigma)),
    rounding_tol(2)
  )
})

test_that("Gauss-Hermite reproduces a closed-form likelihood ratio", {
  # With unit variance, p_0(x) / p_2(x) = exp(2 - 2x), so under Q = N(1, 1) the
  # expectation is exp(2) exactly.
  s <- q_gaussian(d = 1, mean = 1, sigma = matrix(1))
  gh <- resolve_engine(gh_engine(40L), s$Q, s$family)
  expect_equal(ratio_between(gh, 0, 2), exp(2), tolerance = 1e-6)
})

test_that("Gauss-Hermite and Monte Carlo agree on a likelihood ratio", {
  s <- q_gaussian(d = 1, mean = 1, sigma = matrix(1))
  gh <- resolve_engine(gh_engine(40L), s$Q, s$family)
  set.seed(11)
  mc <- resolve_engine(mc_engine(2e5), s$Q, s$family)
  expect_equal(
    ratio_between(mc, 0, 2),
    ratio_between(gh, 0, 2),
    tolerance = 0.05
  )
})

test_that("Gauss-Hermite refuses a non-Gaussian alternative", {
  s <- q_binomial()
  expect_error(
    resolve_engine(gh_engine(10L), s$Q, s$family),
    "needs a Gaussian alternative"
  )
})

test_that("Gauss-Hermite refuses a grid larger than max_nodes", {
  s <- q_gaussian(d = 2)
  expect_error(
    resolve_engine(gh_engine(40L, max_nodes = 100), s$Q, s$family),
    "above `max_nodes`"
  )
})

test_that("Gauss-Hermite works for a Gaussian-prior alternative", {
  fam <- gaussian_family(dim = 1, sigma = matrix(1))
  Q <- induced_distribution(
    fam,
    gaussian_mixing(prior_mean = 0.5, prior_cov = matrix(2))
  )
  eng <- resolve_engine(gh_engine(25L), Q, fam)
  # The induced mixture is N(0.5, 1 + 2), so its variance is 3.
  x <- eng@nodes[, 1]
  expect_lt(abs(expect_q(eng, x) - 0.5), rounding_tol(0.5))
  expect_lt(abs(expect_q(eng, (x - 0.5)^2) - 3), rounding_tol(3))
})

test_that("gauss_hermite reproduces known nodes and weights", {
  gh <- gauss_hermite(3L)
  expect_lt(
    max(abs(gh$nodes - c(-sqrt(3 / 2), 0, sqrt(3 / 2)))),
    rounding_tol(sqrt(3 / 2))
  )
  expect_lt(
    max(abs(gh$weights - sqrt(pi) * c(1 / 6, 2 / 3, 1 / 6))),
    rounding_tol(sqrt(pi))
  )
  expect_lt(
    abs(sum(gauss_hermite(12L)$weights) - sqrt(pi)),
    rounding_tol(sqrt(pi))
  )
})

test_that("gh_engine rejects a non-positive node count", {
  expect_error(gh_engine(0L), "positive")
})
