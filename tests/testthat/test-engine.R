# The exact engine's log-space arithmetic must reproduce a plain brute-force
# weighted sum over the enumerated support.

test_that("exact engine G(theta) matches brute-force enumeration", {
  n <- 6
  alt <- small_alt()
  fam <- multinomial_family(n, 3)
  eng <- exact_engine(fam, as_marginal(alt, fam))
  M <- nrow(support(fam))

  set.seed(10)
  for (rep in 1:3) {
    log_P <- log(as.vector(rmultinom(1, 100, rep(1 / M, M))) + 1)
    log_P <- log_P - log(sum(exp(log_P)))
    theta <- as.vector(rmultinom(1, 20, rep(1 / 3, 3))) / 20
    theta <- pmax(theta, 1e-6)
    theta <- theta / sum(theta)
    via_engine <- ripr:::expect_ratio(eng, log_density(fam, theta), log_P)
    brute <- brute_expect_ratio(alt, fam, theta, log_P)
    expect_equal(via_engine, brute, tolerance = 1e-10)
  }
})

test_that("expect() is the Q-weighted sum and entropy_q is consistent", {
  n <- 5
  alt <- small_alt()
  fam <- multinomial_family(n, 3)
  eng <- exact_engine(fam, as_marginal(alt, fam))

  f <- seq_len(ripr:::n_outcomes(eng)) + 0.5
  expect_equal(expect(eng, f), sum(eng@q_mass * f), tolerance = 1e-12)

  # H(Q) = -sum q log q over the positive-mass outcomes.
  lqm <- eng@log_q_mass[eng@finite_q]
  expect_equal(ripr:::entropy_q(eng), -sum(eng@q_mass_f * lqm), tolerance = 1e-12)
})

test_that("batched expect_ratio equals looping the scalar version", {
  n <- 6
  alt <- small_alt()
  fam <- multinomial_family(n, 3)
  eng <- exact_engine(fam, as_marginal(alt, fam))
  M <- nrow(support(fam))
  log_P <- log(rep(1 / M, M))

  theta_mat <- matrix(c(0.5, 0.3, 0.2, 0.4, 0.4, 0.2, 0.34, 0.33, 0.33), nrow = 3)
  batch <- ripr:::expect_ratio_batch(eng, ripr:::log_density_batch(fam, theta_mat), log_P)
  scalar <- vapply(
    seq_len(ncol(theta_mat)),
    function(j) ripr:::expect_ratio(eng, log_density(fam, theta_mat[, j]), log_P),
    numeric(1L)
  )
  expect_equal(batch, scalar, tolerance = 1e-12)
})
