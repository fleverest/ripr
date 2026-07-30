# Monte Carlo engine + Gaussian family: the continuous-support path. Exercises
# the common-random-numbers draw set, the stochastic standard error, and a
# short scheduler run against a half-space null.

test_that("mc_engine is stochastic with a non-zero standard error", {
  fam <- gaussian_family(2)
  alt <- gaussian_mixing(prior_mean = c(1, 0), prior_cov = diag(2))
  set.seed(1)
  eng <- mc_engine(fam, as_marginal(alt, fam), n_draws = 500L)

  expect_false(ripr:::deterministic(eng))
  expect_equal(ripr:::n_outcomes(eng), 500L)
  # set.seed() reproduces the identical draw set (common random numbers).
  set.seed(1)
  eng2 <- mc_engine(fam, as_marginal(alt, fam), n_draws = 500L)
  expect_equal(eng@outcomes, eng2@outcomes)

  f <- ripr:::eval_log_density(eng, c(0.5, 0.2))
  expect_gt(ripr:::expect_se(eng, exp(f - mean(f))), 0)
})

test_that("Gaussian half-space RIPr runs and certifies a finite gap", {
  fam <- gaussian_family(2)
  # Alternative mean sits outside the null {theta_1 <= theta_2}.
  alt <- gaussian_mixing(prior_mean = c(1, 0), prior_cov = 0.25 * diag(2))
  Q <- as_marginal(alt, fam)
  face <- halfspace_face(v = c(1, -1), c = 0, face_index = 1)
  null <- null_region(faces = list(face))
  set.seed(7)
  eng <- mc_engine(fam, Q, n_draws = 800L)
  prob <- ripr_problem(fam, null, Q, engine = eng)

  init <- init_point(face, alt@prior_mean)
  set.seed(8)
  res <- run_ripr(
    prob,
    init_atoms = matrix(init, ncol = 1L),
    init_atom_faces = 1L,
    fw_iters = 8, em_iters = 3, n_seeds = 40,
    gap_tol = 1e-3, verbose = FALSE,
    certify_ess_min = 0 # certificate quality is not what this test checks
  )

  cert <- res$certificate
  expect_true(is.finite(cert$gap))
  expect_true(is.finite(cert$kl))
  expect_equal(res$e_variable@gap, max(cert$gap_used, 0))
  expect_gte(cert$gap_se, 0)
  # Monte Carlo certificate inflates the gap: gap_used >= raw gap.
  expect_gte(cert$gap_used, cert$gap - 1e-12)
})
