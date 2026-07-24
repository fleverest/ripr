# Objective and certificate behaviour: EM never increases the KL, the scheduler
# converges on a small problem, and the certificate quantities are consistent.

test_that("each EM iteration does not increase the KL objective", {
  n <- 6
  alt <- small_alt()
  faces <- plurality_faces(3)
  fam <- multinomial_family(n, 3)
  prob <- ripr_problem(fam, faces, as_marginal(alt, fam))
  st <- mixture_state(prob$engine, 4L)
  for (k in seq_along(faces)) {
    ripr:::state_add_atom(st, init_point(faces[[k]], alt_mean(alt)), k)
  }
  ripr:::state_set_weights(st, rep(1 / length(faces), length(faces)))
  kl0 <- state_objective(st)$loss

  res <- ripr:::em_step(st, prob, em_iters = 10L, kl_atol = 1e-12, kl_rtol = 1e-10, kl_init = kl0)
  # KL after each EM iteration is non-increasing (up to numerical noise).
  expect_true(all(diff(res$kl_trace) <= 1e-9))
  expect_lte(res$kl, kl0 + 1e-9)
})

test_that("the scheduler converges and the certificate is consistent", {
  n <- 8
  alt <- point_mixing(theta_star = c(0.6, 0.25, 0.15))
  faces <- plurality_faces(3)
  fam <- multinomial_family(n, 3)
  prob <- ripr_problem(fam, faces, as_marginal(alt, fam))

  set.seed(31)
  res <- run_ripr(
    prob,
    init_atoms = init_on_faces(faces, alt@theta_star),
    init_atom_faces = seq_along(faces),
    fw_iters = 25, em_iters = 5, n_seeds = 120,
    gap_tol = 1e-7, verbose = FALSE
  )

  expect_true(res$converged)
  expect_lt(res$gap, 1e-6)

  cert <- res$certificate
  expect_equal(cert$gap_se, 0) # exact engine
  expect_equal(cert$gap_used, cert$gap) # no inflation for exact engines
  # Definitional identities of the certificate quantities.
  expect_equal(cert$kl_lower_bound, cert$kl - cert$gap_used)
  expect_equal(cert$growth_rate, cert$kl - log1p(cert$gap_used))
  expect_gte(cert$gap_used, -1e-4) # gap non-negative up to convergence noise
  expect_true(is.finite(cert$growth_rate))
})

test_that("the best-so-far duality gap trace is non-increasing", {
  n <- 8
  alt <- point_mixing(theta_star = c(0.55, 0.3, 0.15))
  faces <- plurality_faces(3)
  fam <- multinomial_family(n, 3)
  prob <- ripr_problem(fam, faces, as_marginal(alt, fam))

  set.seed(32)
  res <- run_ripr(
    prob,
    init_atoms = init_on_faces(faces, alt@theta_star),
    init_atom_faces = seq_along(faces),
    fw_iters = 20, em_iters = 3, n_seeds = 100,
    gap_tol = 1e-8, verbose = FALSE
  )
  best_gap <- cummin(res$history$gap)
  expect_true(all(diff(best_gap) <= 1e-12))
})

test_that("e_variable floors a negative gap so it never inflates", {
  fam <- multinomial_family(6, 2)
  Q <- point_mixing(theta_star = c(0.6, 0.4))
  Pstar <- finite_mixing(components = matrix(c(0.5, 0.5), ncol = 1), weights = 1)

  ev <- e_variable(
    numerator = as_marginal(Q, fam),
    projection = as_marginal(Pstar, fam),
    gap = -1e-6
  )
  expect_equal(ev@gap, 0)
  # With gap floored to 0, the corrected e-value equals the raw likelihood ratio.
  X <- support(fam)
  expect_equal(
    e_value(ev, X, corrected = TRUE),
    e_value(ev, X, corrected = FALSE),
    tolerance = 1e-12
  )
})

test_that("prune_mixture drops dead atoms and renormalises", {
  m <- finite_mixing(
    components = matrix(c(0.5, 0.5, 0.4, 0.6, 0.3, 0.7), nrow = 2),
    weights = c(0.7, 0.3 - 1e-9, 1e-9)
  )
  p <- prune_mixture(m, threshold = 1e-6)
  expect_equal(ncol(p@components), 2L)
  expect_equal(sum(p@weights), 1, tolerance = 1e-12)
  expect_equal(p@components, m@components[, 1:2], tolerance = 1e-12)
})

test_that("run_ripr returns an e_variable and prune_threshold trims it", {
  n <- 8
  alt <- point_mixing(theta_star = c(0.6, 0.25, 0.15))
  faces <- plurality_faces(3)
  fam <- multinomial_family(n, 3)
  prob <- ripr_problem(fam, faces, as_marginal(alt, fam))
  common <- list(
    prob,
    init_atoms = init_on_faces(faces, alt@theta_star),
    init_atom_faces = seq_along(faces),
    fw_iters = 20, em_iters = 5, n_seeds = 100, gap_tol = 1e-7, verbose = FALSE
  )
  # Same seed => identical trajectory; pruning only drops atoms afterwards.
  set.seed(33)
  res_full <- do.call(run_ripr, common)
  set.seed(33)
  res_pruned <- do.call(run_ripr, c(common, prune_threshold = 1e-6))

  ev <- res_pruned$e_variable
  expect_true(S7::S7_inherits(ev, e_variable))
  expect_true(S7::S7_inherits(ev@projection, marginal))
  expect_true(S7::S7_inherits(ev@projection@mixing, finite_mixing))
  # The projection is surfaced directly on the result too.
  expect_identical(res_pruned$projection, ev@projection)
  expect_identical(ev@numerator@mixing, alt)
  # Pruning cannot enlarge the support, and the projection weights sum to 1.
  expect_lte(ncol(ev@projection@mixing@components), length(res_full$weights))
  expect_equal(sum(ev@projection@mixing@weights), 1, tolerance = 1e-12)
  # The e-variable carries the certified gap, floored at 0.
  expect_equal(ev@gap, max(res_pruned$certificate$gap_used, 0))

  # e_value: corrected divides the raw ratio by (1 + gap); log is consistent.
  X <- support(prob$family)
  raw <- e_value(ev, X, corrected = FALSE)
  cor <- e_value(ev, X, corrected = TRUE)
  expect_equal(cor, raw / (1 + ev@gap), tolerance = 1e-10)
  expect_equal(e_value(ev, X, log = TRUE), log(cor), tolerance = 1e-10)
})
