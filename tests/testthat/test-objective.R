# Objective and certificate behaviour: EM never increases the KL, the scheduler
# converges on a small problem, and the certificate quantities are consistent.

test_that("each EM iteration does not increase the KL objective", {
  n <- 6
  alt <- small_alt()
  faces <- plurality_faces(3)
  prob <- ripr_problem(multinomial_family(n, 3), faces, alt)
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
  alt <- point_alt(theta_star = c(0.6, 0.25, 0.15))
  faces <- plurality_faces(3)
  prob <- ripr_problem(multinomial_family(n, 3), faces, alt)

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
  alt <- point_alt(theta_star = c(0.55, 0.3, 0.15))
  faces <- plurality_faces(3)
  prob <- ripr_problem(multinomial_family(n, 3), faces, alt)

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
