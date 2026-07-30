# The shape of `run_ripr()`'s return. The contract is that each element has a
# single provenance -- deliverable, certification sample, or fit sample -- and
# that nothing straddles two of them.

cone_fit <- function(prune_threshold = 0, checkpoint_iters = NULL, ...) {
  fam <- gaussian_family(dim = 2)
  mu <- c(2, 1)
  Q <- as_marginal(point_mixing(theta_star = mu), fam)
  faces <- list(
    halfspace_face(v = c(1, -1), c = 0, face_index = 1),
    halfspace_face(v = c(1, 1), c = 0, face_index = 2)
  )
  eng <- mc_engine(fam, Q, n_draws = 1000L)
  prob <- ripr_problem(fam, null_region(faces = faces), Q, engine = eng)
  set.seed(2)
  res <- run_ripr(
    prob,
    init_atoms = cbind(init_point(faces[[1]], mu), init_point(faces[[2]], mu)),
    init_atom_faces = c(1, 2),
    fw_iters = 10, em_iters = 3, n_seeds = 30, gap_tol = 1e-9,
    prune_threshold = prune_threshold,
    checkpoint_iters = checkpoint_iters,
    verbose = FALSE,
    # These tests assert the shape of the result, not the quality of the
    # certificate; the ESS warning is unrelated noise here.
    certify_ess_min = 0,
    ...
  )
  list(res = res, prob = prob)
}

test_that("the return value is exactly the six documented elements", {
  res <- cone_fit()$res
  expect_setequal(
    names(res),
    c("projection", "e_variable", "certificate", "history", "checkpoints",
      "converged")
  )
  # The fields folded into `certificate` / `checkpoints` are gone from the top.
  for (dead in c("atoms", "weights", "atom_face_idx", "gap", "gap_raw",
                 "gap_used", "kl", "kl_ulb", "oracle_theta", "metrics",
                 "kl_trace", "state")) {
    expect_false(dead %in% names(res), info = dead)
  }
})

test_that("the final checkpoint is post-prune and aligned with the projection", {
  fit <- cone_fit(prune_threshold = 0.01, checkpoint_iters = c(0, 3))
  res <- fit$res
  fin <- res$checkpoints$final

  expect_true("final" %in% names(res$checkpoints))
  expect_identical(fin$iter, NA_integer_)
  # No fit-sample sweep was run on the pruned mixture.
  expect_true(is.na(fin$oracle_theta))

  n <- ncol(res$projection@mixing@components)
  expect_length(fin$atoms, n)
  expect_length(fin$weights, n)
  expect_length(fin$atom_face_idx, n)
  expect_equal(fin$weights, res$projection@mixing@weights)
  expect_equal(do.call(cbind, fin$atoms), res$projection@mixing@components)

  # Iteration checkpoints are pre-prune, so they may be wider -- and are named
  # so they cannot collide with `final`.
  expect_true(all(c("iter_0", "iter_3") %in% names(res$checkpoints)))
  expect_gte(length(res$checkpoints$iter_3$atoms), n)
})

test_that("the final checkpoint round-trips back into run_ripr()", {
  fit <- cone_fit(prune_threshold = 0.01)
  fin <- fit$res$checkpoints$final
  expect_no_error(
    run_ripr(
      fit$prob,
      init_atoms = do.call(cbind, fin$atoms),
      init_atom_faces = fin$atom_face_idx,
      init_weights = fin$weights,
      fw_iters = 2, em_iters = 2, n_seeds = 20, verbose = FALSE,
      certify_ess_min = 0
    )
  )
})

test_that("resuming with init_weights starts from the previous fit, not uniform", {
  fit <- cone_fit()
  fin <- fit$res$checkpoints$final
  kl_reached <- tail(fit$res$history$kl_after_em, 1)

  resume <- function(w) {
    set.seed(77)
    run_ripr(
      fit$prob,
      init_atoms = do.call(cbind, fin$atoms),
      init_atom_faces = fin$atom_face_idx,
      init_weights = w,
      fw_iters = 0L, em_iters = 0L, n_seeds = 20, verbose = FALSE,
      certify_ess_min = 0
    )
  }
  # With no FW and no EM, the starting KL is whatever the weights imply.
  kept <- resume(fin$weights)$history$kl_trace[[1]]$kl[1]
  unif <- resume(NULL)$history$kl_trace[[1]]$kl[1]

  # Carrying the weights over reproduces the KL the previous run reached ...
  expect_equal(kept, kl_reached, tolerance = 1e-8)
  # ... and is strictly better than the uniform restart it replaces.
  expect_lt(kept, unif)
})

test_that("init_weights is validated at the boundary", {
  fit <- cone_fit()
  fin <- fit$res$checkpoints$final
  A <- do.call(cbind, fin$atoms)
  f <- fin$atom_face_idx
  bad <- function(w) {
    run_ripr(fit$prob, A, f, fw_iters = 0L, em_iters = 0L,
             init_weights = w, verbose = FALSE, certify_ess_min = 0)
  }
  n <- length(f)
  expect_error(bad(rep(1 / n, n - 1L)), "must equal ncol")
  expect_error(bad(c(-0.5, rep(0.5, n - 1L))), "non-negative")
  expect_error(bad(rep(0, n)), "positive sum")
  expect_error(bad(c(NA_real_, rep(0.5, n - 1L))), "finite")
  # Unnormalised weights are accepted and renormalised.
  expect_no_error(bad(rep(2, n)))
})

test_that("`final` is present even with no checkpoints requested and fw_iters = 0", {
  fam <- gaussian_family(dim = 2)
  mu <- c(2, 1)
  Q <- as_marginal(point_mixing(theta_star = mu), fam)
  face <- halfspace_face(v = c(1, -1), c = 0, face_index = 1)
  eng <- mc_engine(fam, Q, n_draws = 500L)
  prob <- ripr_problem(fam, null_region(faces = list(face)), Q, engine = eng)
  set.seed(4)
  res <- run_ripr(
    prob,
    init_atoms = matrix(init_point(face, mu), ncol = 1L),
    init_atom_faces = 1L,
    fw_iters = 0L, em_iters = 5L, n_seeds = 20, verbose = FALSE,
    certify_ess_min = 0
  )
  expect_identical(names(res$checkpoints), "final")
  expect_length(res$checkpoints$final$atom_face_idx, 1L)
})

test_that("history is one row per outer iteration with nested list columns", {
  res <- cone_fit()$res
  h <- res$history

  expect_s3_class(h, "data.frame")
  expect_identical(nrow(h), length(unique(h$iter)))
  for (col in c("gap", "gap_se", "support_size", "kl_after_em", "kl_ulb",
                "face_idx", "elapsed_s", "oracle_theta", "kl_trace")) {
    expect_true(col %in% names(h), info = col)
  }

  # oracle_theta: one length-d parameter vector per iteration.
  expect_type(h$oracle_theta, "list")
  expect_true(all(lengths(h$oracle_theta) == 2L))

  # kl_trace: the inner init/FW/EM steps belonging to that iteration.
  expect_type(h$kl_trace, "list")
  expect_s3_class(h$kl_trace[[1]], "data.frame")
  expect_setequal(names(h$kl_trace[[1]]), c("step_type", "n_atoms", "kl"))
  expect_true("init" %in% h$kl_trace[[1]]$step_type)
  expect_true(all(vapply(h$kl_trace[-1], function(d) "fw" %in% d$step_type,
                         logical(1L))))
})

test_that("the certificate records the settings that produced it", {
  fit <- cone_fit(certify_draws = 4000L)
  cert <- fit$res$certificate

  expect_identical(cert$n_draws, 4000L)
  # Certification inherits the run's n_seeds, scaled by 10 like the draws.
  expect_identical(cert$n_seeds, 300L)

  # Default certify_draws is 10x the engine's fit sample.
  cert2 <- cone_fit()$res$certificate
  expect_identical(cert2$n_draws, 10000L)
})

test_that("exact engines report NA draws and a zero standard error", {
  fam <- multinomial_family(n_trials = 8, k = 2)
  Q <- as_marginal(point_mixing(theta_star = c(0.75, 0.25)), fam)
  face <- polytope_face(vertices = cbind(c(0, 1), c(0.5, 0.5)), face_index = 1)
  prob <- ripr_problem(fam, null_region(faces = list(face)), Q)
  set.seed(5)
  res <- run_ripr(
    prob,
    init_atoms = matrix(c(0.25, 0.75), ncol = 1L),
    init_atom_faces = 1L,
    fw_iters = 5L, em_iters = 3L, n_seeds = 40, verbose = FALSE
  )
  expect_identical(res$certificate$n_draws, NA_integer_)
  expect_identical(res$certificate$n_seeds, 400L)
  expect_equal(res$certificate$gap_se, 0)
  expect_true("final" %in% names(res$checkpoints))
})
