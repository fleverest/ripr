# Properties of the weight solvers and of atom merging.
#
# The exact-Hessian SQP must reach the same optimum as SLSQP, but to a KKT
# residual near rounding (G = 1 on the support, G <= 1 off it) with exact
# zeros, and a solve that stops short must fall back to SLSQP. Merging must
# remove duplicate atoms without changing the mixture.

weights_fixture <- function(
  n = 10L,
  q = c(0.4, 0.35, 0.25),
  control = ripr_control()
) {
  fam <- multinomial_family(n_trials = n, k = 3L)
  parts <- lapply(2:3, function(j) {
    vertices <- diag(3)
    vertices[1L, ] <- replace(numeric(3), c(1L, j), 0.5)
    simplex_region(vertices = vertices)
  })
  set.seed(1)
  atoms <- lapply(parts, function(p) {
    ch <- chart(p)
    ch$to_theta_batch(ch$seed(15L))
  })
  ripr_init(
    fam(q),
    null_model(fam, parts),
    exact_engine(),
    atoms = atoms,
    control = control
  )
}

test_that("the SQP solve reaches SLSQP's optimum, with KKT residual near rounding", {
  # SQP and SLSQP should agree
  st <- weights_fixture()
  ld_all <- compile_engine(st@engine)(st@mixing@atoms)
  w0 <- st@mixing@weights
  sqp <- sqp_weights(ld_all, w0, st@engine, tol = 1e-12, max_iter = 500L)
  slsqp <- slsqp_weights(
    ld_all,
    w0,
    st@engine,
    tol = 1e-10,
    max_iter = 500L
  )
  expect_true(slsqp$converged)
  kl <- \(w) kl_at(st@engine, mixture_log_p(ld_all, w))
  expect_lte(sqp$kkt, 1e-12)
  expect_lte(weight_kkt(ld_all, sqp$weights, st@engine), 1e-11)
  expect_lte(kl(sqp$weights), kl(slsqp$weights) + 1e-14)
  expect_equal(sum(sqp$weights), 1)
  expect_true(any(sqp$weights == 0))
  expect_true(all(sqp$weights >= 0))
})

test_that("a solve that stops short falls back to SLSQP, from where SQP stopped", {
  st <- weights_fixture()
  ld_all <- compile_engine(st@engine)(st@mixing@atoms)
  w0 <- st@mixing@weights
  starts <- list()
  slsqp <- slsqp_weights
  local_mocked_bindings(
    slsqp_weights = function(ld_all, w, ...) {
      starts[[length(starts) + 1L]] <<- w
      slsqp(ld_all, w, ...)
    },
    .package = "ripr"
  )
  kl <- \(w) kl_at(st@engine, mixture_log_p(ld_all, w))

  # SQP reaches its tolerance: no fallback.
  solve_weights(ld_all, w0, st@engine, ripr_control())
  expect_length(starts, 0L)

  # Two iterations cannot: SLSQP runs once, from SQP's end point, and the
  # better of the two is kept.
  ctl <- ripr_control(weight_max_iter = 2L) # Inhibits SQP but not SLSQP.
  short <- sqp_weights(
    ld_all,
    w0,
    st@engine,
    tol = ctl$weight_sqp_tol,
    max_iter = 2L
  )
  expect_gt(short$kkt, ctl$weight_sqp_tol)
  w <- solve_weights(ld_all, w0, st@engine, ctl)$weights
  expect_length(starts, 1L)
  expect_equal(starts[[1L]], short$weights)
  expect_lte(kl(w), kl(short$weights))
  expect_equal(sum(w), 1)
})

test_that("weight_step merges duplicate atoms before solving, leaving the mixture unchanged", {
  st <- weights_fixture()
  a <- st@mixing@atoms
  dup <- ripr_init(
    st@alternative,
    st@null,
    exact_engine(),
    atoms = list(rbind(a[1L, ], a[1L, ], a[2L, ]), a[16L, , drop = FALSE]),
    weights = list(c(0.2, 0.3, 0.25), 0.25)
  )
  merged <- merge_state(dup)
  expect_identical(n_atoms(merged@mixing), 3L)
  expect_equal(exact_kl(merged), exact_kl(dup), tolerance = 1e-14)
  solved <- weight_step(dup)
  expect_lte(n_atoms(solved@mixing), 3L)
})

# --- Solve diagnostics --------------------------------------------------------

test_that("each corrective step records its weight solve's KKT residual", {
  # `solve_kkt` is at most `weight_sqp_tol` whenever the solve converged, and
  # `solve_converged` says so; both are `NA` on steps that solve nothing. A slide's re-solve is exact too: a KL
  # that rises only by rounding must not make it keep the slide's weights.
  tol <- ripr_control()$weight_sqp_tol
  set.seed(1)
  st <- weights_fixture()
  plain <- fw_step(st, 1L)
  expect_true(is.na(utils::tail(plain@trace$solve_kkt, 1L)))
  expect_true(is.na(utils::tail(plain@trace$solve_converged, 1L)))
  solved <- utils::tail(weight_step(st)@trace, 1L)
  expect_lte(solved$solve_kkt, tol)
  expect_true(solved$solve_converged)
  corrected <- utils::tail(fw_step(st, 2L, correct_weights = TRUE)@trace, 2L)
  expect_true(all(corrected$solve_kkt <= tol))
  expect_true(all(corrected$solve_converged))
  slid <- sfw_step(sfw_step(fw_step(st, 2L)), 2L)
  rows <- slid@trace[slid@trace$phase == "sfw", ]
  expect_identical(nrow(rows), 3L)
  expect_true(all(rows$solve_kkt <= tol))
  expect_true(all(rows$solve_converged))
})

test_that("a weight solve that stops short warns, once per call", {
  # One iteration cannot reach the tolerance, nor one evaluation let SLSQP
  # converge.
  capped <- weights_fixture(control = ripr_control(weight_max_iter = 1L))
  expect_warning(weight_step(capped), class = "ripr_weight_solve_warning")
  expect_warning(
    ripr_finish(capped, correct_weights = TRUE),
    class = "ripr_weight_solve_warning"
  )
  w <- testthat::capture_warnings(weight_step(capped, times = 3L))
  expect_length(w, 1L)
  expect_match(w, "3 of 3 weight solves stopped short", fixed = TRUE)
  slsqp <- weights_fixture(
    control = ripr_control(weight_solver = "slsqp", weight_max_iter = 1L)
  )
  expect_warning(
    short <- weight_step(slsqp),
    "stopped before SLSQP converged",
    class = "ripr_weight_solve_warning"
  )
  expect_false(utils::tail(short@trace$solve_converged, 1L))
  # Within budget, SLSQP converges and says nothing.
  roomy <- weights_fixture(control = ripr_control(weight_solver = "slsqp"))
  expect_no_warning(weight_step(roomy))
})
