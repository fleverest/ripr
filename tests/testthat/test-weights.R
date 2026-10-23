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

# KKT residual of `w` over the atoms of `ld_all`, by the engine's G.
kkt_residual <- function(ld_all, w, engine) {
  g <- ripr:::atom_g(ld_all, ripr:::mixture_log_p(ld_all, w), engine)
  on <- w > 0
  max(c(abs(g[on] - 1), g[!on] - 1, 0))
}

test_that("the SQP solve reaches SLSQP's optimum, with KKT residual near rounding", {
  # SQP and SLSQP should agree
  st <- weights_fixture()
  ld_all <- ripr:::compile_engine(st@engine)(st@mixing@atoms)
  w0 <- st@mixing@weights
  sqp <- ripr:::sqp_weights(ld_all, w0, st@engine, tol = 1e-12, max_iter = 500L)
  slsqp <- ripr:::slsqp_weights(
    ld_all,
    w0,
    st@engine,
    tol = 1e-10,
    max_iter = 500L
  )
  kl <- \(w) ripr:::kl_at(st@engine, ripr:::mixture_log_p(ld_all, w))
  expect_lte(sqp$kkt, 1e-12)
  expect_lte(kkt_residual(ld_all, sqp$weights, st@engine), 1e-11)
  expect_lte(kl(sqp$weights), kl(slsqp) + 1e-14)
  expect_equal(sum(sqp$weights), 1)
  expect_true(any(sqp$weights == 0))
  expect_true(all(sqp$weights >= 0))
})

test_that("a solve that stops short falls back to SLSQP, from where SQP stopped", {
  st <- weights_fixture()
  ld_all <- ripr:::compile_engine(st@engine)(st@mixing@atoms)
  w0 <- st@mixing@weights
  starts <- list()
  slsqp <- ripr:::slsqp_weights
  local_mocked_bindings(
    slsqp_weights = function(ld_all, w, ...) {
      starts[[length(starts) + 1L]] <<- w
      slsqp(ld_all, w, ...)
    },
    .package = "ripr"
  )
  kl <- \(w) ripr:::kl_at(st@engine, ripr:::mixture_log_p(ld_all, w))

  # SQP reaches its tolerance: no fallback.
  ripr:::solve_weights(ld_all, w0, st@engine, ripr_control())
  expect_length(starts, 0L)

  # Two iterations cannot: SLSQP runs once, from SQP's end point, and the
  # better of the two is kept.
  ctl <- ripr_control(fc_max_iter = 2L) # Inhibits SQP but not SLSQP.
  short <- ripr:::sqp_weights(
    ld_all,
    w0,
    st@engine,
    tol = ctl$kkt_tol,
    max_iter = 2L
  )
  expect_gt(short$kkt, ctl$kkt_tol)
  w <- ripr:::solve_weights(ld_all, w0, st@engine, ctl)
  expect_length(starts, 1L)
  expect_equal(starts[[1L]], short$weights)
  expect_lte(kl(w), kl(short$weights))
  expect_equal(sum(w), 1)
})

test_that("weight_step gives a near-zero support gap with SQP", {
  sqp <- weight_step(weights_fixture())
  slsqp <- weight_step(weights_fixture(
    control = ripr_control(weight_solver = "slsqp")
  ))
  expect_true(support_gap_below(1e-11)(sqp))
  # Should not have worse KL than SLSQP
  expect_lte(
    utils::tail(sqp@trace$kl, 1L),
    utils::tail(slsqp@trace$kl, 1L) + 1e-14
  )
})

test_that("ripr_control verifies the solver settings", {
  expect_identical(ripr_control()$weight_solver, "sqp")
  expect_identical(ripr_control(weight_solver = "slsqp")$weight_solver, "slsqp")
  expect_error(ripr_control(weight_solver = "other"), "must be one of")
  expect_error(ripr_control(kkt_tol = -1), "larger than or equal to 0")
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
  merged <- ripr:::merge_state(dup)
  expect_identical(n_atoms(merged@mixing), 3L)
  expect_equal(exact_kl(merged), exact_kl(dup), tolerance = 1e-14)
  solved <- weight_step(dup)
  expect_lte(n_atoms(solved@mixing), 3L)
})
