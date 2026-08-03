# Properties of R/scheduler.R.
#
# The identity `sum_j w_j G(theta_j) = 1` is the workhorse here. It is algebraic,
# not a convergence condition, so it holds at every iterate and fails only if the
# engine, the weights, or the atom indexing is wrong.

binomial_problem <- function(n = 10, p = 0.75) {
  fam <- multinomial_family(n_trials = n, k = 2)
  Q <- mixture(point_mixing(theta_star = c(p, 1 - p)), fam)
  null <- null_model(
    fam,
    list(
      simplex_null(vertices = cbind(c(0, 1), c(0.5, 0.5)))
    )
  )
  list(family = fam, Q = Q, null = null)
}

plurality_problem <- function(n = 12, k = 4, q = c(0.4, 0.3, 0.2, 0.1)) {
  fam <- multinomial_family(n_trials = n, k = k)
  Q <- mixture(point_mixing(theta_star = q), fam)
  subnulls <- lapply(2:k, function(j) {
    basis <- lapply(setdiff(seq_len(k), 1L), function(i) {
      v <- numeric(k)
      v[i] <- 1
      v
    })
    tie <- numeric(k)
    tie[c(1L, j)] <- 0.5
    simplex_null(vertices = do.call(cbind, c(basis, list(tie))))
  })
  list(family = fam, Q = Q, null = null_model(fam, subnulls))
}

# The identity every iterate must satisfy.
weighted_gain <- function(state) {
  ld <- compile_engine(state@engine)
  log_p <- row_logsumexp(add_by_col(
    ld(flat_atoms(state)),
    log(flat_weights(state))
  ))
  g <- exp(col_logsumexp(ld(flat_atoms(state)) - log_p + state@engine@log_w))
  sum(flat_weights(state) * g)
}

# --- Initialisation -----------------------------------------------------------

test_that("ripr_init places one atom per subnull by default", {
  p <- plurality_problem()
  st <- ripr_init(p$Q, p$null, exact_engine())
  expect_length(st@atoms, 3L)
  expect_true(all(vapply(st@atoms, ncol, integer(1)) == 1L))
  expect_equal(sum(flat_weights(st)), 1)
  expect_equal(st@iter, 0)
})

test_that("initial atoms lie in their own subnulls", {
  p <- plurality_problem()
  st <- ripr_init(p$Q, p$null, exact_engine())
  for (i in seq_along(st@atoms)) {
    expect_true(contains(p$null@subnulls[[i]], st@atoms[[i]][, 1L], tol = 1e-6))
  }
})

test_that("ripr_init rejects an atom list of the wrong length", {
  p <- plurality_problem()
  expect_error(
    ripr_init(p$Q, p$null, exact_engine(), atoms = list(matrix(rep(0.25, 4)))),
    "one element per subnull"
  )
})

test_that("ripr_init accepts empty subnulls as zero-column matrices", {
  p <- plurality_problem()
  atoms <- list(
    matrix(c(0.1, 0.4, 0.3, 0.2), ncol = 1L),
    matrix(numeric(0), nrow = 4L, ncol = 0L),
    matrix(numeric(0), nrow = 4L, ncol = 0L)
  )
  st <- ripr_init(p$Q, p$null, exact_engine(), atoms = atoms)
  expect_equal(vapply(st@atoms, ncol, integer(1)), c(1L, 0L, 0L))
  expect_equal(sum(flat_weights(st)), 1)
})

test_that("ripr_init errors when no subnull carries an atom", {
  p <- plurality_problem()
  empty <- replicate(
    3L,
    matrix(numeric(0), nrow = 4L, ncol = 0L),
    simplify = FALSE
  )
  expect_error(
    ripr_init(p$Q, p$null, exact_engine(), atoms = empty),
    "at least one subnull"
  )
})

# --- The algebraic identity ---------------------------------------------------

test_that("weighted gains sum to 1 at every iterate", {
  # Algebra, not convergence: sum_j w_j E_Q[p_j / P_W] = E_Q[P_W / P_W] = 1.
  # Fails if quadrature weights are unnormalised, or if weights and atoms are
  # misaligned across subnulls.
  p <- plurality_problem()
  set.seed(1)
  st <- ripr_init(
    p$Q,
    p$null,
    exact_engine(),
    control = ripr_control(n_seeds = 30L)
  )
  expect_equal(weighted_gain(st), 1, tolerance = 1e-10)
  for (i in 1:3) {
    st <- ripr_step(st, 1L)
    expect_equal(weighted_gain(st), 1, tolerance = 1e-10)
  }
})

test_that("the reported gap is non-negative", {
  # Follows from the identity: a weighted average is at most the maximum, so
  # sup G >= 1 whenever the atoms lie in the null.
  p <- plurality_problem()
  set.seed(2)
  st <- ripr_init(
    p$Q,
    p$null,
    exact_engine(),
    control = ripr_control(n_seeds = 30L)
  )
  st <- ripr_step(st, 3L)
  gaps <- st@trace$gap[!is.na(st@trace$gap)]
  expect_true(length(gaps) > 0)
  expect_true(all(gaps > -1e-8))
})

# --- Monotonicity -------------------------------------------------------------

test_that("KL never increases under a line search", {
  # Guaranteed because gamma = 0 is in the search interval: the step can always
  # decline to move. Not a property of Frank-Wolfe in general -- see the fixed
  # schedule below.
  p <- plurality_problem()
  set.seed(3)
  st <- ripr_init(
    p$Q,
    p$null,
    exact_engine(),
    control = ripr_control(n_seeds = 30L)
  )
  st <- ripr_step(st, 4L)
  kl <- st@trace$kl
  expect_true(all(diff(kl) <= 1e-9))
})

test_that("the fixed schedule still reduces KL over a run", {
  # The open-loop schedule consults nothing, so an individual step may overshoot
  # and increase KL; only the O(1/k) rate survives. Asserting the weaker
  # end-to-end claim keeps the distinction visible rather than pretending the
  # per-iteration one holds.
  p <- plurality_problem()
  set.seed(35)
  ctl <- ripr_control(fw_step = "fixed", n_seeds = 30L, gap_tol = -Inf)
  st <- ripr_step(ripr_init(p$Q, p$null, exact_engine(), control = ctl), 6L)
  kl <- st@trace$kl
  expect_lt(utils::tail(kl, 1L), kl[1L])
})

test_that("an EM sweep does not increase KL", {
  p <- plurality_problem()
  set.seed(4)
  st <- ripr_init(
    p$Q,
    p$null,
    exact_engine(),
    control = ripr_control(n_seeds = 30L)
  )
  st <- ripr_step(st, 1L)
  em <- st@trace[st@trace$phase == "em", ]
  expect_true(nrow(em) >= 1L)
  expect_true(all(diff(em$kl) <= 1e-12))
})

test_that("EM weights stay on the simplex", {
  p <- plurality_problem()
  set.seed(5)
  st <- ripr_init(
    p$Q,
    p$null,
    exact_engine(),
    control = ripr_control(n_seeds = 30L)
  )
  st <- ripr_step(st, 2L)
  w <- flat_weights(st)
  expect_equal(sum(w), 1, tolerance = 1e-12)
  expect_true(all(w >= 0))
})

# --- Structure ----------------------------------------------------------------

test_that("each new atom lands in the subnull it was found in", {
  p <- plurality_problem()
  set.seed(6)
  st <- ripr_init(
    p$Q,
    p$null,
    exact_engine(),
    control = ripr_control(n_seeds = 30L)
  )
  st <- ripr_step(st, 3L)
  for (i in seq_along(st@atoms)) {
    for (j in seq_len(ncol(st@atoms[[i]]))) {
      expect_true(contains(
        p$null@subnulls[[i]],
        st@atoms[[i]][, j],
        tol = 1e-5
      ))
    }
  }
})

test_that("flat views round-trip through the per-subnull lists", {
  p <- plurality_problem()
  set.seed(7)
  st <- ripr_init(
    p$Q,
    p$null,
    exact_engine(),
    control = ripr_control(n_seeds = 20L)
  )
  st <- ripr_step(st, 2L)
  expect_equal(ncol(flat_atoms(st)), length(flat_weights(st)))
  expect_equal(length(flat_subnull(st)), length(flat_weights(st)))
  expect_equal(
    vapply(unflatten_weights(st, flat_weights(st)), length, integer(1)),
    vapply(st@atoms, ncol, integer(1))
  )
})

# --- Resumption ---------------------------------------------------------------

test_that("stepping twice matches stepping once with twice the count", {
  p <- binomial_problem()
  ctl <- ripr_control(n_seeds = 20L, gap_tol = -Inf)
  set.seed(8)
  a <- ripr_init(p$Q, p$null, exact_engine(), control = ctl) |>
    ripr_step(2L)
  set.seed(8)
  b <- ripr_init(p$Q, p$null, exact_engine(), control = ctl) |>
    ripr_step(1L) |>
    ripr_step(1L)
  expect_equal(flat_weights(a), flat_weights(b), tolerance = 1e-8)
  expect_equal(a@iter, b@iter)
})

# --- The trace ----------------------------------------------------------------

test_that("the trace records init, fw, em and outer phases", {
  p <- plurality_problem()
  set.seed(9)
  st <- ripr_init(
    p$Q,
    p$null,
    exact_engine(),
    control = ripr_control(n_seeds = 20L)
  )
  st <- ripr_step(st, 2L)
  expect_true(all(c("init", "fw", "em", "outer") %in% st@trace$phase))
  expect_true(all(!is.na(st@trace$kl)))
  expect_equal(st@trace$iter[1], 0L)
})

test_that("snapshot control governs how much is stored", {
  p <- binomial_problem()
  set.seed(10)
  run <- function(mode) {
    ctl <- ripr_control(n_seeds = 20L, snapshot = mode, gap_tol = -Inf)
    length(
      ripr_step(
        ripr_init(p$Q, p$null, exact_engine(), control = ctl),
        2L
      )@snapshots
    )
  }
  none <- run("none")
  all <- run("all")
  outer <- run("outer")
  expect_equal(none, 0L)
  expect_gt(all, outer)
  expect_gt(outer, 0L)
})

test_that("snapshots hold the mixture at the recorded moment", {
  p <- binomial_problem()
  set.seed(11)
  ctl <- ripr_control(n_seeds = 20L, snapshot = "all", gap_tol = -Inf)
  st <- ripr_step(ripr_init(p$Q, p$null, exact_engine(), control = ctl), 1L)
  last <- st@snapshots[[length(st@snapshots)]]
  expect_equal(last$atoms, st@atoms)
  expect_equal(last$weights, st@weights)
})

# --- Finishing ----------------------------------------------------------------

test_that("ripr_finish returns a mixing measure and its mixture", {
  p <- plurality_problem()
  set.seed(12)
  st <- ripr_init(
    p$Q,
    p$null,
    exact_engine(),
    control = ripr_control(n_seeds = 30L)
  )
  fit <- ripr_finish(ripr_step(st, 3L))

  expect_true(S7_inherits(fit$W0, finite_mixing))
  expect_true(S7_inherits(fit$P_star, mixture))
  expect_equal(sum(weights(fit$W0)), 1, tolerance = 1e-12)
  expect_equal(
    sum(exp(dist_log_density(fit$P_star, support(p$family)))),
    1,
    tolerance = 1e-10
  )
})

test_that("ripr_finish errors rather than returning empty mixture", {
  p <- binomial_problem()
  set.seed(13)
  st <- ripr_step(
    ripr_init(
      p$Q,
      p$null,
      exact_engine(),
      control = ripr_control(n_seeds = 20L, gap_tol = -Inf)
    ),
    2L
  )
  expect_lte(
    n_atoms(ripr_finish(st, prune = 1e-3)$W0),
    n_atoms(ripr_finish(st, prune = 0)$W0)
  )
  expect_error(ripr_finish(st, prune = 1), "no atom has weight above")
})

# --- The one-dimensional case with a known answer -----------------------------

test_that("the binomial RIPr concentrates at the boundary", {
  # H_0: p <= 1/2 against Q = Bin(10, 0.75). The projection is the point mass at
  # p = 1/2, so the fitted mixture should put essentially all weight there.
  p <- binomial_problem(n = 10, p = 0.75)
  set.seed(14)
  st <- ripr_init(
    p$Q,
    p$null,
    exact_engine(),
    control = ripr_control(n_seeds = 100L, em_max_iter = 20L)
  )
  fit <- ripr_finish(ripr_step(st, 8L), prune = 1e-4)

  heavy <- fit$W0@components[, which.max(weights(fit$W0))]
  expect_equal(heavy, c(0.5, 0.5), tolerance = 1e-3)
  expect_gt(max(weights(fit$W0)), 0.99)
})

# --- Frank-Wolfe variants -----------------------------------------------------

test_that("the fixed schedule uses gamma = 2/(k+2)", {
  p <- plurality_problem()
  set.seed(20)
  ctl <- ripr_control(fw_step = "fixed", n_seeds = 20L, em_max_iter = 2L)
  st <- ripr_step(ripr_init(p$Q, p$null, exact_engine(), control = ctl), 3L)
  steps <- st@trace$step_size[st@trace$phase == "fw"]
  k <- seq_along(steps) - 1L
  expect_equal(steps, 2 / (k + 2), tolerance = 1e-12)
  # The first step is gamma = 1: the open-loop schedule discards the
  # initialisation entirely rather than mixing with it.
  expect_equal(steps[1L], 1)
})

test_that("fw_step = 'none' never grows the support", {
  p <- plurality_problem()
  set.seed(21)
  ctl <- ripr_control(fw_step = "none", n_seeds = 20L, gap_tol = -Inf)
  st0 <- ripr_init(p$Q, p$null, exact_engine(), control = ctl)
  st <- ripr_step(st0, 3L)
  expect_equal(length(flat_weights(st)), length(flat_weights(st0)))
  expect_true(all(st@trace$step_size[st@trace$phase == "fw"] == 0))
})

test_that("Control settings that can never change between iterations warns", {
  expect_warning(
    ripr_control(fw_step = "none", em_max_iter = 0L),
    "nothing can change"
  )
})

test_that("unimplemented variants error with a pointer to what works", {
  p <- plurality_problem()
  for (v in c("pairwise", "away")) {
    ctl <- ripr_control(fw_step = v, n_seeds = 10L)
    st <- ripr_init(p$Q, p$null, exact_engine(), control = ctl)
    expect_error(ripr_step(st, 1L), "not implemented")
  }
})

test_that("unimplemented variants error even when reached directly", {
  p <- binomial_problem()
  ctl <- ripr_control(fw_step = "pairwise", n_seeds = 10L)
  st <- ripr_init(p$Q, p$null, exact_engine(), control = ctl)
  ld <- compile_engine(st@engine)
  expect_error(
    fw_update(
      st,
      c(0.5, 0.5),
      1L,
      gap = 1,
      log_p = log_p_at_nodes(st, ld),
      ld = ld
    ),
    "not implemented"
  )
})

test_that("fully-corrective satisfies the first-order conditions", {
  # At the minimiser over the simplex: G(theta_c) = 1 for every atom carrying
  # weight, and G(theta_c) <= 1 for those the solver is driving out. Atoms with
  # a tiny positive weight are in the second group, not the first.
  p <- plurality_problem()
  set.seed(32)
  ctl <- ripr_control(
    fw_step = "fully-corrective",
    em_max_iter = 0L,
    n_seeds = 20L,
    gap_tol = -Inf
  )
  st <- ripr_step(ripr_init(p$Q, p$null, exact_engine(), control = ctl), 4L)

  ld <- compile_engine(st@engine)
  w <- flat_weights(st)
  log_p <- row_logsumexp(add_by_col(ld(flat_atoms(st)), log(w)))
  g <- exp(col_logsumexp(ld(flat_atoms(st)) - log_p + st@engine@log_w))

  tol <- 1e-3

  expect_true(all(g <= 1 + tol))
  live <- w > tol
  expect_gt(sum(live), 0L)
  expect_equal(g[live], rep(1, sum(live)), tolerance = tol)
  # And the identity, which holds regardless of convergence.
  expect_equal(sum(w * g), 1, tolerance = 1e-10)
})

test_that("fully-corrective records its corrective sweeps", {
  p <- plurality_problem()
  set.seed(33)
  ctl <- ripr_control(
    fw_step = "fully-corrective",
    em_max_iter = 0L,
    n_seeds = 20L,
    gap_tol = -Inf
  )
  st <- ripr_step(ripr_init(p$Q, p$null, exact_engine(), control = ctl), 2L)
  fc <- st@trace[st@trace$phase == "fc", ]
  expect_gt(nrow(fc), 0L)
  expect_true(all(diff(fc$kl[fc$iter == max(fc$iter)]) <= 1e-12))
})


test_that("ripr_control rejects an unknown variant", {
  expect_error(ripr_control(fw_step = "madeup"), "should be one of")
})

test_that("em_max_iter = 0 gives plain Frank-Wolfe with no corrective step", {
  p <- plurality_problem()
  set.seed(31)
  ctl <- ripr_control(em_max_iter = 0L, n_seeds = 20L, gap_tol = -Inf)
  st <- ripr_step(ripr_init(p$Q, p$null, exact_engine(), control = ctl), 3L)
  expect_equal(nrow(st@trace[st@trace$phase == "em", ]), 0L)
  expect_equal(weighted_gain(st), 1, tolerance = 1e-10)
})

# --- EM on atoms --------------------------------------------------------------

test_that("EM moves atoms and keeps them in their subnulls", {
  p <- plurality_problem()
  set.seed(22)
  ctl <- ripr_control(
    fw_step = "none",
    em_max_iter = 3L,
    n_seeds = 10L,
    gap_tol = -Inf
  )
  st0 <- ripr_init(p$Q, p$null, exact_engine(), control = ctl)
  st <- ripr_step(st0, 2L)

  expect_false(isTRUE(all.equal(flat_atoms(st), flat_atoms(st0))))
  for (i in seq_along(st@atoms)) {
    for (j in seq_len(ncol(st@atoms[[i]]))) {
      expect_true(contains(
        p$null@subnulls[[i]],
        st@atoms[[i]][, j],
        tol = 1e-5
      ))
    }
  }
})

test_that("fully-corrective with em_max_iter = 0 leaves the atoms where they were", {
  # The weights-only configuration: the same convex subproblem the EM weight
  # step solves, but run to a convergence test rather than a fixed count.
  p <- plurality_problem()
  set.seed(23)
  ctl <- ripr_control(
    fw_step = "fully-corrective",
    em_max_iter = 0L,
    n_seeds = 10L,
    gap_tol = -Inf
  )
  st0 <- ripr_init(p$Q, p$null, exact_engine(), control = ctl)
  st <- ripr_step(st0, 2L)
  # Every original atom is still exactly where it started; only new ones added.
  for (i in seq_along(st@atoms)) {
    n0 <- ncol(st0@atoms[[i]])
    expect_equal(st@atoms[[i]][, seq_len(n0), drop = FALSE], st0@atoms[[i]])
  }
  expect_equal(nrow(st@trace[st@trace$phase == "em", ]), 0L)
})

test_that("moving atoms does not increase KL beyond search error", {
  # Only approximately monotone: the M-step goes through a heuristic search in
  # the chart, so it can lose the chart round-trip.
  p <- plurality_problem()
  set.seed(24)
  ctl <- ripr_control(
    fw_step = "none",
    em_max_iter = 5L,
    n_seeds = 10L,
    gap_tol = -Inf
  )
  st <- ripr_step(ripr_init(p$Q, p$null, exact_engine(), control = ctl), 3L)
  expect_true(all(diff(st@trace$kl) <= 1e-8))
})

test_that("the atom M-step is deterministic", {
  # No random seeding: the same state must give the same moved atoms whatever
  # the RNG is doing, which is what makes a fit reproducible without a seed.
  p <- plurality_problem()
  set.seed(40)
  ctl <- ripr_control(fw_step = "none", em_max_iter = 1L, gap_tol = -Inf)
  st <- ripr_init(p$Q, p$null, exact_engine(), control = ctl)
  ld <- compile_engine(st@engine)
  set.seed(41)
  a <- em_atom_step(st, ld)
  set.seed(42)
  b <- em_atom_step(st, ld)
  expect_equal(flat_atoms(a), flat_atoms(b))
})

test_that("the identity still holds when atoms move", {
  p <- plurality_problem()
  set.seed(25)
  ctl <- ripr_control(em_max_iter = 3L, n_seeds = 20L)
  st <- ripr_init(p$Q, p$null, exact_engine(), control = ctl)
  for (i in 1:2) {
    st <- ripr_step(st, 1L)
    expect_equal(weighted_gain(st), 1, tolerance = 1e-10)
  }
})

test_that("unflatten_atoms inverts flat_atoms", {
  p <- plurality_problem()
  set.seed(26)
  st <- ripr_step(
    ripr_init(
      p$Q,
      p$null,
      exact_engine(),
      control = ripr_control(n_seeds = 20L)
    ),
    2L
  )
  expect_equal(unflatten_atoms(st, flat_atoms(st)), st@atoms)
})

# --- Convergence flag ---------------------------------------------------------

test_that("a converged iteration adds no atom", {
  # The oracle sweep runs, the gap is recorded, but nothing is appended.
  p <- binomial_problem()
  set.seed(27)
  ctl <- ripr_control(n_seeds = 100L, em_max_iter = 20L, gap_tol = 1e-6)
  st <- ripr_step(ripr_init(p$Q, p$null, exact_engine(), control = ctl), 20L)
  expect_true(st@converged)

  size_before <- length(flat_weights(st))
  st2 <- ripr_step(st, 1L)
  expect_equal(length(flat_weights(st2)), size_before)
})

test_that("converged is FALSE when the budget runs out first", {
  p <- plurality_problem()
  set.seed(28)
  ctl <- ripr_control(n_seeds = 20L, gap_tol = 1e-14)
  st <- ripr_step(ripr_init(p$Q, p$null, exact_engine(), control = ctl), 2L)
  expect_false(st@converged)
  expect_false(ripr_finish(st)$converged)
})

# --- Reference points ---------------------------------------------------------

test_that("initialisation is independent of the engine", {
  # ref comes from W_1, not from the quadrature nodes, so a stochastic engine
  # no longer makes the starting mixture seed-dependent.
  p <- plurality_problem()
  set.seed(29)
  a <- ripr_init(p$Q, p$null, mc_engine(200L))
  set.seed(30)
  b <- ripr_init(p$Q, p$null, mc_engine(200L))
  expect_equal(a@atoms, b@atoms)
  expect_equal(a@atoms, ripr_init(p$Q, p$null, exact_engine())@atoms)
})

test_that("mode_parameter returns a charged parameter", {
  expect_equal(
    mode_parameter(point_mixing(theta_star = c(0.3, 0.7))),
    c(0.3, 0.7)
  )
  f <- finite_mixing(
    components = cbind(c(0.5, 0.5), c(0.1, 0.9)),
    weights = c(0.25, 0.75)
  )
  expect_equal(mode_parameter(f), c(0.1, 0.9))
})

test_that("reference_parameter is a valid interior parameter", {
  expect_equal(reference_parameter(multinomial_family(5, 4)), rep(0.25, 4))
  expect_equal(reference_parameter(gaussian_family(dim = 3)), rep(0, 3))
})
