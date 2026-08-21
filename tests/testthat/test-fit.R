# Properties of R/fit.R.
#
# The public api: build a state, advance it with verbs, convert it to a
# mixture. There is no fixed pipeline, so what is tested here is that each verb
# advances only its own counter, honours `times` and `until`, and leaves the
# state in a condition the next verb can use whatever order they are called in.
#
# The identity `sum_c w_c G(theta_c) = 1` recurs. It is algebraic, so it holds
# after every verb and fails only if the arithmetic and the state have come
# apart.

plurality <- function(k = 4, q = c(0.42, 0.31, 0.16, 0.11), ...) {
  fam <- multinomial_family(n_trials = 12, k = k)
  Q <- mixture(point_mixing(theta_star = q), fam)
  subnulls <- lapply(2:k, function(j) {
    basis <- lapply(setdiff(seq_len(k), 1L), function(i) {
      replace(numeric(k), i, 1)
    })
    tie <- replace(numeric(k), c(1L, j), 0.5)
    simplex_region(vertices = do.call(cbind, c(basis, list(tie))))
  })
  ripr_init(
    Q,
    null_model(fam, subnulls),
    exact_engine(),
    control = ripr_control(n_seeds = 30L, n_restarts = 4L, ...)
  )
}

# H_0: p <= 1/2 against Bin(10, 0.75). The projection is the point mass at 1/2,
# so this is the one problem here with a known answer.
binomial <- function(p = 0.75, ...) {
  fam <- multinomial_family(n_trials = 10, k = 2)
  Q <- mixture(point_mixing(theta_star = c(p, 1 - p)), fam)
  ripr_init(
    Q,
    null_model(fam, list(simplex_region(vertices = cbind(c(0, 1), c(0.5, 0.5))))),
    exact_engine(),
    control = ripr_control(...)
  )
}

weighted_gain <- function(state) {
  ld <- compile_engine(state@engine)
  w <- flat_weights(state)
  ld_all <- ld(flat_atoms(state))
  sum(
    w *
      exp(col_logsumexp(ld_all - mixture_log_p(ld_all, w) + state@engine@log_w))
  )
}

in_own_subnull <- function(state) {
  all(vapply(
    seq_along(state@atoms),
    function(i) {
      atoms <- state@atoms[[i]]
      all(vapply(
        seq_len(ncol(atoms)),
        function(j) contains(state@null@subnulls[[i]], atoms[, j], tol = 1e-5),
        logical(1)
      ))
    },
    logical(1)
  ))
}

# --- Initialisation -----------------------------------------------------------

test_that("ripr_init places one atom per subnull and starts every counter at 0", {
  st <- plurality()
  expect_length(st@atoms, 3L)
  expect_identical(block_sizes(st), rep(1L, 3L))
  expect_equal(sum(flat_weights(st)), 1)
  expect_identical(st@iters, c(fw = 0L, lb = 0L, em = 0L, weight = 0L))
})

test_that("initial atoms lie in their own subnulls", {
  expect_true(in_own_subnull(plurality()))
})

test_that("initialisation does not depend on the engine's randomness", {
  # The reference point comes from the alternative, not from the quadrature
  # nodes, so a stochastic engine must not make the starting mixture depend on
  # the seed.
  fam <- multinomial_family(n_trials = 12, k = 4)
  Q <- mixture(point_mixing(c(0.42, 0.31, 0.16, 0.11)), fam)
  sub <- list(simplex_region(
    vertices = cbind(
      c(0, 1, 0, 0),
      c(0, 0, 1, 0),
      c(0, 0, 0, 1),
      c(0.5, 0.5, 0, 0)
    )
  ))
  H0 <- null_model(fam, sub)
  set.seed(1)
  a <- ripr_init(Q, H0, mc_engine(200L))
  set.seed(2)
  b <- ripr_init(Q, H0, mc_engine(200L))
  expect_equal(a@atoms, b@atoms)
  expect_equal(a@atoms, ripr_init(Q, H0, exact_engine())@atoms)
})

# --- Counters -----------------------------------------------------------------

test_that("each verb advances only its own counter", {
  # Their sum is the total work done, which is what makes `kl` plottable against
  # any one of them.
  st <- plurality() |>
    fw_step(2L) |>
    em_step(3L) |>
    weight_step(4L) |>
    lb_step(1L)
  expect_identical(st@iters, c(fw = 2L, lb = 1L, em = 3L, weight = 4L))
})

test_that("a trace row carries the counters as at that row", {
  st <- fw_step(plurality(), 2L)
  st <- em_step(st, 2L)
  fw_rows <- st@trace[st@trace$phase == "fw", ]
  em_rows <- st@trace[st@trace$phase == "em", ]
  expect_identical(fw_rows$fw, c(1L, 2L))
  expect_identical(em_rows$fw, c(2L, 2L))
  expect_identical(em_rows$em, c(1L, 2L))
})

# --- The algebraic identity ---------------------------------------------------

test_that("the identity survives every verb", {
  # Fails if the weights and atoms come apart, or if a verb writes back weights
  # that do not correspond to the atoms it left behind.
  st <- plurality()
  for (advance in list(
    function(s) fw_step(s, 2L),
    function(s) fw_step(s, 2L, directions = c("forward", "away")),
    function(s) em_step(s, 2L),
    function(s) weight_step(s, 5L),
    function(s) lb_step(s, 1L)
  )) {
    st <- advance(st)
    expect_equal(weighted_gain(st), 1, tolerance = 1e-10)
    expect_equal(sum(flat_weights(st)), 1, tolerance = 1e-12)
    expect_true(all(flat_weights(st) >= 0))
  }
})

test_that("every atom stays in the subnull it was found in", {
  st <- plurality() |> fw_step(3L) |> em_step(3L) |> lb_step(1L)
  expect_true(in_own_subnull(st))
})

# --- Monotonicity -------------------------------------------------------------

test_that("KL never increases under a line search", {
  # Guaranteed because gamma = 0 is in the search interval: the step can always
  # decline to move. Not true of the fixed schedule -- see below.
  st <- plurality() |> fw_step(5L) |> em_step(5L) |> weight_step(5L)
  expect_true(all(diff(st@trace$kl) <= 1e-9))
})

test_that("the fixed schedule does not discard a seeded initialisation", {
  # Frank-Wolfe opens at gamma = 1, replacing the iterate with the atom the
  # oracle just found. Sensible from an empty support, ruinous from a seeded
  # one: `ripr_init` places a considered atom per subnull, while the oracle
  # returns the worst-case theta, which puts near-zero mass where Q has some.
  # Indexing the schedule by the component count rather than the step count
  # avoids it -- with three initial atoms the first step is 2/5, not 1.
  st <- plurality()
  n0 <- length(flat_weights(st))
  stepped <- fw_step(st, 6L, size = "fixed")
  fw <- stepped@trace[stepped@trace$phase == "fw", ]
  expect_equal(fw$step_size[1L], 2 / (n0 + 2))
  # Left at gamma = 1 this reached KL of 36 against a starting 0.37.
  expect_true(max(fw$kl) < 10 * st@trace$kl[1L])
})

test_that("the fixed schedule follows Jaggi's sequence in the component count", {
  # gamma = 2/(k+2) where k is the number of components before the step, so a
  # step that adds one moves the sequence on by exactly one place.
  st <- fw_step(plurality(), 5L, size = "fixed")
  fw <- st@trace[st@trace$phase == "fw", ]
  added <- !is.na(fw$subnull)
  expect_true(all(added))
  expect_equal(fw$step_size, 2 / (fw$support_size + 1))
})

# --- times and until ----------------------------------------------------------

test_that("times is honoured exactly when no predicate is given", {
  expect_identical(fw_step(plurality(), 7L)@iters[["fw"]], 7L)
  expect_identical(em_step(plurality(), 5L)@iters[["em"]], 5L)
})

test_that("until stops early and times remains a ceiling", {
  st <- fw_step(plurality(), 50L, until = function(s) s@iters[["fw"]] >= 3L)
  expect_identical(st@iters[["fw"]], 3L)
})

test_that("gap_below stops on the recorded gap", {
  st <- fw_step(plurality(), 30L, until = gap_below(0.5))
  expect_true(utils::tail(st@trace$gap, 1L) < 0.5)
  expect_true(st@iters[["fw"]] < 30L)
})

test_that("gap_below refuses a stale gap", {
  # `em_step` records no gap by default, so after one the last recorded gap
  # belongs to an earlier iterate. Answering from it would silently describe a
  # mixture that no longer exists.
  st <- em_step(fw_step(plurality(), 2L), 1L)
  expect_error(gap_below(1e-8)(st), "no gap recorded")
})

test_that("record_gap makes a gap available to the predicate", {
  st <- em_step(fw_step(plurality(), 2L), 1L, record_gap = TRUE)
  expect_silent(gap_below(1e-8)(st))
  expect_true(!is.na(utils::tail(st@trace$gap, 1L)))
})

test_that("support_gap_below needs no oracle sweep", {
  # It reads `max_c G(theta_c) - 1` off the current atoms, so it is available
  # whatever the last verb recorded.
  st <- em_step(fw_step(plurality(), 2L), 1L)
  expect_silent(support_gap_below(1e-8)(st))
  expect_true(is.logical(support_gap_below(1e-8)(st)))
})

test_that("kl_flat needs two rows before it can fire", {
  st <- plurality()
  expect_false(kl_flat(1e-6)(st))
  expect_true(is.logical(kl_flat(1e-300)(fw_step(st, 2L))))
})

# --- Composition --------------------------------------------------------------

test_that("splitting a call in two matches taking it in one", {
  # Nothing is carried between iterations that a fresh call would not rebuild,
  # so a fit can be resumed.
  set.seed(8)
  a <- fw_step(plurality(), 3L)
  set.seed(8)
  b <- fw_step(fw_step(plurality(), 1L), 2L)
  expect_equal(flat_weights(a), flat_weights(b), tolerance = 1e-8)
  expect_equal(flat_atoms(a), flat_atoms(b), tolerance = 1e-8)
  expect_identical(a@iters, b@iters)
})

test_that("verbs compose in any order", {
  # There is no fixed pipeline, so an EM sweep before any Frank-Wolfe step, or a
  # weight solve between two of them, must simply work.
  st <- plurality() |>
    em_step(2L) |>
    weight_step(3L) |>
    fw_step(2L) |>
    weight_step(3L) |>
    em_step(2L)
  expect_equal(weighted_gain(st), 1, tolerance = 1e-10)
  expect_true(in_own_subnull(st))
})

test_that("only the stepping verbs can grow the support", {
  # EM moves atoms and the weight solve reweights them; neither adds or drops
  # one.
  st <- fw_step(plurality(), 3L)
  before <- block_sizes(st)
  expect_identical(block_sizes(em_step(st, 3L)), before)
  expect_identical(block_sizes(weight_step(st, 10L)), before)
  expect_true(sum(block_sizes(fw_step(st, 1L))) >= sum(before))
})

# --- Directions ---------------------------------------------------------------

test_that("directions names a choice, and the trace says which was taken", {
  st <- fw_step(plurality(), 6L, directions = c("forward", "away"))
  taken <- unique(stats::na.omit(st@trace$direction))
  expect_true(all(taken %in% c("forward", "away")))
})

test_that("an away step adds no atom and records no subnull", {
  st <- fw_step(plurality(), 8L, directions = c("forward", "away"))
  away <- st@trace[!is.na(st@trace$direction) & st@trace$direction == "away", ]
  if (nrow(away) > 0L) {
    expect_true(all(is.na(away$subnull)))
  }
  expect_true(all(
    !is.na(st@trace$subnull[
      !is.na(st@trace$direction) & st@trace$direction == "forward"
    ])
  ))
})

test_that("a misspelt direction is caught with a suggestion", {
  expect_error(fw_step(plurality(), directions = "awya"), "Did you mean")
  expect_error(fw_step(plurality(), size = "linesearch"), "must be one of")
})

# --- Snapshots ----------------------------------------------------------------

test_that("snapshot counts calls under step and iterations under all", {
  expect_length(fw_step(plurality(snapshot = "none"), 4L)@snapshots, 0L)
  expect_length(fw_step(plurality(snapshot = "step"), 4L)@snapshots, 1L)
  expect_length(fw_step(plurality(snapshot = "all"), 4L)@snapshots, 4L)
  # Which means composition is how the granularity is chosen.
  st <- plurality(snapshot = "step")
  expect_length(fw_step(fw_step(st, 1L), 1L)@snapshots, 2L)
})

# --- Finishing ----------------------------------------------------------------

test_that("ripr_finish is a pure conversion by default", {
  # Refining is opt-in: what comes back is what was fitted.
  st <- fw_step(plurality(), 5L)
  fit <- ripr_finish(st)
  expect_equal(fit$kl, kl_divergence(st))
  expect_equal(unname(weights(fit$W0)), unname(flat_weights(st)))
  expect_identical(fit$rounds, 0L)
})

test_that("kl describes the returned mixture, not the state", {
  # Refining and pruning both change the mixture, so a `kl` taken from the state
  # would describe something the caller was not given.
  st <- fw_step(plurality(), 5L)
  fit <- ripr_finish(st, reoptimise = TRUE, identify = TRUE)
  ld <- compile_engine(st@engine)
  expect_equal(
    fit$kl,
    expect_q(
      st@engine,
      st@engine@log_q - mixture_log_p(ld(fit$W0@components), weights(fit$W0))
    )
  )
})

test_that("refining lowers KL and shrinks the support", {
  st <- fw_step(plurality(), 6L)
  plain <- ripr_finish(st)
  refined <- ripr_finish(st, reoptimise = TRUE, identify = TRUE)
  expect_true(refined$kl <= plain$kl)
  expect_true(n_atoms(refined$W0) <= n_atoms(plain$W0))
  expect_true(refined$rounds >= 1L)
})

test_that("prune keeps, drops, or errors", {
  st <- fw_step(plurality(), 5L)
  identified <- ripr_finish(st, identify = TRUE)
  expect_true(
    n_atoms(ripr_finish(st, identify = TRUE, prune = -1)$W0) >=
      n_atoms(identified$W0)
  )
  expect_true(
    n_atoms(ripr_finish(st, prune = 0.2)$W0) <= n_atoms(identified$W0)
  )
  expect_error(ripr_finish(st, prune = 1), "must be below 1")
  expect_error(ripr_finish(st, prune = 0.99), "no atom has weight above")
})

test_that("the returned mixture is a distribution", {
  fit <- ripr_finish(fw_step(plurality(), 4L))
  expect_true(S7_inherits(fit$W0, finite_mixing))
  expect_true(S7_inherits(fit$P_star, mixture))
  expect_equal(sum(weights(fit$W0)), 1, tolerance = 1e-12)
  expect_equal(
    sum(exp(dist_log_density(fit$P_star, enumerate_space(fit$P_star@family@sample_space)))),
    1,
    tolerance = 1e-10
  )
})

test_that("gap_final measures the returned mixture, gap_fit the state", {
  # `gap_fit` cannot show what finishing did -- it was recorded on the way in.
  st <- fw_step(plurality(), 5L)
  fit <- ripr_finish(st, reoptimise = TRUE, identify = TRUE, record_gap = TRUE)
  expect_equal(fit$gap_fit, utils::tail(st@trace$gap, 1L))
  expect_true(!is.na(fit$gap_final))
  expect_true(is.na(ripr_finish(st)$gap_final))
})

# --- The one problem with a known answer --------------------------------------

test_that("the binomial projection concentrates at the boundary", {
  # H_0: p <= 1/2 against Bin(10, 0.75). The projection is the point mass at
  # p = 1/2, so nearly all the weight belongs on the tie.
  fit <- binomial(p = 0.75, n_seeds = 100L) |>
    fw_step(8L) |>
    em_step(20L) |>
    ripr_finish(reoptimise = TRUE, identify = TRUE)
  heavy <- fit$W0@components[, which.max(weights(fit$W0))]
  expect_equal(heavy, c(0.5, 0.5), tolerance = 1e-3)
  expect_true(max(weights(fit$W0)) > 0.99)
})
