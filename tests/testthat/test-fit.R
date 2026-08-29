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
  Q <- induced_distribution(fam, dirac(theta = q))
  parts <- lapply(2:k, function(j) {
    basis <- lapply(setdiff(seq_len(k), 1L), function(i) {
      replace(numeric(k), i, 1)
    })
    tie <- replace(numeric(k), c(1L, j), 0.5)
    simplex_region(vertices = do.call(cbind, c(basis, list(tie))))
  })
  ripr_init(
    Q,
    null_model(fam, parts),
    exact_engine(),
    control = ripr_control(n_seeds = 30L, n_restarts = 4L, ...)
  )
}

# H_0: p <= 1/2 against Bin(10, 0.75). The projection is the point mass at 1/2,
# so this is the one problem here with a known answer.
binomial <- function(p = 0.75, ...) {
  fam <- multinomial_family(n_trials = 10, k = 2)
  Q <- induced_distribution(fam, dirac(theta = c(p, 1 - p)))
  ripr_init(
    Q,
    null_model(
      fam,
      list(simplex_region(vertices = cbind(c(0, 1), c(0.5, 0.5))))
    ),
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

in_own_part <- function(state) {
  all(vapply(
    seq_along(state@atoms),
    function(i) {
      atoms <- state@atoms[[i]]
      all(vapply(
        seq_len(ncol(atoms)),
        function(j) {
          contains(parts(state@null@region)[[i]], atoms[, j], tol = 1e-5)
        },
        logical(1)
      ))
    },
    logical(1)
  ))
}

# --- Initialisation -----------------------------------------------------------

test_that("initial atoms lie in their own parts", {
  expect_true(in_own_part(plurality()))
})

test_that("initialisation does not depend on the engine's randomness", {
  # The reference point comes from the alternative, not from the quadrature
  # nodes, so a stochastic engine must not make the starting mixture depend on
  # the seed.
  fam <- multinomial_family(n_trials = 12, k = 4)
  Q <- induced_distribution(fam, dirac(c(0.42, 0.31, 0.16, 0.11)))
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
    expect_equal(weighted_gain(st), 1, tolerance = rounding_tol(1))
    expect_equal(sum(flat_weights(st)), 1, tolerance = rounding_tol(1))
    expect_true(all(flat_weights(st) >= 0))
  }
})

test_that("every atom stays in the part it was found in", {
  st <- plurality() |> fw_step(3L) |> em_step(3L) |> lb_step(1L)
  expect_true(in_own_part(st))
})

# --- Monotonicity -------------------------------------------------------------

test_that("KL never increases under a line search", {
  # Guaranteed because gamma = 0 is in the search interval: the step can always
  # decline to move. Not true of the fixed schedule -- see below.
  st <- plurality() |> fw_step(5L) |> em_step(5L) |> weight_step(5L)
  expect_true(all(diff(st@trace$kl) <= rounding_tol(1)))
})

test_that("the fixed schedule does not discard a seeded initialisation", {
  # Frank-Wolfe opens at gamma = 1, replacing the iterate with the atom the
  # oracle just found. Sensible from an empty support, ruinous from a seeded
  # one: `ripr_init` places a considered atom per part, while the oracle
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
  added <- !is.na(fw$part)
  expect_true(all(added))
  expect_equal(fw$step_size, 2 / (fw$support_size + 1))
})

# --- times and until ----------------------------------------------------------

test_that("until stops early and times remains a ceiling", {
  st <- fw_step(plurality(), 50L, until = function(s) s@iters[["fw"]] >= 3L)
  expect_identical(st@iters[["fw"]], 3L)
})

test_that("gap_below stops on the recorded gap", {
  st <- fw_step(plurality(), 30L, record_gap = TRUE, until = gap_below(0.5))
  expect_true(utils::tail(st@trace$gap, 1L) < 0.5)
  expect_true(st@iters[["fw"]] < 30L)
})

test_that("gap_below refuses a stale gap", {
  # No verb records a gap unasked, so the last recorded one belongs to an
  # earlier iterate. Answering from it would silently describe a mixture that
  # no longer exists.
  st <- em_step(fw_step(plurality(), 2L), 1L)
  expect_error(gap_below(1e-8)(st), "no gap recorded")
})

test_that("record_gap makes a gap available to the predicate", {
  st <- em_step(fw_step(plurality(), 2L), 1L, record_gap = TRUE)
  expect_silent(gap_below(1e-8)(st))
  expect_true(!is.na(utils::tail(st@trace$gap, 1L)))
  # And through the weight verb too, which sweeps after its step.
  st <- weight_step(st, 1L, record_gap = TRUE)
  expect_true(!is.na(utils::tail(st@trace$gap, 1L)))
})


test_that("ripr_init refuses an atoms list that mismatches the parts", {
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  Q <- induced_distribution(fam, dirac(theta = c(0.5, 0.3, 0.2)))
  null <- null_model(
    fam,
    list(
      simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
      simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
    )
  )
  expect_error(
    ripr_init(Q, null, atoms = list(matrix(c(0.2, 0.5, 0.3), ncol = 1L))),
    "one element per part"
  )
  expect_error(
    ripr_init(
      Q,
      null,
      atoms = list(
        matrix(numeric(0), nrow = 3L, ncol = 0L),
        matrix(numeric(0), nrow = 3L, ncol = 0L)
      )
    ),
    "at least one part must carry an atom"
  )
})

# --- What a trace row measures ------------------------------------------------
#
# A row spans a step, so it touches two mixtures. `oracle_value`/`oracle_theta`
# describe the one it started from; `kl`/`gap`/`gap_theta` the one it produced.
# Conflating the two makes `kl - log1p(gap)` -- the guaranteed log-growth rate
# of the resulting e-variable -- a statement about no mixture at all.

test_that("an oracle row's gap measures the mixture the step produced", {
  st <- fw_step(plurality(), 4L, record_gap = TRUE)
  tr <- st@trace[st@trace$phase == "fw", ]

  # The pre-step gap is `oracle_value - 1` and is a different number. Equality
  # would mean `gap` had been copied off the oracle rather than swept after it.
  expect_false(isTRUE(all.equal(tr$gap, tr$oracle_value - 1)))

  # Swept independently at the mixture the last row produced, it agrees.
  ld <- compile_engine(st@engine)
  fresh <- linear_gap(st, log_p_at_nodes(st, ld), ld, flat_atoms(st))
  expect_equal(utils::tail(tr$gap, 1L), fresh$gap, tolerance = rounding_tol(1))
})

test_that("an lb row's gap is swept, not carried over from the row before", {
  # The old pre-step sweep made an `lb` row's gap identical to the previous
  # row's: it paid for a full oracle sweep to recompute a number already there.
  st <- em_step(fw_step(plurality(), 2L), 1L, record_gap = TRUE)
  st <- lb_step(st, 1L, record_gap = TRUE)
  gaps <- stats::na.omit(st@trace$gap)
  expect_false(isTRUE(all.equal(gaps[length(gaps)], gaps[length(gaps) - 1L])))
})

test_that("oracle_theta is the atom the step added", {
  st <- fw_step(plurality(), 4L)
  rows <- st@trace[st@trace$phase == "fw" & !is.na(st@trace$part), ]
  expect_true(nrow(rows) > 0L)

  # `part` names the block the atom entered, so the atom must be in it. The
  # atoms have not moved: only an em or weight sweep would shift them.
  for (i in seq_len(nrow(rows))) {
    block <- st@atoms[[rows$part[i]]]
    gaps <- sqrt(colSums((block - rows$oracle_theta[[i]])^2))
    expect_identical(min(gaps), 0)
  }
})

test_that("an away step records where the oracle looked but adds nothing", {
  # `oracle_theta` says what the oracle proposed; `part` says whether the
  # step took it. An away step proposes a point and then moves the other way.
  st <- fw_step(plurality(), 12L, directions = c("forward", "away"))
  away <- st@trace[!is.na(st@trace$direction) & st@trace$direction == "away", ]
  skip_if(nrow(away) == 0L, "no away step was taken")
  expect_true(all(is.na(away$part)))
  expect_true(all(!is.na(away$oracle_theta)))
})

test_that("theta columns hold the parameter itself, one per row", {
  st <- em_step(fw_step(plurality(k = 4), 1L, record_gap = TRUE), 1L)
  tr <- st@trace
  expect_true(is.list(tr$gap_theta))
  expect_true(is.list(tr$oracle_theta))

  # Whole parameters, not flattened coordinates: a family free to make the
  # parameter something other than a length-4 vector needs no schema change.
  expect_identical(lengths(tr$oracle_theta[tr$phase == "fw"]), 4L)

  # A row that recorded no such point says so the way the rest of the trace
  # does, so `is.na()` reads these columns like any other.
  expect_true(all(is.na(tr$gap_theta[tr$phase == "init"])))
  expect_true(all(is.na(tr$oracle_theta[tr$phase == "em"])))
  expect_true(all(!is.na(tr$oracle_theta[tr$phase == "fw"])))
})

test_that("a trace row survives a parameter that is not a vector", {
  # The list column exists for this: nothing in the trace's type says a
  # parameter is a numeric vector, so a matrix-valued one records unchanged.
  st <- fw_step(plurality(), 1L)
  covariance <- matrix(c(2, 0.5, 0.5, 1), nrow = 2L)
  st <- record(st, phase = "fw", kl = 1, oracle_theta = covariance)
  expect_identical(utils::tail(st@trace$oracle_theta, 1L)[[1L]], covariance)
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
  expect_equal(flat_weights(a), flat_weights(b), tolerance = rounding_tol(1))
  expect_equal(flat_atoms(a), flat_atoms(b), tolerance = rounding_tol(1))
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
  expect_equal(weighted_gain(st), 1, tolerance = rounding_tol(1))
  expect_true(in_own_part(st))
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
  expect_equal(sum(weights(fit$W0)), 1, tolerance = rounding_tol(1))
  expect_equal(
    sum(exp(log_density(
      fit$P_star,
      enumerate_space(fit$P_star@family@sample_space)
    ))),
    1,
    tolerance = rounding_tol(1)
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


# --- A null that is convex but not a simplex ----------------------------------

test_that("a polytope null fits and certifies end to end", {
  # Nothing the caller writes here is a simplex: the null is the polytope
  # `{theta_1 <= 1/2, theta_2 <= 1/2}` in the 2-simplex, stated as its four
  # vertices. `cells()` cuts it into two triangles.
  set.seed(11)
  fam <- multinomial_family(n_trials = 8L, k = 3L)
  square <- polytope_region(
    vertices = cbind(
      c(0.5, 0.5, 0),
      c(0, 0.5, 0.5),
      c(0, 0, 1),
      c(0.5, 0, 0.5)
    )
  )
  null <- null_model(fam, list(square))
  Q <- induced_distribution(fam, dirac(theta = c(0.7, 0.2, 0.1)))

  fit <- ripr_init(
    Q,
    null,
    exact_engine(),
    control = ripr_control(n_seeds = 50L, n_restarts = 5L)
  ) |>
    fw_step(8L, record_gap = TRUE) |>
    em_step(8L) |>
    weight_step(8L) |>
    ripr_finish(record_gap = TRUE)

  # Every atom is in the null, and filed under the only part there is.
  expect_true(all(fit$part == 1L))
  for (j in seq_len(ncol(fit$W0@components))) {
    expect_true(contains(square, fit$W0@components[, j], tol = 1e-6))
  }

  X <- likelihood(Q) / likelihood(fit$P_star)
  cert <- certify(X, null, tol = 1e-9)
  expect_true(all(cert$converged))
  expect_gte(cert$sup_ub, cert$sup_lb)

  # The same identity the simplex nulls satisfy: the certified bound lands at
  # `1 + gap`, with `gap` the duality gap the fit stopped on. Nothing about the
  # decomposition disturbs it.
  expect_equal(cert$sup_ub - 1, fit$gap_final, tolerance = 1e-4)
})
