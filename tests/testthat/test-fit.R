# Properties of R/fit.R.
#
# The public api: build a state, advance it with verbs, convert it to a
# mixture. There is no fixed pipeline, so what is tested here is that each verb
# records one row per step, honours `times` and `until`, and leaves the
# state in a condition the next verb can use whatever order they are called in.
#
# The identity `sum_c w_c G(theta_c) = 1` recurs. It is algebraic, so it holds
# after every verb and fails only if the arithmetic and the state have come
# apart.

plurality <- function(
  k = 4,
  q = c(0.42, 0.31, 0.16, 0.11),
  record_gap = FALSE,
  atoms = NULL,
  ...
) {
  fam <- multinomial_family(n_trials = 12, k = k)
  Q <- mixture(fam, dirac(theta = q))
  parts <- lapply(2:k, function(j) {
    basis <- lapply(setdiff(seq_len(k), 1L), function(i) {
      replace(numeric(k), i, 1)
    })
    tie <- replace(numeric(k), c(1L, j), 0.5)
    simplex_region(vertices = do.call(rbind, c(basis, list(tie))))
  })
  ripr_init(
    Q,
    null_model(fam, parts),
    exact_engine(),
    atoms = atoms,
    record_gap = record_gap,
    control = ripr_control(n_seeds = 30L, n_restarts = 4L, ...)
  )
}

# H_0: p <= 1/2 against Bin(10, 0.75). The projection is the point mass at 1/2,
# so this is the one problem here with a known answer.
binomial <- function(p = 0.75, ...) {
  fam <- multinomial_family(n_trials = 10, k = 2)
  Q <- mixture(fam, dirac(theta = c(p, 1 - p)))
  ripr_init(
    Q,
    null_model(
      fam,
      list(simplex_region(vertices = rbind(c(0, 1), c(0.5, 0.5))))
    ),
    exact_engine(),
    control = ripr_control(...)
  )
}

# The identity, with `G` computed by brute force (helper-fit.R) rather than by
# the engine, so it checks the state against the densities it claims.
weighted_gain <- function(state) {
  sum(weights(state@mixing) * exact_g(state, atoms(state@mixing)))
}

n_fw <- function(state) sum(state@trace$phase == "fw")

# --- Initialisation -----------------------------------------------------------

test_that("initial atoms lie in their own parts", {
  expect_true(in_own_part(plurality()))
})

test_that("initialisation does not depend on the engine's randomness", {
  # The reference point comes from the alternative, not from the quadrature
  # nodes, so a stochastic engine must not make the starting mixture depend on
  # the seed.
  fam <- multinomial_family(n_trials = 12, k = 4)
  Q <- mixture(fam, dirac(c(0.42, 0.31, 0.16, 0.11)))
  sub <- list(simplex_region(
    vertices = rbind(
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
  expect_equal(a@mixing@atoms, b@mixing@atoms)
  expect_equal(a@mixing@atoms, ripr_init(Q, H0, exact_engine())@mixing@atoms)
})

# --- Step index ---------------------------------------------------------------

test_that("each verb writes one row per step, tagged with its phase", {
  # The row count is the total work done, and the phases say which verb did
  # each part of it, which is what makes `kl` plottable against either.
  st <- plurality() |>
    fw_step(2L) |>
    em_step(3L) |>
    weight_step(4L) |>
    lb_step(1L)
  expect_identical(
    verb_counts(st),
    c(fw = 2L, lb = 1L, em = 3L, weight = 4L)
  )
  expect_identical(st@trace$phase[1L], "init")
})

test_that("a trace row carries its step index", {
  st <- fw_step(plurality(), 2L)
  st <- em_step(st, 2L)
  expect_identical(st@trace$step, 0:4)
  expect_identical(st@trace$phase, c("init", "fw", "fw", "em", "em"))
  # A per-verb count as at each row is a cumulative count over the phases.
  expect_identical(cumsum(st@trace$phase == "fw"), c(0L, 1L, 2L, 2L, 2L))
})

test_that("every row records the wall-clock time it took", {
  # The clock is the fair axis between step rules whose iterations differ in
  # cost by an order of magnitude, and it can only be measured from inside a
  # `times` loop. Per row rather than cumulative, so a subset still reads.
  st <- plurality()
  outer <- system.time(
    st <- st |> fw_step(2L) |> em_step(2L) |> weight_step(1L) |> lb_step(1L)
  )[["elapsed"]]
  tr <- st@trace
  expect_type(tr$elapsed, "double")
  expect_false(anyNA(tr$elapsed))
  expect_true(all(tr$elapsed >= 0))
  # Each row's clock ran inside the call that produced it.
  expect_lte(sum(tr$elapsed[tr$phase != "init"]), outer + 1e-6)
})

test_that("elapsed excludes the until predicate", {
  slow <- function(state, candidate) {
    Sys.sleep(0.1)
    FALSE
  }
  st <- em_step(plurality(), 2L, until = slow)
  expect_true(all(st@trace$elapsed[st@trace$phase == "em"] < 0.1))
})

test_that("an fw step that uses a cached oracle still prices its search", {
  # The search ran inside the previous verb's `record_gap` sweep; using its
  # cached result must not make the fw row look cheaper than its rule.
  st <- em_step(plurality(), 1L, record_gap = TRUE)
  swept <- utils::tail(st@trace$gap_after_elapsed, 1L)
  expect_gt(swept, 0)
  stepped <- fw_step(st, 1L)
  expect_gte(utils::tail(stepped@trace$elapsed, 1L), swept)
})

test_that("`correct` makes one fw call one fully-corrective iteration", {
  # After the step every weight is re-solved, so the weights satisfy the
  # optimality conditions over the atoms: G <= 1 everywhere, = 1 where weighted.
  st <- fw_step(plurality(), 3L)
  one <- fw_step(st, 1L, correct = TRUE)
  g <- exact_g(one, atoms(one@mixing))
  expect_lt(max(g), 1 + 1e-6)
  expect_equal(g[weights(one@mixing) > 0], rep(1, n_atoms(one@mixing)), tolerance = 1e-6)

  added <- one@trace[seq(nrow(st@trace) + 1L, nrow(one@trace)), ]
  expect_identical(nrow(added), 1L)
  expect_identical(added$phase, "fw")
  expect_identical(verb_counts(one)[["weight"]], 0L)
  # The row's KL is the corrected mixture's, not the stepped one's.
  expect_equal(added$kl, exact_kl(one), tolerance = 1e-12)
  expect_lt(added$kl, utils::tail(fw_step(st, 1L)@trace$kl, 1L))
})

# --- The algebraic identity ---------------------------------------------------

test_that("the identity survives every verb", {
  # Fails if the weights and atoms come apart, or if a verb writes back weights
  # that do not correspond to the atoms it left behind.
  st <- plurality()
  for (advance in list(
    function(s) fw_step(s, 2L),
    function(s) fw_step(s, 2L, variant = "away-step"),
    function(s) em_step(s, 2L),
    function(s) weight_step(s, 5L),
    function(s) lb_step(s, 1L)
  )) {
    st <- advance(st)
    expect_equal(weighted_gain(st), 1, tolerance = rounding_tol(1))
    # And the KL the verb recorded is that of the mixture it left.
    expect_equal(
      utils::tail(st@trace$kl, 1L),
      exact_kl(st),
      tolerance = rounding_tol(1)
    )
    expect_equal(sum(weights(st@mixing)), 1, tolerance = rounding_tol(1))
    expect_true(all(weights(st@mixing) >= 0))
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
  n0 <- length(weights(st@mixing))
  stepped <- fw_step(st, 6L, size = "fixed")
  fw <- stepped@trace[stepped@trace$phase == "fw", ]
  expect_equal(fw$step_size[1L], 2 / (n0 + 2))
  # Left at gamma = 1 this reached KL of 36 against a starting 0.37.
  expect_true(max(fw$kl) < 10 * st@trace$kl[1L])
})

test_that("the fixed schedule follows Jaggi's sequence in the step count", {
  # gamma = 2/(k+2), k advancing by one per oracle step. While every step adds
  # an atom that is also the component count, which is what this checks.
  st <- fw_step(plurality(), 5L, size = "fixed")
  fw <- st@trace[st@trace$phase == "fw", ]
  added <- !is.na(fw$part)
  expect_true(all(added))
  expect_equal(fw$step_size, 2 / (fw$support_size + 1))
})

# --- times and until ----------------------------------------------------------

test_that("until stops early and times remains a ceiling", {
  st <- fw_step(
    plurality(),
    50L,
    until = function(state, candidate) n_fw(state) >= 3L
  )
  expect_identical(n_fw(st), 3L)
})

test_that("until must take the state and the candidate", {
  # Refusing a one-argument predicate up front beats an "unused argument"
  # error after the first step's search has run.
  expect_error(
    fw_step(plurality(), 2L, until = function(s) TRUE),
    "function of `\\(state, candidate\\)`"
  )
  expect_error(em_step(plurality(), 2L, until = TRUE), "must be a function")
  # `...` can absorb the candidate.
  expect_no_error(em_step(plurality(), 1L, until = function(...) FALSE))
})

test_that("gap_below stops on the recorded gap", {
  st <- fw_step(plurality(), 30L, until = gap_below(0.5))
  expect_true(utils::tail(st@trace$gap_after, 1L) < 0.5)
  expect_true(st@oracle$value - 1 < 0.5)
  expect_true(n_fw(st) < 30L)
})

test_that("a verb whose `until` already holds takes no step", {
  # Tested after the step instead, a call would step whatever it was handed:
  # the verb would have no fixed point, so composing it with itself would keep
  # stepping and a loop around it would never end.
  p <- function(state, candidate) n_fw(state) >= 3L
  a <- fw_step(plurality(), 10L, until = p)
  expect_identical(n_fw(a), 3L)
  b <- fw_step(a, 10L, until = p)
  expect_identical(n_fw(b), 3L)
  expect_identical(nrow(b@trace), nrow(a@trace))
  # `b`'s oracle filled the last row and stays cached, so a third call uses
  # it: the trace is unchanged and no search draws from the RNG stream.
  expect_false(is.null(b@oracle))
  seed <- .Random.seed
  c <- fw_step(b, 10L, until = p)
  expect_identical(c@trace, b@trace)
  expect_identical(.Random.seed, seed)
})

test_that("a predicate sees the state the step would make", {
  # The pending row appended, so any property of the transition can be read
  # off the pair. An EM sweep is deterministic, so the candidate matches the
  # state an unconditional call makes.
  st <- plurality()
  one <- em_step(st, 1L)
  seen <- NULL
  none <- em_step(st, 5L, until = function(state, candidate) {
    seen <<- candidate
    TRUE
  })
  expect_identical(none@trace, st@trace)
  expect_identical(nrow(seen@trace), nrow(st@trace) + 1L)
  expect_identical(utils::tail(seen@trace$phase, 1L), "em")
  expect_identical(seen@trace$kl, one@trace$kl)
  expect_identical(seen@mixing, one@mixing)
})

test_that("a stalled sweep does not block a verb that can still progress", {
  # EM cannot grow the support: with atoms in one part only, its trace goes
  # flat while a Frank--Wolfe step still adds the atom another part needs.
  # Reseeding before the measured and the asserted runs makes their oracle
  # searches identical, so `first` is the pending KL the predicate sees.
  empty <- matrix(numeric(0), nrow = 0L, ncol = 4L)
  st <- plurality(atoms = list(
    rbind(c(0.125, 0.375, 0.25, 0.25), c(0.25, 0.75, 0, 0)),
    empty,
    empty
  ))
  set.seed(42)
  before <- em_step(st, 59L)
  st <- em_step(before, 1L)
  stalled <- abs(diff(utils::tail(st@trace$kl, 2L)))
  set.seed(7)
  first <- abs(diff(utils::tail(fw_step(st, 1L)@trace$kl, 2L)))
  expect_gt(first, stalled)

  tol <- (stalled + first) / 2
  expect_true(kl_flat(tol)(before, st))
  set.seed(7)
  expect_gt(n_fw(fw_step(st, 5L, until = kl_flat(tol))), n_fw(st))
})

test_that("gap_below refuses a stale gap", {
  # The recorded gaps belong to earlier iterates; answering from one would
  # silently describe a mixture that no longer exists.
  st <- em_step(fw_step(plurality(), 2L), 1L)
  expect_null(st@oracle)
  expect_error(gap_below(1e-8)(st), "no gap recorded")
})

test_that("record_gap makes a gap available to the predicate", {
  st <- em_step(fw_step(plurality(), 2L), 1L, record_gap = TRUE)
  expect_silent(gap_below(1e-8)(st))
  expect_true(!is.na(utils::tail(st@trace$gap_after, 1L)))
  # The cache and the row it filled are one number.
  expect_identical(st@oracle$value - 1, utils::tail(st@trace$gap_after, 1L))
  # And through the weight verb too, which sweeps after its step.
  st <- weight_step(st, 1L, record_gap = TRUE)
  expect_true(!is.na(utils::tail(st@trace$gap_after, 1L)))
})

test_that("ripr_init can record the starting mixture's gap", {
  # So `record_gap = TRUE` throughout leaves no row without one, and the
  # predicate can be asked before any step is taken.
  st <- plurality(record_gap = TRUE)
  expect_false(is.na(st@trace$gap_after))
  expect_false(anyNA(st@trace$gap_after_theta[[1L]]))
  expect_silent(gap_below(1e-8)(st))
  # Unasked, the init row is like any other: no sweep, no gap, no cache.
  expect_true(is.na(plurality()@trace$gap_after))
  expect_null(plurality()@oracle)
})

test_that("every step clears the cached oracle, and record_gap resets it", {
  # The cache belongs to the current mixture, so a verb that moves the
  # mixture without measuring it must leave nothing behind for `fw_step()` or
  # `gap_below()` to mistake for current.
  st <- fw_step(plurality(), 1L, record_gap = TRUE)
  expect_false(is.null(st@oracle))
  for (advance in list(
    function(s) fw_step(s, 1L),
    function(s) lb_step(s, 1L),
    function(s) em_step(s, 1L),
    function(s) weight_step(s, 1L)
  )) {
    expect_null(advance(st)@oracle)
  }
  swept <- weight_step(st, 1L, record_gap = TRUE)
  expect_false(is.null(swept@oracle))
  expect_identical(
    swept@oracle$theta,
    utils::tail(swept@trace$gap_after_theta, 1L)[[1L]]
  )
})

test_that("fw_step uses a cached oracle instead of searching again", {
  # One search per step: the sweep that `record_gap` ran is the one the step
  # consumes, so no search draws from the RNG stream and the atom goes where
  # the sweep said.
  st <- em_step(plurality(), 1L, record_gap = TRUE)
  cached <- st@oracle
  seed <- .Random.seed
  stepped <- fw_step(st, 1L)
  expect_identical(.Random.seed, seed)
  expect_identical(utils::tail(stepped@trace$oracle_theta, 1L)[[1L]], cached$theta)
})


test_that("ripr_init refuses an atoms list that mismatches the parts", {
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  Q <- mixture(fam, dirac(theta = c(0.5, 0.3, 0.2)))
  null <- null_model(
    fam,
    list(
      simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
      simplex_region(vertices = rbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
    )
  )
  expect_error(
    ripr_init(Q, null, atoms = list(matrix(c(0.2, 0.5, 0.3), nrow = 1L))),
    "one element per part"
  )
  expect_error(
    ripr_init(
      Q,
      null,
      atoms = list(
        matrix(numeric(0), nrow = 0L, ncol = 3L),
        matrix(numeric(0), nrow = 0L, ncol = 3L)
      )
    ),
    "at least one part must carry an atom"
  )
})

# --- What a trace row measures ------------------------------------------------
#
# A row spans a step, so it touches two mixtures. `oracle_value`/`oracle_theta`
# describe the one it started from; `kl`/`gap_after`/`gap_after_theta` the one it produced.
# Conflating the two makes `kl - log1p(gap)` -- the guaranteed log-growth rate
# of the resulting e-variable -- a statement about no mixture at all.

test_that("a row's gap is the next step's oracle, not a second search", {
  # The oracle of an fw step maximised `G` over exactly the mixture the row
  # before it produced, so its maximiser is that row's `gap_after_theta` and
  # its value that row's `gap_after + 1`. Storing the value again on the fw
  # row would say nothing new, so `oracle_value` stays NA there.
  st <- fw_step(plurality(), 4L)
  tr <- st@trace
  expect_false(anyNA(tr$gap_after[-nrow(tr)]))
  expect_identical(tr$gap_after_theta[-nrow(tr)], tr$oracle_theta[-1L])
  expect_true(all(is.na(tr$oracle_value[tr$phase == "fw"])))
  # Nothing followed the last row, so nothing has measured it.
  expect_true(is.na(utils::tail(tr$gap_after, 1L)))
})

test_that("record_gap on fw_step completes the final row", {
  st <- fw_step(plurality(), 3L, record_gap = TRUE)
  fw <- st@trace[st@trace$phase == "fw", ]
  expect_false(anyNA(fw$gap_after))
  # Each sweep is the next step's oracle, used from the cache rather than
  # searched twice, and each later row prices the search it consumed.
  expect_identical(fw$gap_after_theta[-nrow(fw)], fw$oracle_theta[-1L])
  expect_true(all(fw$elapsed[-1L] >= fw$gap_after_elapsed[-nrow(fw)]))
})

test_that("a filled gap is the gap of the mixture the row produced", {
  # Measured independently of the bookkeeping that wrote it. `gap_below(Inf)`
  # holds everywhere, so the second call fills the last row and steps nowhere,
  # leaving the mixture the row's gap describes as the current one.
  # `ripr_finish(record_gap = TRUE)` searches afresh over the mixture it
  # returns, which without refinement is the state's own, regrouped by part.
  set.seed(1)
  st <- fw_step(plurality(), 3L)
  st <- fw_step(st, 1L, until = gap_below(Inf))
  fresh <- ripr_finish(st, record_gap = TRUE)@gap_final
  # Loose enough for two independent multi-start searches; a wrong-mixture
  # fill would miss by orders more.
  expect_equal(utils::tail(st@trace$gap_after, 1L), fresh, tolerance = 1e-8)
  # And both are `sup G - 1`, attained where the row says.
  expect_equal(
    utils::tail(st@trace$gap_after, 1L),
    exact_g(st, utils::tail(st@trace$gap_after_theta, 1L)[[1L]]) - 1
  )
})

test_that("an lb row's gap is swept, not carried over from the row before", {
  # The old pre-step sweep made an `lb` row's gap identical to the previous
  # row's: it paid for a full oracle sweep to recompute a number already there.
  st <- em_step(fw_step(plurality(), 2L), 1L, record_gap = TRUE)
  st <- lb_step(st, 1L, record_gap = TRUE)
  gaps <- stats::na.omit(st@trace$gap_after)
  expect_false(isTRUE(all.equal(gaps[length(gaps)], gaps[length(gaps) - 1L])))
})

test_that("oracle_theta is the atom the step added", {
  st <- fw_step(plurality(), 4L)
  rows <- st@trace[st@trace$phase == "fw" & !is.na(st@trace$part), ]
  expect_true(nrow(rows) > 0L)

  # `part` names the part the atom entered, so the atom must be filed under it.
  # The atoms have not moved: only an em or weight sweep would shift them.
  for (i in seq_len(nrow(rows))) {
    block <- st@mixing@atoms[st@part == rows$part[i], , drop = FALSE]
    gaps <- sqrt(rowSums(sweep(block, 2L, rows$oracle_theta[[i]])^2))
    expect_identical(min(gaps), 0)
  }
})

test_that("an away step records where the oracle looked but adds nothing", {
  # `oracle_theta` says what the oracle proposed; `part` says whether the
  # step took it. An away step proposes a point and then moves the other way.
  st <- fw_step(plurality(), 12L, variant = "away-step")
  away <- st@trace[!is.na(st@trace$direction) & st@trace$direction == "away", ]
  skip_if(nrow(away) == 0L, "no away step was taken")
  expect_true(all(is.na(away$part)))
  expect_true(all(!is.na(away$oracle_theta)))
})

test_that("theta columns hold the parameter itself, one per row", {
  st <- em_step(fw_step(plurality(k = 4), 1L), 1L)
  tr <- st@trace
  expect_true(is.list(tr$gap_after_theta))
  expect_true(is.list(tr$oracle_theta))

  # Whole parameters, not flattened coordinates: a family free to make the
  # parameter something other than a length-4 vector needs no schema change.
  expect_identical(lengths(tr$oracle_theta[tr$phase == "fw"]), 4L)

  # A row that recorded no such point says so the way the rest of the trace
  # does, so `is.na()` reads these columns like any other. Nothing stepped
  # after the last row, so nothing has measured its gap.
  expect_true(is.na(tr$gap_after_theta[[nrow(tr)]]))
  expect_true(all(is.na(tr$oracle_theta[tr$phase == "em"])))
  expect_true(all(!is.na(tr$oracle_theta[tr$phase == "fw"])))
})

test_that("support_gap_below needs no oracle sweep", {
  # It reads `max_c G(theta_c) - 1` off the current atoms, so it is available
  # whatever the last verb recorded.
  st <- em_step(fw_step(plurality(), 2L), 1L)
  expect_silent(support_gap_below(1e-8)(st))
  expect_true(is.logical(support_gap_below(1e-8)(st)))
})

test_that("kl_flat compares the pending step's KL with the current one", {
  st <- plurality()
  stepped <- fw_step(st, 1L)
  change <- abs(diff(stepped@trace$kl))
  expect_true(kl_flat(2 * change)(st, stepped))
  expect_false(kl_flat(change / 2)(st, stepped))
})

# --- Composition --------------------------------------------------------------

test_that("splitting a call in two matches taking it in one", {
  # Nothing is carried between iterations that a fresh call would not rebuild,
  # so a fit can be resumed.
  set.seed(8)
  a <- fw_step(plurality(), 3L)
  set.seed(8)
  b <- fw_step(fw_step(plurality(), 1L), 2L)
  expect_equal(weights(a@mixing), weights(b@mixing), tolerance = rounding_tol(1))
  expect_equal(atoms(a@mixing), atoms(b@mixing), tolerance = rounding_tol(1))
  expect_identical(a@trace$step, b@trace$step)
  expect_identical(a@trace$phase, b@trace$phase)
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
  expect_identical(em_step(st, 3L)@part, st@part)
  expect_identical(weight_step(st, 10L)@part, st@part)
  expect_true(length(fw_step(st, 1L)@part) >= length(st@part))
})

# --- Directions ---------------------------------------------------------------

test_that("a misspelt direction is caught with a suggestion", {
  expect_error(fw_step(plurality(), variant = "awaystep"), "Did you mean")
  expect_error(fw_step(plurality(), size = "linesearch"), "must be one of")
})

# --- Snapshots ----------------------------------------------------------------

test_that("snapshot counts calls under step and iterations under all", {
  # `ripr_init()` is one call and one iteration, so it contributes one
  # snapshot under either, and none under "none".
  expect_length(fw_step(plurality(snapshot = "none"), 4L)@snapshots, 0L)
  expect_length(fw_step(plurality(snapshot = "step"), 4L)@snapshots, 2L)
  expect_length(fw_step(plurality(snapshot = "all"), 4L)@snapshots, 5L)
  # Which means composition is how the granularity is chosen.
  st <- plurality(snapshot = "step")
  expect_length(fw_step(fw_step(st, 1L), 1L)@snapshots, 3L)
})

test_that("the first snapshot is the starting mixture", {
  # A snapshot sequence that began at the first step could not show where the
  # fit started from.
  st <- plurality(snapshot = "step")
  expect_length(st@snapshots, 1L)
  first <- st@snapshots[[1L]]
  expect_identical(first$phase, "init")
  expect_identical(first$mixing, st@mixing)
  expect_identical(first$part, st@part)
  expect_identical(first$step, 0L)
})

# --- Finishing ----------------------------------------------------------------

test_that("ripr_finish returns a ripr_fit that keeps its state", {
  st <- fw_step(plurality(), 3L)
  fit <- ripr_finish(st)
  # By class name: the class object itself is not exported.
  expect_s3_class(fit, "ripr::ripr_fit")
  expect_true(S7::S7_inherits(fit@W0, finite_dist))
  expect_true(S7::S7_inherits(fit@P_star, mixture))
  expect_identical(fit@part, sort(st@part))
  expect_identical(fit@rounds, 0L)
  # The history is the state's, reached through it rather than copied out.
  expect_identical(fit@state, st)
  expect_identical(fit@state@trace, st@trace)

  out <- paste(capture.output(print(fit)), collapse = "\n")
  expect_match(out, "<ripr_fit>", fixed = TRUE)
  expect_match(out, paste("mixed over", n_atoms(fit@W0), "atoms"), fixed = TRUE)
  expect_match(out, "refined  no", fixed = TRUE)
  expect_match(out, "4 rows, 0 snapshots", fixed = TRUE)
  expect_no_match(out, "@")
  expect_match(format(fit), paste0("^ripr_fit: ", n_atoms(fit@W0), " atoms, KL "))
})

test_that("reoptimising in ripr_finish solves the weights to optimality", {
  # Over the fitted atoms, the optimal weights satisfy G <= 1 at every atom,
  # with equality wherever an atom keeps weight; and KL can only fall.
  st <- fw_step(plurality(), 6L)
  fit <- ripr_finish(st, reoptimise = TRUE)
  g <- exact_g(st, atoms(st@mixing), W = fit@W0)
  expect_lt(max(g), 1 + 1e-6)
  expect_equal(exact_g(st, atoms(fit@W0), W = fit@W0), rep(1, n_atoms(fit@W0)), tolerance = 1e-6)
  expect_lte(fit@kl, exact_kl(st))
})

test_that("ripr_finish is a pure conversion by default", {
  # Refining is opt-in: what comes back is what was fitted.
  st <- fw_step(plurality(), 5L)
  fit <- ripr_finish(st)
  expect_equal(fit@kl, utils::tail(st@trace$kl, 1L))
  expect_equal(fit@kl, exact_kl(st), tolerance = rounding_tol(1))
  # The same atoms and weights, reported grouped by part.
  by_part <- order(st@part)
  expect_equal(unname(weights(fit@W0)), unname(weights(st@mixing))[by_part])
  expect_equal(fit@W0@atoms, unname(atoms(st@mixing))[by_part, ])
  expect_identical(fit@part, sort(st@part))
})

test_that("kl describes the returned mixture, not the state", {
  # Refining and pruning both change the mixture, so a `kl` taken from the state
  # would describe something the caller was not given.
  st <- fw_step(plurality(), 5L)
  fit <- ripr_finish(st, reoptimise = TRUE, identify = TRUE)
  expect_equal(fit@kl, exact_kl(st, fit@W0), tolerance = rounding_tol(1))
})

test_that("refining lowers KL and shrinks the support", {
  st <- fw_step(plurality(), 6L)
  plain <- ripr_finish(st)
  refined <- ripr_finish(st, reoptimise = TRUE, identify = TRUE)
  expect_true(refined@kl <= plain@kl)
  expect_true(n_atoms(refined@W0) <= n_atoms(plain@W0))
  expect_true(refined@rounds >= 1L)
})

test_that("threshold keeps, drops, or errors", {
  st <- fw_step(plurality(), 5L)
  identified <- ripr_finish(st, identify = TRUE)
  expect_true(
    n_atoms(ripr_finish(st, identify = TRUE, threshold = -1)@W0) >=
      n_atoms(identified@W0)
  )
  expect_true(
    n_atoms(ripr_finish(st, threshold = 0.2)@W0) <= n_atoms(identified@W0)
  )
  expect_error(ripr_finish(st, threshold = 1), "`threshold` must be below 1")
  expect_error(ripr_finish(st, threshold = 0.99), "no atom has weight above")
})

test_that("the returned mixture is a distribution", {
  fit <- ripr_finish(fw_step(plurality(), 4L))
  expect_equal(sum(weights(fit@W0)), 1, tolerance = rounding_tol(1))
  expect_equal(
    sum(exp(log_density(
      fit@P_star,
      enumerate_space(fit@P_star@family@sample_space)
    ))),
    1,
    tolerance = rounding_tol(1)
  )
})

test_that("gap_final measures the returned mixture, gap_fit the state", {
  # `gap_fit` cannot show what finishing did -- it was recorded on the way in.
  st <- fw_step(plurality(), 5L)
  fit <- ripr_finish(st, reoptimise = TRUE, identify = TRUE, record_gap = TRUE)
  gaps <- st@trace$gap_after[!is.na(st@trace$gap_after)]
  expect_equal(fit@gap_fit, utils::tail(gaps, 1L))
  expect_true(is.na(utils::tail(st@trace$gap_after, 1L)))
  expect_true(!is.na(fit@gap_final))
  expect_true(is.na(ripr_finish(st)@gap_final))
})

# --- The one problem with a known answer --------------------------------------

test_that("the binomial projection concentrates at the boundary", {
  # H_0: p <= 1/2 against Bin(10, 0.75). The projection is the point mass at
  # p = 1/2, so nearly all the weight belongs on the tie.
  fit <- binomial(p = 0.75, n_seeds = 100L) |>
    fw_step(8L) |>
    em_step(20L) |>
    ripr_finish(reoptimise = TRUE, identify = TRUE)
  heavy <- fit@W0@atoms[which.max(weights(fit@W0)), ]
  expect_equal(heavy, c(0.5, 0.5), tolerance = 1e-3)
  expect_true(max(weights(fit@W0)) > 0.99)
})


# --- A null that is convex but not a simplex ----------------------------------

test_that("a polytope null fits and certifies end to end", {
  # Nothing the caller writes here is a simplex: the null is the polytope
  # `{theta_1 <= 1/2, theta_2 <= 1/2}` in the 2-simplex, stated as its four
  # vertices. `cells()` cuts it into two triangles.
  set.seed(11)
  fam <- multinomial_family(n_trials = 8L, k = 3L)
  square <- polytope_region(
    vertices = rbind(
      c(0.5, 0.5, 0),
      c(0, 0.5, 0.5),
      c(0, 0, 1),
      c(0.5, 0, 0.5)
    )
  )
  null <- null_model(fam, list(square))
  Q <- mixture(fam, dirac(theta = c(0.7, 0.2, 0.1)))

  fit <- ripr_init(
    Q,
    null,
    exact_engine(),
    control = ripr_control(n_seeds = 50L, n_restarts = 5L)
  ) |>
    fw_step(8L) |>
    em_step(8L) |>
    weight_step(8L) |>
    ripr_finish(record_gap = TRUE)

  # Every atom is in the null, and filed under the only part there is.
  expect_true(all(fit@part == 1L))
  for (j in seq_len(nrow(fit@W0@atoms))) {
    expect_true(contains(square, fit@W0@atoms[j, ], tol = 1e-6))
  }

  X <- likelihood(Q) / likelihood(fit@P_star)
  cert <- certify(X, null, tol = 1e-9)
  expect_true(all(cert@converged))
  expect_gte(cert@sup_ub, cert@sup_lb)

  # The same identity the simplex nulls satisfy: the certified bound lands at
  # `1 + gap`, with `gap` the duality gap the fit stopped on. Nothing about the
  # decomposition disturbs it.
  expect_equal(cert@sup_ub - 1, fit@gap_final, tolerance = 1e-4)
})
