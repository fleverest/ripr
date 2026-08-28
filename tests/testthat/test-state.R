# Properties of R/state.R.
#
# The state carries atoms and weights as lists indexed by part, and every
# piece of arithmetic works on flat vectors instead. Almost everything here is
# about that boundary holding: the flat views must agree with each other and
# invert cleanly, or an atom silently changes which chart it is searched in.

fixture <- function(snapshot = "none") {
  fam <- multinomial_family(n_trials = 8, k = 3)
  alternative <- induced_distribution(fam, point_mixing(c(0.5, 0.3, 0.2)))
  parts <- lapply(2:3, function(j) {
    simplex_region(
      vertices = cbind(
        c(0, 1, 0),
        c(0, 0, 1),
        replace(numeric(3), c(1L, j), 0.5)
      )
    )
  })
  ripr_state(
    atoms = list(
      cbind(c(0.40, 0.35, 0.25), c(0.30, 0.45, 0.25)),
      cbind(c(0.35, 0.25, 0.40))
    ),
    weights = list(c(0.5, 0.2), 0.3),
    alternative = alternative,
    null = null_model(fam, parts),
    engine = resolve_engine(exact_engine(), alternative, fam),
    control = ripr_control(snapshot = snapshot),
    trace = empty_trace(),
    snapshots = list(),
    iters = c(fw = 0L, lb = 0L, em = 0L, weight = 0L)
  )
}

# --- Validator ----------------------------------------------------------------

test_that("atoms and weights must agree, block by block", {
  # A parallel index vector would make this an invariant to remember; the list
  # makes it structural, and the validator makes a mismatch impossible.
  st <- fixture()
  expect_error(
    S7::set_props(st, weights = list(c(0.5, 0.2))),
    "one element per element of"
  )
  expect_error(
    S7::set_props(st, weights = list(c(0.3, 0.2, 0.2), 0.3)),
    "match the column count"
  )
})

test_that("there must be one atom block per part", {
  st <- fixture()
  expect_error(
    S7::set_props(
      st,
      atoms = st@atoms[1],
      weights = list(c(0.5, 0.5))
    ),
    "one element per part"
  )
})

test_that("weights are normalised across the whole list, not within a block", {
  # Each block sums to less than one; only the total is constrained.
  st <- fixture()
  expect_equal(sum(unlist(st@weights)), 1)
  expect_error(
    S7::set_props(st, weights = list(c(0.5, 0.2), 0.5)),
    "sum to 1"
  )
})

test_that("iters must be named", {
  # Named rather than a fixed set, so adding a verb later does not touch the
  # class. The names are what `bump()` and `record()` index by.
  st <- fixture()
  expect_error(
    S7::set_props(st, iters = c(0L, 0L, 0L, 0L)),
    "named integer vector"
  )
})

# --- Flat views ---------------------------------------------------------------

test_that("the flat views agree with one another", {
  st <- fixture()
  expect_identical(ncol(flat_atoms(st)), 3L)
  expect_length(flat_weights(st), 3L)
  expect_identical(flat_part(st), c(1L, 1L, 2L))
})

test_that("flat ordering is by part block, not chronological", {
  # An atom appended to part 1 lands before every atom of part 2, however
  # late it was added. Anything reading the flat order as a history is wrong.
  st <- fixture()
  st <- add_atom(st, c(0.45, 0.30, 0.25), 1L, c(0.4, 0.2, 0.1, 0.3))
  expect_identical(flat_part(st), c(1L, 1L, 1L, 2L))
  expect_equal(flat_atoms(st)[, 3L], c(0.45, 0.30, 0.25))
})

test_that("insert_index names the slot add_atom will use", {
  # The step layer inserts the candidate at this index from the outset, which is
  # what saves it from re-indexing the weights afterwards.
  st <- fixture()
  expect_identical(insert_index(st, 1L), 3L)
  expect_identical(insert_index(st, 2L), 4L)

  at <- insert_index(st, 1L)
  st2 <- add_atom(st, c(0.45, 0.30, 0.25), 1L, c(0.4, 0.2, 0.1, 0.3))
  expect_equal(flat_atoms(st2)[, at], c(0.45, 0.30, 0.25))
  expect_equal(flat_weights(st2)[at], 0.1)
})

test_that("flattening and unflattening round-trip", {
  st <- fixture()
  expect_identical(unflatten_weights(st, flat_weights(st)), st@weights)
  expect_equal(unflatten_atoms(st, flat_atoms(st)), st@atoms)
})

test_that("empty blocks survive the round trip", {
  # A part with no atoms must stay in the list as a zero-column matrix, or
  # the block-to-part correspondence shifts.
  expect_identical(
    split_by_sizes(c(0.4, 0.6), c(1L, 0L, 1L)),
    list(0.4, numeric(0), 0.6)
  )
})

test_that("set_weights leaves the atoms alone", {
  st <- fixture()
  st2 <- set_weights(st, c(0.2, 0.3, 0.5))
  expect_identical(st2@atoms, st@atoms)
  expect_equal(flat_weights(st2), c(0.2, 0.3, 0.5))
})

# --- Core quantities ----------------------------------------------------------

test_that("KL is zero when the mixture is the alternative", {
  # The one value of KL that is known in closed form without integrating
  # anything, so it catches a mis-signed or mis-weighted reduction.
  fam <- multinomial_family(n_trials = 8, k = 3)
  theta <- c(0.5, 0.3, 0.2)
  alternative <- induced_distribution(fam, point_mixing(theta))
  st <- ripr_state(
    atoms = list(matrix(theta, ncol = 1L)),
    weights = list(1),
    alternative = alternative,
    null = null_model(fam, list(simplex_region(vertices = diag(3)))),
    engine = resolve_engine(exact_engine(), alternative, fam),
    control = ripr_control(),
    trace = empty_trace(),
    snapshots = list(),
    iters = c(fw = 0L, lb = 0L, em = 0L, weight = 0L)
  )
  expect_equal(kl_divergence(st), 0, tolerance = rounding_tol(1))
})

test_that("kl_divergence accepts a precomputed log_p", {
  # Recomputing log_p is an O(MC) reduction.
  st <- fixture()
  ld <- compile_engine(st@engine)
  expect_equal(
    kl_divergence(st, log_p = log_p_at_nodes(st, ld)),
    kl_divergence(st)
  )
})

# --- Trace --------------------------------------------------------------------

test_that("a trace row carries every step counter", {
  # So `kl` can be plotted against any one of them, or against their sum.
  st <- fixture()
  st <- bump(st, "fw")
  st <- bump(st, "em", by = 3L)
  st <- record(st, phase = "fw", kl = 0.25)
  row <- st@trace[1L, ]
  expect_identical(c(row$fw, row$lb, row$em, row$weight), c(1L, 0L, 3L, 0L))
  expect_identical(row$phase, "fw")
  expect_equal(row$kl, 0.25)
})

test_that("record derives support size and max weight from the state", {
  st <- record(fixture(), phase = "init", kl = 0.1)
  expect_identical(st@trace$support_size, 3L)
  expect_equal(st@trace$max_weight, 0.5)
})

test_that("record requires kl", {
  # The old default silently recomputed it, hiding an O(MC) computation
  expect_error(record(fixture(), phase = "fw"), "kl")
})

test_that("unfilled columns come back as NA of the right type", {
  st <- record(fixture(), phase = "em", kl = 0.1)
  expect_identical(st@trace$part, NA_integer_)
  expect_identical(st@trace$direction, NA_character_)
  expect_identical(st@trace$gap, NA_real_)
})

# --- Snapshots ----------------------------------------------------------------

test_that("wants_snapshot distinguishes per-call from per-iteration", {
  # "step" counts calls, so it fires only on the last iteration of one; "all"
  # counts iterations. `record()` cannot decide this itself -- it does not know
  # where in a `times` loop it sits.
  expect_false(wants_snapshot(fixture("none"), last = TRUE))
  expect_true(wants_snapshot(fixture("step"), last = TRUE))
  expect_false(wants_snapshot(fixture("step"), last = FALSE))
  expect_true(wants_snapshot(fixture("all"), last = FALSE))
})

test_that("a snapshot keeps the mixture and the counters", {
  st <- bump(fixture(), "fw", by = 2L)
  st <- snapshot_state(st, "fw")
  snap <- st@snapshots[[1L]]
  expect_identical(snap$iters[["fw"]], 2L)
  expect_identical(snap$atoms, st@atoms)
  expect_identical(snap$weights, st@weights)
})

test_that("record does not snapshot", {
  # The two differ in frequency within one call and in cost by orders of
  # magnitude, so the verb owns the decision.
  st <- record(fixture("all"), phase = "fw", kl = 0.1)
  expect_length(st@snapshots, 0L)
})
