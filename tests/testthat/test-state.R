# Properties of R/state.R.
#
# The state carries its mixture as a `finite_dist` with a part per atom, so
# almost everything here is about the two staying aligned, or an atom silently
# changes which chart it is searched in. The constructor is not exported, so
# every state here comes from `ripr_init()` and the verbs, and what is tested
# is what a caller can see: the validator, the trace, the snapshots, the print.

fixture <- function(snapshot = "none", record_gap = FALSE) {
  fam <- multinomial_family(n_trials = 8, k = 3)
  alternative <- mixture(fam, dirac(c(0.5, 0.3, 0.2)))
  parts <- lapply(2:3, function(j) {
    simplex_region(
      vertices = rbind(
        c(0, 1, 0),
        c(0, 0, 1),
        replace(numeric(3), c(1L, j), 0.5)
      )
    )
  })
  # Two atoms in the first part, `theta_1 <= theta_2`, and one in the second,
  # `theta_1 <= theta_3`. Neither part's weights sum to one on its own.
  ripr_init(
    alternative,
    null_model(fam, parts),
    atoms = list(
      rbind(c(0.30, 0.45, 0.25), c(0.35, 0.40, 0.25)),
      rbind(c(0.30, 0.25, 0.45))
    ),
    weights = list(c(0.5, 0.2), 0.3),
    record_gap = record_gap,
    control = ripr_control(n_seeds = 20L, n_restarts = 3L, snapshot = snapshot)
  )
}

# --- Validator ----------------------------------------------------------------

test_that("the mixture and its parts must agree", {
  st <- fixture()
  expect_identical(st@part, c(1L, 1L, 2L))
  expect_error(st@part <- c(1L, 2L), "one element per atom")
  # Replacing the mixture alone, when the atom count changes, is refused.
  expect_error(
    st@mixing <- dirac(c(0.4, 0.35, 0.25)),
    "one element per atom"
  )
  # And a weight for every atom, whichever part it is in.
  expect_error(
    ripr_init(
      st@alternative,
      st@null,
      atoms = list(atoms(st@mixing)[1:2, ], atoms(st@mixing)[3L, ]),
      weights = list(0.5, 0.5)
    ),
    "one entry per row"
  )
})

test_that("every atom's part is a part of the null", {
  st <- fixture()
  expect_error(st@part <- c(1L, 1L, 3L), "index the parts")
  expect_error(st@part <- c(1L, NA, 2L), "index the parts")
})

test_that("weights are normalised across all atoms, not within a part", {
  # Each part's weights sum to less than one; only the total is constrained.
  st <- fixture()
  expect_equal(weights(st@mixing), c(0.5, 0.2, 0.3))
  again <- function(w) {
    ripr_init(
      st@alternative,
      st@null,
      atoms = list(atoms(st@mixing)[1:2, ], atoms(st@mixing)[3L, ]),
      weights = w
    )
  }
  expect_error(again(list(c(0.5, 0.2), 0.5)), "sum to 1")
  expect_error(again(list(c(1.2, -0.2), 0)), "non-negative")
})

test_that("the oracle cache is empty or a complete oracle result", {
  # `fw_step()` reads `value`, `theta`, `part` and `elapsed` off it without
  # checking, so a partial cache must be impossible to build.
  st <- fixture()
  expect_null(st@oracle)
  expect_error(
    st@oracle <- list(value = 1.2),
    "list\\(value, theta, part, elapsed\\)"
  )
  full <- list(value = 1.2, theta = c(0.4, 0.3, 0.3), part = 1L, elapsed = 0)
  st@oracle <- full
  expect_identical(st@oracle, full)
  # What a sweep caches is complete by construction.
  swept <- fixture(record_gap = TRUE)@oracle
  expect_true(all(c("value", "theta", "part", "elapsed") %in% names(swept)))
})

# --- Adding atoms -------------------------------------------------------------

test_that("a new atom is appended last, with its part", {
  # The step layer carries its candidate as the last row for this reason:
  # the weights it computes need no re-indexing when written back.
  set.seed(1)
  st <- fixture()
  stepped <- fw_step(st, 1L)
  row <- utils::tail(stepped@trace, 1L)
  expect_false(is.na(row$part))
  n <- n_atoms(st@mixing)
  expect_identical(n_atoms(stepped@mixing), n + 1L)
  expect_identical(stepped@part, c(st@part, row$part))
  expect_identical(atoms(stepped@mixing)[n + 1L, ], row$oracle_theta[[1L]])
  # The incumbents are rescaled, not moved or reordered.
  expect_identical(atoms(stepped@mixing)[seq_len(n), ], atoms(st@mixing))
  expect_equal(
    weights(stepped@mixing)[seq_len(n)],
    weights(st@mixing) * (1 - row$step_size)
  )
  expect_equal(weights(stepped@mixing)[n + 1L], row$step_size)
})

# --- Core quantities ----------------------------------------------------------

test_that("KL is zero when the mixture is the alternative", {
  # The one value of KL that is known in closed form without integrating
  # anything, so it catches a mis-signed or mis-weighted reduction. The
  # default start projects the alternative's own parameter, which is already
  # in the null.
  fam <- multinomial_family(n_trials = 8, k = 3)
  theta <- c(0.5, 0.3, 0.2)
  st <- ripr_init(
    mixture(fam, dirac(theta)),
    null_model(fam, list(simplex_region(vertices = diag(3))))
  )
  expect_equal(atoms(st@mixing)[1L, ], theta)
  expect_equal(st@trace$kl, 0, tolerance = rounding_tol(1))
})

test_that("the init row's KL is the starting mixture's", {
  st <- fixture()
  expect_equal(st@trace$kl, exact_kl(st), tolerance = rounding_tol(1))
})

# --- Trace --------------------------------------------------------------------

test_that("a trace row carries its step index and phase", {
  # One index for every verb, counting from the init row; a per-verb count is
  # a count over `phase`, so no counter can drift out of step with the rows.
  st <- em_step(fw_step(fixture(), 1L), 1L)
  expect_identical(st@trace$step, 0:2)
  expect_identical(st@trace$phase, c("init", "fw", "em"))
  expect_identical(verb_counts(st), c(fw = 1L, lb = 0L, em = 1L, weight = 0L))
})

test_that("the trace is a read-only data frame over the stored rows", {
  # Rows are kept as a list, so appending one does not copy the others; the
  # data frame is what users and plotting code read.
  st <- fixture()
  expect_s3_class(st@trace, "data.frame")
  expect_identical(nrow(st@trace), 1L)
  stepped <- em_step(st, 2L)
  # Every row has every column, whatever the verb that wrote it filled in.
  expect_identical(
    vapply(st@trace, class, ""),
    vapply(stepped@trace, class, "")
  )
  expect_type(stepped@trace$oracle_theta, "list")
  expect_type(stepped@trace$gap_after_theta, "list")
  expect_error(stepped@trace <- data.frame(), "read-only")
  expect_length(stepped@trace_rows, 3L)
  expect_identical(nrow(stepped@trace), 3L)
  kl <- utils::tail(stepped@trace$kl, 1L)
  expect_identical(
    format(stepped),
    paste0("ripr_state: 3 atoms, KL ", format(signif(kl, 6L)))
  )
})

test_that("a new row clears the oracle cache", {
  # A new row means a new mixture, and the cache described the old one.
  st <- fixture(record_gap = TRUE)
  expect_false(is.null(st@oracle))
  expect_null(em_step(st, 1L)@oracle)
})

test_that("a row's gap is written once, from the oracle, never overwritten", {
  st <- fixture(record_gap = TRUE)
  row <- st@trace
  found <- st@oracle
  expect_identical(row$gap_after, found$value - 1)
  expect_identical(row$gap_after_part, found$part)
  expect_identical(row$gap_after_theta[[1L]], found$theta)
  expect_identical(row$gap_after_elapsed, found$elapsed)

  # With the cache dropped, the next step searches again and would fill the
  # same row a second time; the first measurement stands.
  st@oracle <- NULL
  set.seed(99)
  stepped <- fw_step(st, 1L)
  cols <- c(
    "gap_after",
    "gap_after_theta",
    "gap_after_part",
    "gap_after_elapsed"
  )
  expect_identical(stepped@trace[1L, cols], st@trace[1L, cols])
})

test_that("a row derives support size and max weight from the mixture", {
  st <- fixture()
  expect_identical(st@trace$support_size, 3L)
  expect_equal(st@trace$max_weight, 0.5)
  stepped <- weight_step(st, 1L)
  expect_equal(
    utils::tail(stepped@trace$max_weight, 1L),
    max(weights(stepped@mixing))
  )
})

test_that("columns a verb does not fill come back as NA of the right type", {
  st <- em_step(fixture(), 1L)
  row <- utils::tail(st@trace, 1L)
  expect_identical(row$part, NA_integer_)
  expect_identical(row$direction, NA_character_)
  expect_identical(row$step_size, NA_real_)
  expect_identical(row$oracle_value, NA_real_)
  expect_identical(row$gap_after, NA_real_)
  expect_identical(row$gap_after_part, NA_integer_)
  expect_identical(row$gap_after_elapsed, NA_real_)
  expect_true(is.na(row$oracle_theta[[1L]]))
  # The clock is always filled: every verb times its own step.
  expect_false(anyNA(st@trace$elapsed))
})

# --- Snapshots ----------------------------------------------------------------

test_that("a snapshot keeps the mixture and parts with its step", {
  st <- fw_step(fixture("all"), 1L)
  expect_length(st@snapshots, 2L)
  snap <- st@snapshots[[2L]]
  expect_named(snap, c("step", "phase", "mixing", "part"))
  expect_identical(snap$step, 1L)
  expect_identical(snap$phase, "fw")
  expect_identical(snap$mixing, st@mixing)
  expect_identical(snap$part, st@part)
})

# --- Printing -----------------------------------------------------------------

test_that("a state prints a summary, not a property dump", {
  st <- fw_step(fixture(), 1L)
  kl <- format(signif(utils::tail(st@trace$kl, 1L), 6L))
  out <- paste(capture.output(print(st)), collapse = "\n")
  expect_match(out, "<ripr_state>", fixed = TRUE)
  expect_match(out, paste0(n_atoms(st@mixing), " atoms, KL ", kl), fixed = TRUE)
  expect_match(out, "2 parts", fixed = TRUE)
  expect_match(out, "fw 1, lb 0, em 0, weight 0", fixed = TRUE)
  expect_match(out, "2 rows, 0 snapshots", fixed = TRUE)
  expect_no_match(out, "@")
  # The step's search filled the init row, and nothing has measured the
  # mixture the step made: that gap belongs to an earlier iterate.
  expect_match(out, "at an earlier iterate", fixed = TRUE)
  expect_match(
    paste(capture.output(print(fixture())), collapse = "\n"),
    "gap      none recorded",
    fixed = TRUE
  )

  # A cached oracle describes the current mixture and prints as current; once
  # another step clears it, the recorded gap is flagged as belonging to an
  # earlier iterate.
  measured <- fw_step(fixture(), 1L, record_gap = TRUE)
  gap <- format(signif(measured@oracle$value - 1, 3L))
  expect_match(
    paste(capture.output(print(measured)), collapse = "\n"),
    paste0("gap      ", gap, "\n"),
    fixed = TRUE
  )
  expect_match(
    paste(capture.output(print(em_step(measured, 1L))), collapse = "\n"),
    paste(gap, "at an earlier iterate"),
    fixed = TRUE
  )
})
