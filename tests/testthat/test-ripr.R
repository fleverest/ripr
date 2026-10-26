# Properties of ripr(), the one-call front door.
#
# It is shorthand for a verb pipeline, so what is tested is that it runs the
# pipeline it documents, stops on the gap or the budget (and warns on the
# latter), and resumes a state or fit where it left off.

front_door_problem <- function(n = 8L, q = c(0.4, 0.35, 0.25)) {
  fam <- multinomial_family(n_trials = n, k = 3L)
  parts <- lapply(2:3, function(j) {
    vertices <- diag(3)
    vertices[1L, ] <- replace(numeric(3), c(1L, j), 0.5)
    simplex_region(vertices = vertices)
  })
  list(Q = fam(q), null = null_model(fam, parts))
}

quick <- function(...) ripr_control(n_seeds = 30L, n_restarts = 4L, ...)

test_that("each method stops on the gap and returns a finished fit", {
  p <- front_door_problem()
  for (method in c("sfw", "pdap", "fcfw")) {
    set.seed(1)
    expect_no_warning(
      fit <- ripr(p$Q, p$null, method, tol = 1e-6, control = quick())
    )
    expect_s3_class(fit, "ripr::ripr_fit")
    expect_lt(fit@gap_final, 1e-6)
    expect_lt(nrow(fit@state@trace) - 1L, 100L)
    expect_true(all(
      fit@state@trace$phase[-1L] ==
        c(
          sfw = "sfw",
          pdap = "fw",
          fcfw = "fw"
        )[[method]]
    ))
  }
})

test_that("ripr() is the pipeline it documents", {
  p <- front_door_problem()
  set.seed(2)
  fit <- ripr(p$Q, p$null, tol = 1e-8, control = quick())
  set.seed(2)
  by_hand <- ripr_init(p$Q, p$null, control = quick()) |>
    sfw_step(100L, until = gap_below(1e-8)) |>
    ripr_finish(record_gap = TRUE)
  expect_identical(fit@W0, by_hand@W0)
  expect_identical(fit@gap_final, by_hand@gap_final)
  expect_identical(fit@state@trace$kl, by_hand@state@trace$kl)
})

test_that("running out of steps warns, and a fit resumes where it stopped", {
  p <- front_door_problem()
  set.seed(3)
  expect_warning(
    short <- ripr(
      p$Q,
      p$null,
      "fcfw",
      tol = 1e-9,
      max_steps = 1L,
      control = quick()
    ),
    "after 1 step,",
    class = "ripr_budget_warning"
  )
  expect_identical(nrow(short@state@trace), 2L)
  expect_gte(short@gap_final, 1e-9)

  # The trace carries on from the fit's, and `max_steps` counts new steps.
  resumed <- ripr(short, method = "fcfw", tol = 1e-9, max_steps = 60L)
  tr <- resumed@state@trace
  expect_identical(tr$kl[1:2], short@state@trace$kl)
  expect_lte(nrow(tr), 62L)
  expect_lt(resumed@gap_final, 1e-9)
  # A state resumes the same way.
  expect_no_warning(
    ripr(resumed@state, method = "fcfw", tol = 1e-6, max_steps = 1L)
  )
})

test_that("resuming keeps the state's null and engine, but takes new control", {
  p <- front_door_problem()
  set.seed(4)
  st <- ripr_init(p$Q, p$null, control = quick())
  expect_error(ripr(st, p$null), "keeps its own `null` and `engine`")
  expect_error(ripr(st, engine = exact_engine()), "keeps its own")
  expect_error(ripr(st, "pdap"), 'method = "pdap"', fixed = TRUE)
  fit <- suppressWarnings(
    ripr(st, max_steps = 1L, control = quick(snapshot = "all"))
  )
  expect_identical(fit@state@control$snapshot, "all")
  expect_length(fit@state@snapshots, 1L)
})

test_that("ripr() refuses what it cannot fit", {
  p <- front_door_problem()
  expect_error(ripr(p$Q), "`null` is required")
  expect_error(ripr(p$null, p$null), "must be the alternative")
  expect_error(ripr(p$Q, p$null, method = "em"), "must be one of")
  expect_error(ripr(p$Q, p$null, max_steps = 0), "whole number")
  expect_error(ripr(p$Q, p$null, tol = -1), "larger than or equal to 0")
})
