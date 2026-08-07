# Properties of R/control.R.
#
# `ripr_control()` holds computational knobs only: nothing in it may change what
# the fit converges to, only how well it gets there. The tests below are mostly
# about validation, because these values are read deep inside oracle calls and a
# bad one surfaces far from its cause.

# --- Shape --------------------------------------------------------------------

test_that("counts are coerced to integer", {
  # `seq_len()` and friends downstream want integers; a double that happens to
  # be whole would work until it did not.
  ctl <- ripr_control(n_seeds = 10, n_restarts = 2, lb_fc_max_iter = 5)
  expect_type(ctl$n_seeds, "integer")
  expect_type(ctl$n_restarts, "integer")
  expect_type(ctl$lb_fc_max_iter, "integer")
  expect_type(ctl$lb_fc_tol, "double")
})

# --- snapshot -----------------------------------------------------------------

test_that("snapshot defaults to none and accepts the three levels", {
  expect_identical(ripr_control()$snapshot, "none")
  for (level in c("none", "step", "all")) {
    expect_identical(ripr_control(snapshot = level)$snapshot, level)
  }
})

test_that("snapshot rejects anything else, NULL included", {
  # Base `match.arg()` silently maps NULL to the first choice, which would turn
  # an accidental NULL into "record nothing" -- the wrong way to fail for a
  # deliberate diagnostic setting. `rlang::arg_match()` errors instead.
  expect_error(ripr_control(snapshot = "outer"), "must be one of")
  expect_error(ripr_control(snapshot = NULL), "must be")
  expect_error(ripr_control(snapshot = 5), "must be")
  expect_error(ripr_control(snapshot = c("none", "all")), "must be")
})

test_that("a near-miss level suggests the intended one", {
  expect_error(ripr_control(snapshot = "setp"), "Did you mean")
})

# --- Numeric validation -------------------------------------------------------

test_that("counts must be whole numbers", {
  # `as.integer()` truncates, so 1.99 would silently become 1.
  expect_error(ripr_control(n_seeds = 2.7), "whole number")
  expect_error(ripr_control(n_restarts = 1.99), "whole number")
  expect_error(ripr_control(lb_fc_max_iter = 10.5), "whole number")
})

test_that("counts respect their lower bounds", {
  # `n_seeds = 0` is allowed: the current atoms are always seeded separately, so
  # a search with no random starts is still a search. Refining zero of them is
  # not, so `n_restarts` starts at 1.
  expect_silent(ripr_control(n_seeds = 0))
  expect_error(ripr_control(n_seeds = -1), "whole number")
  expect_error(ripr_control(n_restarts = 0), "whole number")
  expect_error(ripr_control(lb_fc_max_iter = 0), "whole number")
})

test_that("the tolerance must be a single non-negative number", {
  expect_silent(ripr_control(lb_fc_tol = 0))
  expect_error(ripr_control(lb_fc_tol = -1), "larger than or equal to 0")
  expect_error(ripr_control(lb_fc_tol = c(1e-8, 1e-9)), "must be a number")
})

test_that("non-numeric input is rejected", {
  expect_error(ripr_control(n_seeds = "10"), "whole number")
  expect_error(ripr_control(n_seeds = NA), "whole number")
  expect_error(ripr_control(lb_fc_tol = "small"), "must be a number")
})
