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
  ctl <- ripr_control(
    n_seeds = 10,
    n_restarts = 2,
    search_max_eval = 3,
    slide_max_eval = 4,
    weight_max_iter = 5
  )
  for (count in c(
    "n_seeds",
    "n_restarts",
    "search_max_eval",
    "slide_max_eval",
    "weight_max_iter"
  )) {
    expect_type(ctl[[count]], "integer")
  }
  expect_type(ctl$weight_slsqp_tol, "double")
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

# --- Numeric validation -------------------------------------------------------

test_that("counts must be whole numbers", {
  # `as.integer()` truncates, so 1.99 would silently become 1.
  expect_error(ripr_control(n_seeds = 2.7), "whole number")
  expect_error(ripr_control(n_restarts = 1.99), "whole number")
  expect_error(ripr_control(weight_max_iter = 10.5), "whole number")
  expect_error(ripr_control(search_max_eval = 1.5), "whole number")
  expect_error(ripr_control(slide_max_eval = 2.5), "whole number")
})

test_that("counts respect their lower bounds", {
  # `n_seeds = 0` is allowed: the current atoms are always seeded separately, so
  # a search with no random starts is still a search. Refining zero of them is
  # not, so `n_restarts` starts at 1.
  expect_silent(ripr_control(n_seeds = 0))
  expect_error(ripr_control(n_seeds = -1), "whole number")
  expect_error(ripr_control(n_restarts = 0), "whole number")
  expect_error(ripr_control(weight_max_iter = 0), "whole number")
  expect_error(ripr_control(search_max_eval = 0), "whole number")
  expect_error(ripr_control(slide_max_eval = 0), "whole number")
})

test_that("each tolerance must be a single non-negative number", {
  # The error names the offending setting, and ripr_control() as the caller.
  for (tol in c(
    "search_xtol",
    "slide_xtol",
    "weight_sqp_tol",
    "weight_slsqp_tol",
    "merge_atol",
    "merge_rtol",
    "near_tol"
  )) {
    expect_silent(do.call(ripr_control, stats::setNames(list(0), tol)))
    bad <- function(v) do.call(ripr_control, stats::setNames(list(v), tol))
    expect_error(bad(-1), paste0("`", tol, "` must be"))
    expect_error(bad(c(1e-8, 1e-9)), "must be a number")
  }
  expect_error(ripr_control(slide_xtol = -1), class = "rlang_error")
  err <- tryCatch(ripr_control(slide_xtol = -1), error = identity)
  expect_match(conditionMessage(err), "larger than or equal to 0")
  expect_identical(deparse(conditionCall(err)[[1L]]), "ripr_control")
})

test_that("non-numeric input is rejected", {
  expect_error(ripr_control(n_seeds = "10"), "whole number")
  expect_error(ripr_control(n_seeds = NA), "whole number")
  expect_error(ripr_control(weight_slsqp_tol = "small"), "must be a number")
})
