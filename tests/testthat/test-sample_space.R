# Properties of R/sample_space.R.
#
# A sample space answers three questions and nothing else: how wide one outcome
# is, whether a given object is one, and -- when it can -- what all of them are.
# Two things are worth guarding. That membership is refused precisely, since
# every random variable's input check routes through here and a silently
# reshaped outcome returns a number rather than a complaint. And that spaces
# compare by value, since that is what lets two unrelated families over the same
# space share random variables.

# --- Dimension ----------------------------------------------------------------

test_that("dimension is the width of one outcome", {
  expect_equal(space_dim(count_space(n = 4L, k = 3L)), 3L)
  expect_equal(space_dim(count_space(n = 0L, k = 1L)), 1L)
  expect_equal(space_dim(real_space(2L)), 2L)
  expect_type(space_dim(real_space(2L)), "integer")
})

test_that("spaces refuse nonsensical shapes at construction", {
  expect_error(count_space(n = -1L, k = 3L), "non-negative integer")
  expect_error(count_space(n = 4L, k = 0L), ">= 1")
  expect_error(count_space(n = c(1L, 2L), k = 3L), "single")
  expect_error(real_space(0L), "positive integer")
  expect_error(real_space(c(1L, 2L)), "single")
})

# --- Enumeration --------------------------------------------------------------

test_that("enumerate_counts produces every count vector exactly once", {
  for (n in c(0L, 1L, 5L)) {
    for (k in c(1L, 2L, 4L)) {
      x <- enumerate_counts(n, k)
      expect_equal(nrow(x), choose(n + k - 1, k - 1))
      expect_equal(ncol(x), k)
      expect_true(all(rowSums(x) == n))
      expect_true(all(x >= 0))
      expect_equal(nrow(unique(x)), nrow(x))
    }
  }
})

test_that("a count space enumerates itself", {
  space <- count_space(n = 3L, k = 2L)
  expect_true(is_finite_space(space))
  expect_equal(enumerate_space(space), enumerate_counts(3L, 2L))
  expect_equal(ncol(enumerate_space(space)), space_dim(space))
})

test_that("every enumerated outcome is an element of the space it came from", {
  # The two halves of the space -- what it lists, and what it accepts -- are
  # written independently, so nothing forces them to agree.
  space <- count_space(n = 5L, k = 3L)
  x <- enumerate_space(space)
  expect_identical(validate_outcome(space, x), x)
  for (i in seq_len(nrow(x))) {
    expect_silent(validate_outcome(space, x[i, ]))
  }
})

test_that("an infinite space says so, and refuses to be enumerated", {
  space <- real_space(2L)
  expect_false(is_finite_space(space))
  expect_error(enumerate_space(space), "cannot be enumerated")
  expect_error(enumerate_space(space), "real_space", fixed = TRUE)
})

# --- Membership: shape checks shared by every space ---------------------------

test_that("one outcome may be a vector, and many a matrix", {
  space <- real_space(3L)
  expect_identical(dim(validate_outcome(space, c(1, 2, 3))), c(1L, 3L))
  x <- rbind(c(1, 2, 3), c(4, 5, 6))
  expect_identical(validate_outcome(space, x), x)
})

test_that("the wrong shape is named, not reshaped", {
  space <- real_space(3L)
  expect_error(validate_outcome(space, c(1, 2)), "length-3 vector")
  expect_error(validate_outcome(space, c(1, 2, 3, 4)), "length-3 vector")
  expect_error(validate_outcome(space, matrix(1, 2L, 4L)), "3 columns")
})

test_that("outcomes must be numeric and present", {
  space <- real_space(3L)
  expect_error(validate_outcome(space, c("a", "b", "c")), "must be numeric")
  expect_error(validate_outcome(space, c(1, NA, 3)), "must not be missing")
  expect_error(validate_outcome(space, c(1, NaN, 3)), "must not be missing")
})

# --- Membership: what each space adds ----------------------------------------

test_that("counts must be whole, non-negative, and sum to the total", {
  # The message names the number rather than calling it `n_trials`: the space
  # cannot know what the family reading it calls that total.
  space <- count_space(n = 8L, k = 3L)
  expect_silent(validate_outcome(space, c(8, 0, 0)))
  expect_error(validate_outcome(space, c(4, 2, 1)), "summing to 8")
  expect_error(validate_outcome(space, c(9, -1, 0)), "non-negative whole")
  expect_error(validate_outcome(space, c(4.5, 2, 1.5)), "non-negative whole")
})

test_that("reals must be finite", {
  # An infinite outcome has zero density under every parameter here, so a
  # likelihood ratio there is `0 / 0`, and a `NaN` is a worse answer than a
  # complaint.
  space <- real_space(2L)
  expect_silent(validate_outcome(space, c(-3, 40)))
  expect_error(validate_outcome(space, c(Inf, 0)), "finite")
  expect_error(validate_outcome(space, c(0, -Inf)), "finite")
})

# --- Identity -----------------------------------------------------------------

test_that("spaces built separately from the same description are identical", {
  # What makes a random variable from one family combinable with one from
  # another over the same space: `shared_space()` compares by value.
  expect_identical(count_space(n = 20L, k = 3L), count_space(n = 20L, k = 3L))
  expect_identical(count_space(20, 3), count_space(n = 20L, k = 3L))
  expect_identical(real_space(2L), real_space(2))
  expect_false(identical(count_space(20, 3), count_space(20, 4)))
  expect_false(identical(real_space(1L), count_space(1, 1)))
})
