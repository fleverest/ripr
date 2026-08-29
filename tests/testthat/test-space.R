# Properties of R/space.R and R/count_space.R.
#
# A space answers three questions and nothing else: how wide one outcome
# is, whether a given object is one, and -- when it can -- what all of them are.
# Two things are worth guarding. That membership is refused precisely, since
# every random variable's input check routes through here and a silently
# reshaped outcome returns a number rather than a complaint. And that spaces
# compare by value, since that is what lets two unrelated families over the same
# space share random variables.

# --- Dimension ----------------------------------------------------------------

test_that("spaces refuse nonsensical shapes at construction", {
  expect_error(count_space(n = -1L, k = 3L), "non-negative integer")
  expect_error(count_space(n = 4L, k = 0L), ">= 1")
  expect_error(count_space(n = c(1L, 2L), k = 3L), "single")
  expect_error(real_region(0L), "positive integer")
  expect_error(real_region(c(1L, 2L)), "single")
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
  space <- real_region(2L)
  expect_false(is_finite_space(space))
  expect_error(enumerate_space(space), "cannot be enumerated")
})

# --- Membership: shape checks shared by every space ---------------------------

test_that("one outcome may be a vector, and many a matrix", {
  space <- real_region(3L)
  expect_identical(dim(validate_outcome(space, c(1, 2, 3))), c(1L, 3L))
  x <- rbind(c(1, 2, 3), c(4, 5, 6))
  expect_identical(validate_outcome(space, x), x)
})

test_that("the wrong shape is named, not reshaped", {
  space <- real_region(3L)
  expect_error(validate_outcome(space, c(1, 2)), "length-3 vector")
  expect_error(validate_outcome(space, matrix(1, 2L, 4L)), "3 columns")
})

test_that("outcomes must be numeric and present", {
  space <- real_region(3L)
  expect_error(validate_outcome(space, c("a", "b", "c")), "must be numeric")
  expect_error(validate_outcome(space, c(1, NA, 3)), "must not be missing")
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
  space <- real_region(2L)
  expect_silent(validate_outcome(space, c(-3, 40)))
  expect_error(validate_outcome(space, c(Inf, 0)), "finite")
})

# --- Identity -----------------------------------------------------------------

test_that("spaces built separately from the same description are identical", {
  # What makes a random variable from one family combinable with one from
  # another over the same space: `shared_space()` compares by value.
  expect_identical(count_space(n = 20L, k = 3L), count_space(n = 20L, k = 3L))
  expect_identical(count_space(20, 3), count_space(n = 20L, k = 3L))
  expect_identical(real_region(2L), real_region(2))
  expect_false(identical(count_space(20, 3), count_space(20, 4)))
  expect_false(identical(real_region(1L), count_space(1, 1)))
})


# --- The shared `space` parent ------------------------------------------------

test_that("both branches of the hierarchy are spaces", {
  # Both branches hang directly off it: the finite, listable one and the
  # geometric one. There is no longer an abstract `sample_space` between them.
  expect_identical(count_space@parent, space)
  expect_identical(region@parent, space)
})

test_that("a count space answers contains(), agreeing with validate_outcome", {
  s <- count_space(n = 4L, k = 3L)
  expect_true(contains(s, c(2L, 1L, 1L)))
  expect_true(contains(s, c(4L, 0L, 0L)))

  expect_false(contains(s, c(2L, 1L, 0L))) # sums to 3, not 4
  expect_false(contains(s, c(5L, 0L, -1L))) # negative
  expect_false(contains(s, c(2.5, 1, 0.5))) # not whole
  expect_false(contains(s, c(2L, 2L))) # wrong length

  # `contains()` and `validate_outcome()` restate the same three conditions
  # rather than sharing one implementation: the predicate returns a single
  # boolean, so layering it under validate_outcome would cost the specific
  # messages that say *which* condition failed. The duplication is only safe
  # while the two agree, so that is pinned here in both directions.
  agrees <- function(x) {
    valid <- !inherits(try(validate_outcome(s, x), silent = TRUE), "try-error")
    contains(s, x) == valid
  }
  expect_true(all(apply(enumerate_space(s), 1L, agrees)))
  for (bad in list(
    c(2L, 1L, 0L),      # sums to 3, not 4
    c(5L, 0L, -1L),     # negative
    c(2.5, 1, 0.5),     # not whole
    c(0L, 0L, 0L),      # sums to 0
    c(4L, 1L, -1L)      # sums correctly, but negative
  )) {
    expect_true(agrees(bad))
  }
})


test_that("a region can serve as a sample space", {
  # The property is typed `space`, not `sample_space`, so what a family draws
  # outcomes from and what it takes parameters in are the same kind of thing.
  # Nothing yet relies on this; it is what lets `real_region` be retired.
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  swapped <- parametric_family(
    sample_space = simplex_region(vertices = diag(3)),
    parameter_space = fam@parameter_space
  )
  expect_true(S7::S7_inherits(swapped@sample_space, region))
  expect_equal(space_dim(swapped@sample_space), 3L)
})
