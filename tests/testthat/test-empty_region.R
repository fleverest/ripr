# Tests for `empty_region`, the zero of the set algebra.

simplex3 <- function() simplex_region(vertices = diag(3))

test_that("empty_region() carries a dimension and nothing else", {
  e <- empty_region(3L)
  expect_true(S7_inherits(e, region))
  expect_false(S7_inherits(e, convex_region))
  expect_identical(space_dim(e), 3L)

  expect_error(empty_region(0L), "positive integer")
  expect_error(empty_region(c(2L, 3L)), "positive integer")
  expect_error(empty_region("three"), "positive integer")
})


test_that("an empty region acts as a list of length zero", {
  e <- empty_region(3L)
  expect_identical(parts(e), list())
  expect_identical(cells(e), list())
  expect_identical(n_parts(e), 0L)
  expect_identical(n_cells(e), 0L)
  expect_identical(length(e), 0L)
  expect_identical(as.list(e), list())
  expect_error(e[[1L]], "subscript out of bounds")
  # Any subset of nothing is nothing.
  expect_identical(e[integer(0)], e)
  expect_identical(e[1L], e)
})


test_that("the empty region answers the predicates an empty set should", {
  e <- empty_region(3L)
  expect_true(is_empty(e))
  expect_true(is_bounded(e))
  expect_false(contains(e, c(1, 0, 0)))
})


test_that("an empty set has no representation to take", {
  e <- empty_region(2L)
  expect_error(h_rep(e), "not defined for an `empty_region`")
  expect_error(v_rep(e), "not defined for an `empty_region`")
  expect_error(q_hrep(e), "not defined for an `empty_region`")
  expect_error(q_vrep(e), "not defined for an `empty_region`")
})


test_that("union treats empty as its identity", {
  s <- simplex3()
  e <- empty_region(3L)
  expect_identical(union(e, s), s)
  expect_identical(union(s, e), s)
  expect_identical(union_region(e, e), e)
  expect_identical(union_region(e), e)

  # A mismatched dimension is still a mismatch, not silently dropped.
  expect_error(union(empty_region(2L), s), "same ambient dimension")
})


test_that("intersection absorbs to empty and setdiff subtracts nothing", {
  s <- simplex3()
  e <- empty_region(3L)

  expect_true(S7_inherits(intersect(s, e), empty_region))
  expect_true(S7_inherits(intersect(e, s), empty_region))
  expect_true(S7_inherits(setdiff(e, s), empty_region))
  # Subtracting nothing gives the set back.
  expect_true(setequal(setdiff(s, e), s))
  # A convex region minus itself is empty: every facet is ambient-implied,
  # so the whole ambient is covered and nothing remains.
  expect_true(S7_inherits(setdiff(s, s), empty_region))
})


test_that("setequal knows the empty set from an occupied one", {
  e <- empty_region(3L)
  expect_true(setequal(e, empty_region(3L)))
  expect_false(setequal(e, empty_region(2L)))
  expect_false(setequal(e, simplex3()))
  expect_false(setequal(simplex3(), e))
})


test_that("disjoin passes an empty region through", {
  e <- empty_region(3L)
  expect_identical(disjoin(e), e)
})


test_that("the algebra is closed through an empty intermediate", {
  s <- simplex3()
  nothing <- intersect(
    point_region(theta = c(1, 0, 0)),
    point_region(theta = c(0, 1, 0))
  )
  expect_true(S7_inherits(nothing, empty_region))
  # An empty result chains straight back into every verb.
  expect_identical(union(nothing, s), s)
  expect_true(S7_inherits(intersect(nothing, s), empty_region))
  expect_true(setequal(setdiff(s, nothing), s))
})


test_that("a null model refuses an empty region", {
  fam <- multinomial_family(n_trials = 2L, k = 3L)
  expect_error(
    null_model(fam, empty_region(3L)),
    "the null is empty"
  )
})


test_that("an empty region formats as what it is", {
  expect_match(format(empty_region(3L)), "empty region in dimension 3")
  expect_output(print(empty_region(3L)), "empty region in dimension 3")
})
