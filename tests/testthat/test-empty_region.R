# Tests for `empty_region`, the zero of the set algebra.

simplex3 <- function() simplex_region(vertices = diag(3))


test_that("an empty region has no parts and no cells", {
  e <- empty_region()
  expect_identical(parts(e), list())
  expect_identical(cells(e), list())
})


test_that("the empty region answers the predicates an empty set should", {
  e <- empty_region()
  expect_true(is_empty(e))
  expect_true(is_bounded(e))
  expect_false(contains(e, c(1, 0, 0)))
})


test_that("an empty set has no representation or dimension to take", {
  # No facets, no generators and no chart: there is nothing to describe. Nor
  # an ambient dimension, since the empty set is the same set in all of them.
  e <- empty_region()
  expect_error(e@facets, "Can't find property")
  expect_error(e@generators, "Can't find property")
  expect_error(chart(e), "Can't find method")
  expect_error(space_dim(e), "same set in every ambient dimension")
})


test_that("`|` treats empty as its identity, in any dimension", {
  s <- simplex3()
  e <- empty_region()
  expect_identical(e | s, s)
  expect_identical(s | e, s)
  expect_identical(union_region(e, e), e)
  expect_identical(union_region(e), e)

  # The same empty region is compatible with a region of any dimension.
  square <- polytope_region(
    vertices = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  )
  expect_identical(e | square, square)
})


test_that("`&` absorbs to empty and `-` subtracts nothing", {
  s <- simplex3()
  e <- empty_region()

  expect_true(S7_inherits(s & e, empty_region))
  expect_true(S7_inherits(e & s, empty_region))
  expect_true(S7_inherits(e - s, empty_region))
  # Subtracting nothing gives the set back, untouched.
  expect_identical(s - e, s)
  # A convex region minus itself is empty: every facet is ambient-implied,
  # so the whole ambient is covered and nothing remains.
  expect_true(S7_inherits(s - s, empty_region))
})


test_that("`==` knows the empty set from an occupied one", {
  e <- empty_region()
  expect_true(e == empty_region())
  expect_false(e == simplex3())
  expect_false(simplex3() == e)
})


test_that("disjoin passes an empty region through", {
  e <- empty_region()
  expect_identical(disjoin(e), e)
})


test_that("the algebra is closed through an empty intermediate", {
  s <- simplex3()
  nothing <- point_region(theta = c(1, 0, 0)) &
    point_region(theta = c(0, 1, 0))
  expect_true(S7_inherits(nothing, empty_region))
  # An empty result chains straight back into every operator.
  expect_identical(nothing | s, s)
  expect_true(S7_inherits(nothing & s, empty_region))
  expect_true((s - nothing) == s)
})


test_that("a null model refuses an empty region", {
  fam <- multinomial_family(n_trials = 2L, k = 3L)
  expect_error(
    null_model(fam, empty_region()),
    "the null is empty"
  )
})


test_that("an empty region formats as what it is", {
  expect_match(format(empty_region()), "the empty region")
  expect_output(print(empty_region()), "the empty region")
})
