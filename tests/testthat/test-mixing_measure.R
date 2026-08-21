# Properties of R/mixing_measure.R.
#
# A mixing measure is a law over the parameter space and knows nothing of a
# sample space, so nothing here evaluates a density---that is test-mixture.R.

# --- Construction and validation ----------------------------------------------

test_that("finite_mixing rejects weights that are not a probability vector", {
  comp <- cbind(c(0.5, 0.5), c(0.2, 0.8))
  expect_error(
    finite_mixing(components = comp, weights = c(0.5, 0.6)),
    "sum to 1"
  )
  expect_error(
    finite_mixing(components = comp, weights = c(1.2, -0.2)),
    "non-negative"
  )
  expect_error(
    finite_mixing(components = comp, weights = 1),
    "one entry per column"
  )
})

test_that("finite_mixing rejects components that are not a matrix", {
  expect_error(
    finite_mixing(components = c(0.5, 0.5), weights = 1),
    "must be a matrix"
  )
})

test_that("finite_mixing accepts a zero weight on a live atom", {
  # Frank-Wolfe drives weights to zero without removing atoms, so this state is
  # reachable mid-fit and must not error.
  m <- finite_mixing(
    components = cbind(c(0.5, 0.5), c(0.1, 0.9)),
    weights = c(1, 0)
  )
  expect_equal(n_atoms(m), 2L)
})

# --- Accessors ----------------------------------------------------------------

test_that("atoms, weights and n_atoms agree across mixing-measure kinds", {
  p <- point_mixing(theta_star = c(0.5, 0.5))
  f <- finite_mixing(
    components = cbind(c(0.5, 0.5), c(0.1, 0.9)),
    weights = c(0.3, 0.7)
  )
  expect_equal(n_atoms(p), 1L)
  expect_equal(dim(atoms(p)), c(2L, 1L))
  expect_equal(weights(p), 1)
  expect_equal(n_atoms(f), 2L)
  expect_equal(atoms(f), f@components)
  expect_equal(weights(f), f@weights)
})

test_that("atoms always returns a matrix, including for a point mass", {
  # Callers index atoms by column unconditionally; a dropped dimension here
  # would surface much later as a confusing subscript error.
  p <- point_mixing(theta_star = c(0.25, 0.75))
  expect_true(is.matrix(atoms(p)))
  expect_equal(atoms(p)[, 1L], c(0.25, 0.75))
})

test_that("S7 namespaces the class attribute", {
  # For the weights methods: S3 dispatch sees "ripr::finite_mixing",
  # not "finite_mixing", which is why the methods are registered via S7's
  # method<- rather than defined as weights.finite_mixing.
  f <- finite_mixing(
    components = cbind(c(0.5, 0.5), c(0.1, 0.9)),
    weights = c(0.3, 0.7)
  )
  expect_true("ripr::finite_mixing" %in% class(f))
  expect_true("ripr::mixing_measure" %in% class(f))
})

test_that("weights dispatches on mixing measures without breaking stats", {
  f <- finite_mixing(
    components = cbind(c(0.5, 0.5), c(0.1, 0.9)),
    weights = c(0.3, 0.7)
  )
  expect_identical(weights(f), stats::weights(f))
  expect_equal(weights(point_mixing(theta_star = c(0.5, 0.5))), 1)
  # The generic still works on the objects it was written for.
  expect_null(weights(stats::lm(mpg ~ wt, data = datasets::mtcars)))
})


test_that("the class hierarchy is as declared", {
  # A mixing measure is a law over the parameter space and a `distribution` is
  # a law over the sample space. They are peers, not branches of a shared root:
  # a mixing measure has atoms and no density, an outcome law the reverse, so
  # there is no interface for a common supertype to carry.
  f <- finite_mixing(components = cbind(c(0.5, 0.5)), weights = 1)
  expect_true(S7_inherits(f, mixing_measure))
  expect_false(S7_inherits(f, distribution))
  expect_true(S7_inherits(point_mixing(theta_star = c(1)), mixing_measure))
})

# --- Pruning ------------------------------------------------------------------

test_that("prune drops small atoms and renormalises", {
  m <- finite_mixing(
    components = cbind(c(0.5, 0.5), c(0.1, 0.9), c(0.2, 0.8)),
    weights = c(0.6, 1e-9, 0.4 - 1e-9)
  )
  p <- prune(m, threshold = 1e-6)
  expect_equal(n_atoms(p), 2L)
  expect_equal(sum(weights(p)), 1)
  expect_equal(p@components[, 1], c(0.5, 0.5))
})

test_that("prune keeps everything when nothing is below the threshold", {
  m <- finite_mixing(
    components = cbind(c(0.5, 0.5), c(0.1, 0.9)),
    weights = c(0.4, 0.6)
  )
  expect_equal(prune(m, threshold = 1e-6), m)
})

test_that("prune errors rather than returning an empty mixing measure", {
  m <- finite_mixing(
    components = cbind(c(0.5, 0.5), c(0.1, 0.9)),
    weights = c(0.5, 0.5)
  )
  expect_error(prune(m, threshold = 0.9), "no atom has weight above")
})

test_that("prune never merges coincident atoms", {
  # Pruning is by weight alone. Two atoms at the same parameter both survive,
  # which matters because Frank-Wolfe can add an atom where one already sits.
  theta <- c(0.5, 0.5)
  m <- finite_mixing(
    components = cbind(theta, theta),
    weights = c(0.5, 0.5)
  )
  expect_equal(n_atoms(prune(m, threshold = 1e-6)), 2L)
})
