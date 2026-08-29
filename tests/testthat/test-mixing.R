# Properties of R/mixing.R.
#
# A mixing measure is a law over the parameter space and knows nothing of a
# sample space, so nothing here evaluates a density---that is test-mixture.R.

# --- Construction and validation ----------------------------------------------

test_that("finite_dist rejects weights that are not a probability vector", {
  comp <- cbind(c(0.5, 0.5), c(0.2, 0.8))
  expect_error(
    finite_dist(components = comp, weights = c(0.5, 0.6)),
    "sum to 1"
  )
  expect_error(
    finite_dist(components = comp, weights = c(1.2, -0.2)),
    "non-negative"
  )
  expect_error(
    finite_dist(components = comp, weights = 1),
    "one entry per column"
  )
})

test_that("finite_dist rejects components that are not a matrix", {
  expect_error(
    finite_dist(components = c(0.5, 0.5), weights = 1),
    "must be a matrix"
  )
})

test_that("finite_dist accepts a zero weight on a live atom", {
  # Frank-Wolfe drives weights to zero without removing atoms, so this state is
  # reachable mid-fit and must not error.
  m <- finite_dist(
    components = cbind(c(0.5, 0.5), c(0.1, 0.9)),
    weights = c(1, 0)
  )
  expect_equal(n_atoms(m), 2L)
})

# --- Accessors ----------------------------------------------------------------

test_that("atoms always returns a matrix, including for a point mass", {
  # Callers index atoms by column unconditionally; a dropped dimension here
  # would surface much later as a confusing subscript error.
  p <- dirac(theta = c(0.25, 0.75))
  expect_true(is.matrix(atoms(p)))
  expect_equal(atoms(p)[, 1L], c(0.25, 0.75))
})

test_that("weights dispatches on mixing measures without breaking stats", {
  f <- finite_dist(
    components = cbind(c(0.5, 0.5), c(0.1, 0.9)),
    weights = c(0.3, 0.7)
  )
  expect_identical(weights(f), stats::weights(f))
  expect_equal(weights(dirac(theta = c(0.5, 0.5))), 1)
  # The generic still works on the objects it was written for.
  expect_null(weights(stats::lm(mpg ~ wt, data = datasets::mtcars)))
})

# --- Pruning ------------------------------------------------------------------

test_that("prune drops small atoms and renormalises", {
  m <- finite_dist(
    components = cbind(c(0.5, 0.5), c(0.1, 0.9), c(0.2, 0.8)),
    weights = c(0.6, 1e-9, 0.4 - 1e-9)
  )
  p <- prune(m, threshold = 1e-6)
  expect_equal(n_atoms(p), 2L)
  expect_equal(sum(weights(p)), 1)
  expect_equal(p@components[, 1], c(0.5, 0.5))
})

test_that("prune keeps everything when nothing is below the threshold", {
  m <- finite_dist(
    components = cbind(c(0.5, 0.5), c(0.1, 0.9)),
    weights = c(0.4, 0.6)
  )
  expect_equal(prune(m, threshold = 1e-6), m)
})

test_that("prune errors rather than returning an empty mixing measure", {
  m <- finite_dist(
    components = cbind(c(0.5, 0.5), c(0.1, 0.9)),
    weights = c(0.5, 0.5)
  )
  expect_error(prune(m, threshold = 0.9), "no atom has weight above")
})

test_that("prune never merges coincident atoms", {
  # Pruning is by weight alone. Two atoms at the same parameter both survive,
  # which matters because Frank-Wolfe can add an atom where one already sits.
  theta <- c(0.5, 0.5)
  m <- finite_dist(
    components = cbind(theta, theta),
    weights = c(0.5, 0.5)
  )
  expect_equal(n_atoms(prune(m, threshold = 1e-6)), 2L)
})

# --- reference_point ----------------------------------------------------------

test_that("reference_point lands inside the support, mode or not", {
  simplex <- simplex_region(vertices = diag(3))

  # Concentrations above 1: the interior mode.
  expect_equal(reference_point(dirichlet(alpha = c(4, 3, 2))), c(3, 2, 1) / 6)
  # Below 1 the mode is on the boundary or undefined, so it falls back to
  # the mean.
  expect_equal(
    reference_point(dirichlet(alpha = c(0.5, 0.5, 0.5))),
    c(0.5, 0.5, 0.5) / 1.5
  )

  for (w in list(
    dirichlet(alpha = c(4, 3, 2)),
    dirichlet(alpha = c(1, 1, 1)),
    dirichlet(alpha = c(1, 5, 1))
  )) {
    expect_true(contains(simplex, reference_point(w)))
  }

  # A finite_dist returns one of its atoms.
  w <- finite_dist(
    components = cbind(c(0.6, 0.4), c(0.2, 0.8)),
    weights = c(0.3, 0.7)
  )
  expect_true(any(apply(atoms(w), 2L, identical, reference_point(w))))
  expect_equal(reference_point(dirac(theta = c(0.5, 0.5))), c(0.5, 0.5))
})

test_that("reference_point works for a family too, so ripr_init always has one", {
  # When the alternative is not a mixture there is no measure to take a point
  # from, and the family's own space answers instead.
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  p <- reference_point(fam)
  expect_true(contains(fam@parameter_space, p))
  expect_equal(p, rep(1 / 3, 3))

  g <- gaussian_family(dim = 2L)
  expect_true(contains(g@parameter_space, reference_point(g)))
})
