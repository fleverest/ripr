# Properties of R/mixture.R.
#
# A mixture is the pushforward of a mixing measure through a family, so these
# check the integral and its sampler agree with each other and with the family
# they are built from.

# --- Induced density ----------------------------------------------------------

test_that("a point mixing induces the family at its own parameter", {
  fam <- multinomial_family(n_trials = 10, k = 2)
  theta <- c(0.75, 0.25)
  p <- mixture(fam, dirac(theta = theta))
  x <- enumerate_space(fam@sample_space)
  expect_equal(log_density(p, x), kernel_loglik(fam, theta, x))
})

test_that("an induced mixture is itself a probability distribution", {
  fam <- multinomial_family(n_trials = 9, k = 3)
  comp <- cbind(c(0.6, 0.3, 0.1), c(0.2, 0.2, 0.6), c(1 / 3, 1 / 3, 1 / 3))
  p <- mixture(
    fam,
    finite_dist(components = comp, weights = c(0.2, 0.5, 0.3))
  )
  expect_equal(
    sum(exp(log_density(p, enumerate_space(fam@sample_space)))),
    1,
    tolerance = rounding_tol(1)
  )
})

test_that("a finite mixture is the weighted sum of its components", {
  fam <- multinomial_family(n_trials = 6, k = 3)
  w <- c(0.4, 0.6)
  comp <- cbind(c(0.5, 0.3, 0.2), c(0.1, 0.1, 0.8))
  x <- enumerate_space(fam@sample_space)
  p <- mixture(fam, finite_dist(components = comp, weights = w))

  manual <- log(
    w[1] *
      exp(kernel_loglik(fam, comp[, 1], x)) +
      w[2] * exp(kernel_loglik(fam, comp[, 2], x))
  )
  expect_equal(log_density(p, x), manual, tolerance = rounding_tol(1))
})

test_that("a degenerate finite mixing agrees with the point mixing", {
  fam <- multinomial_family(n_trials = 5, k = 3)
  theta <- c(0.5, 0.3, 0.2)
  x <- enumerate_space(fam@sample_space)
  a <- mixture(fam, dirac(theta = theta))
  b <- mixture(
    fam,
    finite_dist(components = matrix(theta, ncol = 1L), weights = 1)
  )
  expect_equal(log_density(a, x), log_density(b, x))
})

test_that("a zero-weight atom contributes nothing", {
  fam <- multinomial_family(n_trials = 5, k = 3)
  x <- enumerate_space(fam@sample_space)
  live <- c(0.5, 0.3, 0.2)
  with_dead <- mixture(
    fam,
    finite_dist(
      components = cbind(live, c(0.1, 0.1, 0.8)),
      weights = c(1, 0)
    )
  )
  expect_equal(log_density(with_dead, x), kernel_loglik(fam, live, x))
})

test_that("a mixture carrying a boundary atom stays finite where it should", {
  # P_W must be strictly positive wherever Q is, or the ratio is undefined; a
  # boundary atom in the mixture is the usual way that gets violated.
  fam <- multinomial_family(n_trials = 5, k = 3)
  p <- mixture(
    fam,
    finite_dist(
      components = cbind(c(0.5, 0.5, 0), c(1 / 3, 1 / 3, 1 / 3)),
      weights = c(0.5, 0.5)
    )
  )
  ld <- log_density(p, enumerate_space(fam@sample_space))
  expect_true(all(is.finite(ld)))
  expect_equal(sum(exp(ld)), 1, tolerance = rounding_tol(1))
})

# --- Sampling -----------------------------------------------------------------

test_that("draw returns the right shape and respects the trial total", {
  fam <- multinomial_family(n_trials = 7, k = 3)
  p <- mixture(
    fam,
    finite_dist(
      components = cbind(c(0.6, 0.3, 0.1), c(0.1, 0.1, 0.8)),
      weights = c(0.5, 0.5)
    )
  )
  set.seed(2)
  d <- draw(p, 200L)
  expect_equal(dim(d), c(200L, 3L))
  expect_true(all(rowSums(d) == 7))
})

test_that("draw handles a mixture where one component is never selected", {
  # The per-component loop must not assume every component draws at least once.
  fam <- multinomial_family(n_trials = 4, k = 2)
  p <- mixture(
    fam,
    finite_dist(
      components = cbind(c(0.5, 0.5), c(0.9, 0.1)),
      weights = c(1, 0)
    )
  )
  set.seed(5)
  d <- draw(p, 30L)
  expect_equal(dim(d), c(30L, 2L))
  expect_false(anyNA(d))
})

test_that("draws from a mixture are consistent with its density", {
  # A slow-but-honest check that mixture_draw samples the law
  # mixture_log_density describes: empirical frequencies converge to the pmf.
  skip_on_cran()
  fam <- multinomial_family(n_trials = 4, k = 2)
  p <- mixture(
    fam,
    finite_dist(
      components = cbind(c(0.75, 0.25), c(0.25, 0.75)),
      weights = c(0.5, 0.5)
    )
  )
  set.seed(42)
  d <- draw(p, 2e5)
  x <- enumerate_space(fam@sample_space)
  empirical <- vapply(
    seq_len(nrow(x)),
    \(i) mean(d[, 1] == x[i, 1]),
    numeric(1)
  )
  expect_equal(empirical, exp(log_density(p, x)), tolerance = 0.01)
})

# --- Family recovery ----------------------------------------------------------

test_that("mixture rejects components of the wrong type", {
  fam <- multinomial_family(n_trials = 3, k = 2)
  expect_error(mixture(fam, fam))
  expect_error(mixture("not a family", dirac(theta = c(0.5, 0.5))))
})

# --- Callable families --------------------------------------------------------

test_that("calling a family is the map theta -> p_theta", {
  # A family IS the kernel, so applying it to a parameter should be spelled the
  # way applying a function is.
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  theta <- c(0.5, 0.3, 0.2)
  expect_identical(fam(theta), mixture(fam, theta))
  expect_true(S7_inherits(fam(theta), distribution))
})

test_that("a kernel extends from points to measures, so fam(W) is the same map", {
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  w <- finite_dist(
    components = cbind(c(0.6, 0.2, 0.2), c(0.2, 0.6, 0.2)),
    weights = c(0.3, 0.7)
  )
  expect_identical(fam(w), mixture(fam, w))

  # The degenerate case agrees with the general one rather than being special.
  theta <- c(0.5, 0.3, 0.2)
  x <- enumerate_space(fam@sample_space)
  expect_equal(
    log_density(fam(theta), x),
    log_density(fam(dirac(theta = theta)), x)
  )
})

test_that("the callable closure captures nothing", {
  # Built once at namespace level rather than per constructor. A closure made
  # inside a constructor would close over that frame and carry a second copy of
  # the family's own properties into every serialised family, for a reference
  # nothing reads. Checked by where the closure lives rather than by measuring
  # a size, since S7 class metadata dominates either number.
  fam <- multinomial_family(n_trials = 60L, k = 4L)
  expect_identical(environment(fam), asNamespace("ripr"))
  expect_identical(environment(fam), environment(gaussian_family(dim = 2L)))
})

test_that("a family survives a serialisation round trip and stays callable", {
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  back <- unserialize(serialize(fam, NULL))
  expect_identical(back(c(0.5, 0.3, 0.2)), fam(c(0.5, 0.3, 0.2)))
})

test_that("a family is callable through the usual indirections", {
  # `sys.function()` has to recover the family with its S7 attributes intact
  # however the call was written.
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  theta <- c(0.5, 0.3, 0.2)
  expect_identical(do.call(fam, list(theta)), fam(theta))
  expect_identical(lapply(list(theta), fam)[[1L]], fam(theta))
  expect_identical((function(f, t) f(t))(fam, theta), fam(theta))
})

# --- Printing -----------------------------------------------------------------

test_that("a family prints its two spaces, not its closure", {
  # The parent is `class_function`, so the defaults reach `deparse()` and dump
  # the shared closure plus every attribute.
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  expect_equal(format(fam), "multinomial_family: simplex_region -> count_space")
  out <- capture.output(print(fam))
  expect_match(out[[1L]], "multinomial_family")
  expect_true(any(grepl("simplex_region, dimension 3", out, fixed = TRUE)))
  expect_true(any(grepl("count_space, dimension 3", out, fixed = TRUE)))
  expect_false(any(grepl("sys.function", out, fixed = TRUE)))
})

test_that("a distribution's print distinguishes a point mass from a mixture", {
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  expect_match(
    format(fam(c(0.5, 0.3, 0.2))),
    "at theta = (0.5, 0.3, 0.2)",
    fixed = TRUE
  )
  expect_match(
    format(fam(finite_dist(
      components = cbind(c(0.6, 0.2, 0.2), c(0.2, 0.6, 0.2)),
      weights = c(0.5, 0.5)
    ))),
    "mixed over 2 atoms",
    fixed = TRUE
  )
})


test_that("draw samples every distribution over a parameter space", {
  set.seed(1)
  p <- draw(dirac(theta = c(0.5, 0.3, 0.2)), 4L)
  expect_equal(dim(p), c(4L, 3L))
  expect_true(all(apply(p, 1L, identical, c(0.5, 0.3, 0.2))))

  # Weights govern how often each atom is drawn, and repeats stay in place.
  w <- finite_dist(
    components = cbind(c(1, 0, 0), c(0, 1, 0)),
    weights = c(0.25, 0.75)
  )
  d <- draw(w, 4000L)
  expect_equal(dim(d), c(4000L, 3L))
  expect_equal(mean(d[, 2L]), 0.75, tolerance = 0.03)
})

test_that("one mixture_draw method serves every mixing measure", {
  # The collapse must not have changed any measure's behaviour: a point mass
  # still gives its own kernel, a finite mixture still respects its weights.
  fam <- multinomial_family(n_trials = 30L, k = 3L)

  set.seed(1)
  x <- draw(fam(dirac(theta = c(1, 0, 0))), 5L)
  expect_equal(x, matrix(rep(c(30L, 0L, 0L), each = 5L), nrow = 5L))

  set.seed(1)
  degenerate <- finite_dist(
    components = cbind(c(1, 0, 0), c(0, 0, 1)),
    weights = c(0.5, 0.5)
  )
  y <- draw(fam(degenerate), 2000L)
  # Every draw is one atom or the other, in roughly equal numbers.
  expect_true(all(y[, 2L] == 0L))
  expect_equal(mean(y[, 1L] == 30L), 0.5, tolerance = 0.05)
})

test_that("a closed-form pairing still overrides the general method", {
  # A Gaussian prior through a Gaussian kernel is N(m, Sigma + V) directly and
  # never samples a parameter, so it must beat the generic two-step method.
  fam <- gaussian_family(dim = 1L, sigma = matrix(1))
  prior <- gaussian_dist(prior_mean = 0, prior_cov = matrix(3))
  set.seed(1)
  x <- draw(fam(prior), 2e5)
  expect_equal(var(as.vector(x)), 4, tolerance = 0.05)
})


# --- What a mixture refuses ---------------------------------------------------

test_that("a mixture refuses a measure that lives outside the parameters", {
  fam <- multinomial_family(n_trials = 6L, k = 3L)

  expect_error(fam(c(0.5, 0.3, 0.9)), "outside the family's parameter space")
  expect_error(fam(dirac(theta = c(0.5, 0.5))), "over 2 dimensions")
  expect_error(
    fam(finite_dist(
      components = cbind(c(0.5, 0.3, 0.2), c(0.5, 0.3, 0.9)),
      weights = c(0.5, 0.5)
    )),
    "outside the family's parameter space"
  )

  # A Gaussian prior is over all of R^3, which is not inside the simplex. This
  # pairing was silently accepted before, and fails at the density instead.
  expect_error(
    fam(gaussian_dist(prior_mean = rep(0, 3), prior_cov = diag(3))),
    "outside the family's parameter space"
  )
})

test_that("a mixture accepts every measure that does live there", {
  fam <- multinomial_family(n_trials = 6L, k = 3L)
  simplex <- fam@parameter_space

  for (w in list(
    dirac(theta = c(0.5, 0.3, 0.2)),
    finite_dist(
      components = cbind(c(.6, .2, .2), c(.2, .6, .2)),
      weights = c(.5, .5)
    ),
    dirichlet(alpha = c(2, 2, 2))
  )) {
    expect_true(S7::S7_inherits(fam(w), mixture))
  }

  # And a truncated Dirichlet, whose declared space is a strict subset -- the
  # case an equality check would have rejected.
  medial <- simplex_region(
    vertices = cbind(c(.5, .5, 0), c(.5, 0, .5), c(0, .5, .5))
  )
  w <- truncated_dirichlet(alpha = c(2, 2, 2), region = medial)
  expect_false(setequal(w@sample_space, simplex))
  expect_true(S7::S7_inherits(fam(w), mixture))
})

test_that("the support check is exact for atoms, containment for the rest", {
  simplex <- simplex_region(vertices = diag(3))
  inside <- finite_dist(components = cbind(c(.6, .2, .2)), weights = 1)
  outside <- finite_dist(components = cbind(c(.6, .2, .9)), weights = 1)

  expect_true(ripr:::supported_in(inside, simplex))
  expect_false(ripr:::supported_in(outside, simplex))

  expect_true(S7::S7_inherits(outside@sample_space, real_region))
  expect_true(ripr:::supported_in(dirichlet(alpha = c(2, 2, 2)), simplex))
})
