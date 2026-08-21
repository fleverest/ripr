# Properties of R/mixture.R.
#
# A mixture is the pushforward of a mixing measure through a family, so these
# check the integral and its sampler agree with each other and with the family
# they are built from.

# --- Induced density ----------------------------------------------------------

test_that("a point mixing induces the family at its own parameter", {
  fam <- multinomial_family(n_trials = 10, k = 2)
  theta <- c(0.75, 0.25)
  p <- induced_distribution(fam, point_mixing(theta_star = theta))
  x <- enumerate_space(fam@sample_space)
  expect_equal(log_density(p, x), kernel_loglik(fam, theta, x))
})

test_that("an induced mixture is itself a probability distribution", {
  fam <- multinomial_family(n_trials = 9, k = 3)
  comp <- cbind(c(0.6, 0.3, 0.1), c(0.2, 0.2, 0.6), c(1 / 3, 1 / 3, 1 / 3))
  p <- induced_distribution(fam, finite_mixing(components = comp, weights = c(0.2, 0.5, 0.3)))
  expect_equal(
    sum(exp(log_density(p, enumerate_space(fam@sample_space)))),
    1,
    tolerance = 1e-12
  )
})

test_that("a finite mixture is the weighted sum of its components", {
  fam <- multinomial_family(n_trials = 6, k = 3)
  w <- c(0.4, 0.6)
  comp <- cbind(c(0.5, 0.3, 0.2), c(0.1, 0.1, 0.8))
  x <- enumerate_space(fam@sample_space)
  p <- induced_distribution(fam, finite_mixing(components = comp, weights = w))

  manual <- log(
    w[1] *
      exp(kernel_loglik(fam, comp[, 1], x)) +
      w[2] * exp(kernel_loglik(fam, comp[, 2], x))
  )
  expect_equal(log_density(p, x), manual, tolerance = 1e-12)
})

test_that("a degenerate finite mixing agrees with the point mixing", {
  fam <- multinomial_family(n_trials = 5, k = 3)
  theta <- c(0.5, 0.3, 0.2)
  x <- enumerate_space(fam@sample_space)
  a <- induced_distribution(fam, point_mixing(theta_star = theta))
  b <- induced_distribution(fam, finite_mixing(components = matrix(theta, ncol = 1L), weights = 1))
  expect_equal(log_density(a, x), log_density(b, x))
})

test_that("a zero-weight atom contributes nothing", {
  fam <- multinomial_family(n_trials = 5, k = 3)
  x <- enumerate_space(fam@sample_space)
  live <- c(0.5, 0.3, 0.2)
  with_dead <- induced_distribution(fam, finite_mixing(
      components = cbind(live, c(0.1, 0.1, 0.8)),
      weights = c(1, 0)
    ))
  expect_equal(log_density(with_dead, x), kernel_loglik(fam, live, x))
})

test_that("a mixture carrying a boundary atom stays finite where it should", {
  # P_W must be strictly positive wherever Q is, or the ratio is undefined; a
  # boundary atom in the mixture is the usual way that gets violated.
  fam <- multinomial_family(n_trials = 5, k = 3)
  p <- induced_distribution(fam, finite_mixing(
      components = cbind(c(0.5, 0.5, 0), c(1 / 3, 1 / 3, 1 / 3)),
      weights = c(0.5, 0.5)
    ))
  ld <- log_density(p, enumerate_space(fam@sample_space))
  expect_true(all(is.finite(ld)))
  expect_equal(sum(exp(ld)), 1, tolerance = 1e-12)
})

test_that("log_density accepts a bare vector as one outcome", {
  fam <- multinomial_family(n_trials = 10, k = 2)
  p <- induced_distribution(fam, point_mixing(theta_star = c(0.75, 0.25)))
  expect_length(log_density(p, c(8, 2)), 1L)
})

# --- Sampling -----------------------------------------------------------------

test_that("draw returns the right shape and respects the trial total", {
  fam <- multinomial_family(n_trials = 7, k = 3)
  p <- induced_distribution(fam, finite_mixing(
      components = cbind(c(0.6, 0.3, 0.1), c(0.1, 0.1, 0.8)),
      weights = c(0.5, 0.5)
    ))
  set.seed(2)
  d <- draw(p, 200L)
  expect_equal(dim(d), c(200L, 3L))
  expect_true(all(rowSums(d) == 7))
})

test_that("draw handles a mixture where one component is never selected", {
  # The per-component loop must not assume every component draws at least once.
  fam <- multinomial_family(n_trials = 4, k = 2)
  p <- induced_distribution(fam, finite_mixing(
      components = cbind(c(0.5, 0.5), c(0.9, 0.1)),
      weights = c(1, 0)
    ))
  set.seed(5)
  d <- draw(p, 30L)
  expect_equal(dim(d), c(30L, 2L))
  expect_false(anyNA(d))
})

test_that("draws from a mixture are consistent with its density", {
  # A slow-but-honest check that induced_draw samples the law
  # induced_log_density describes: empirical frequencies converge to the pmf.
  skip_on_cran()
  fam <- multinomial_family(n_trials = 4, k = 2)
  p <- induced_distribution(fam, finite_mixing(
      components = cbind(c(0.75, 0.25), c(0.25, 0.75)),
      weights = c(0.5, 0.5)
    ))
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

test_that("family reports the family for a mixture and NULL otherwise", {
  fam <- multinomial_family(n_trials = 3, k = 2)
  p <- induced_distribution(fam, point_mixing(theta_star = c(0.5, 0.5)))
  expect_identical(family(p), fam)

  direct <- new_class("direct_law", parent = distribution)
  expect_null(family(direct(sample_space = fam@sample_space)))
})

test_that("mixture rejects components of the wrong type", {
  fam <- multinomial_family(n_trials = 3, k = 2)
  expect_error(induced_distribution(fam, fam))
  expect_error(induced_distribution("not a family", point_mixing(theta_star = c(0.5, 0.5))))
})

# --- Callable families --------------------------------------------------------

test_that("calling a family is the map theta -> p_theta", {
  # A family IS the kernel, so applying it to a parameter should be spelled the
  # way applying a function is.
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  theta <- c(0.5, 0.3, 0.2)
  expect_identical(fam(theta), induced_distribution(fam, theta))
  expect_true(S7_inherits(fam(theta), distribution))
  expect_equal(
    log_density(fam(theta), c(2L, 1L, 1L)),
    kernel_loglik(fam, theta, c(2L, 1L, 1L))
  )
})

test_that("a kernel extends from points to measures, so fam(W) is the same map", {
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  w <- finite_mixing(
    components = cbind(c(0.6, 0.2, 0.2), c(0.2, 0.6, 0.2)),
    weights = c(0.3, 0.7)
  )
  expect_identical(fam(w), induced_distribution(fam, w))

  # The degenerate case agrees with the general one rather than being special.
  theta <- c(0.5, 0.3, 0.2)
  x <- enumerate_space(fam@sample_space)
  expect_equal(
    log_density(fam(theta), x),
    log_density(fam(point_mixing(theta_star = theta)), x)
  )
})

test_that("being callable costs a family none of its other behaviour", {
  fam <- gaussian_family(dim = 2L, sigma = diag(c(1, 4)))
  expect_true(S7_inherits(fam, parametric_family))
  expect_equal(fam@sigma, diag(c(1, 4)))
  expect_equal(space_dim(fam@parameter_space), 2L)
  expect_true(is.function(compile_loglik(fam, rbind(c(0, 0)))))
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
  expect_false(any(c("n_trials", "k", "space") %in% ls(environment(fam))))
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
  expect_match(format(fam(c(0.5, 0.3, 0.2))), "at theta = (0.5, 0.3, 0.2)", fixed = TRUE)
  expect_match(
    format(fam(finite_mixing(
      components = cbind(c(0.6, 0.2, 0.2), c(0.2, 0.6, 0.2)),
      weights = c(0.5, 0.5)
    ))),
    "mixed over 2 atoms",
    fixed = TRUE
  )
})

# --- The sample space a distribution is a law over ----------------------------

test_that("a distribution carries the sample space it is a law over", {
  # So a bespoke Q with no family behind it still tells a random variable what
  # it may be evaluated on.
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  expect_identical(fam(c(0.5, 0.3, 0.2))@sample_space, fam@sample_space)
  expect_identical(likelihood(fam(c(0.5, 0.3, 0.2)))@sample_space, fam@sample_space)
})
