# The certification contract: `gap_used` is a conservative one-sided bound built
# from a single fresh sample, and `estimate = TRUE` only ever *adds* a
# diagnostic. Nothing about the diagnostic may reach the rescaling constant.

# README Example 2: X ~ N(mu, I_2), null = {mu_1 <= mu_2} U {mu_1 <= -mu_2},
# alternative mu = (2, 1).
cone_problem <- function(n_draws = 1000L) {
  fam <- gaussian_family(dim = 2)
  Q <- as_marginal(point_mixing(theta_star = c(2, 1)), fam)
  faces <- list(
    halfspace_face(v = c(1, -1), c = 0, face_index = 1),
    halfspace_face(v = c(1, 1), c = 0, face_index = 2)
  )
  eng <- mc_engine(fam, Q, n_draws = n_draws)
  ripr_problem(fam, null_region(faces = faces), Q, engine = eng)
}

# The exact RIPr for that problem, from an independent 2-D Gauss-Hermite solve:
# KL* = 0.2313753560, and the duality gap at this mixture is below 1e-7.
cone_optimum <- function() {
  finite_mixing(
    components = cbind(c(1.548311, 1.548311), c(1.267288, -1.267288)),
    weights = c(0.95099, 0.04901)
  )
}

test_that("`estimate` does not move `gap_used`", {
  prob <- cone_problem()
  proj <- cone_optimum()

  set.seed(404)
  a <- certify(proj, prob, n_draws = 2000L, n_seeds = 20, estimate = FALSE)
  set.seed(404)
  b <- certify(proj, prob, n_draws = 2000L, n_seeds = 20, estimate = TRUE)

  # The rescaling constant, and everything derived from it, must be untouched.
  expect_identical(a$gap_used, b$gap_used)
  expect_identical(a$gap, b$gap)
  expect_identical(a$gap_se, b$gap_se)
  expect_identical(a$kl_lower_bound, b$kl_lower_bound)
  expect_identical(a$growth_rate, b$growth_rate)
  expect_identical(a$oracle_theta, b$oracle_theta)
})

test_that("`gap_est` is NA unless requested, and finite when requested", {
  prob <- cone_problem()
  proj <- cone_optimum()

  set.seed(11)
  off <- certify(proj, prob, n_draws = 1000L, n_seeds = 20)
  expect_true(is.na(off$gap_est))
  expect_identical(off$gap_est, NA_real_)

  set.seed(11)
  on <- certify(proj, prob, n_draws = 1000L, n_seeds = 20, estimate = TRUE)
  expect_false(is.na(on$gap_est))
  expect_true(is.finite(on$gap_est))

  # Exact engines never estimate: no sample to draw.
  fam <- multinomial_family(n_trials = 6, k = 2)
  alt <- point_mixing(theta_star = c(0.75, 0.25))
  face <- polytope_face(
    vertices = cbind(c(0, 1), c(0.5, 0.5)),
    face_index = 1
  )
  eprob <- ripr_problem(fam, null_region(faces = list(face)), as_marginal(alt, fam))
  ecert <- certify(
    finite_mixing(components = matrix(c(0.5, 0.5), ncol = 1), weights = 1),
    eprob,
    estimate = TRUE
  )
  expect_identical(ecert$gap_est, NA_real_)
  expect_equal(ecert$gap_se, 0)
})

test_that("the Monte Carlo certificate errs upward at the known optimum", {
  prob <- cone_problem()
  proj <- cone_optimum()

  set.seed(2024)
  # `ess_min = 0`: this seed happens to land in the noise-driven regime the
  # accuracy section describes (ess ~ 5). The direction assertions below hold
  # in both regimes, which is the point of this test; the collapsed-ESS path
  # has its own test.
  cert <- certify(proj, prob, n_draws = 40000L, n_seeds = 60, ess_min = 0)

  # Inflation is one-sided.
  expect_gte(cert$gap_used, cert$gap)
  expect_gt(cert$gap_se, 0)

  # The true gap at this mixture is < 1e-7; the single-sample maximum plus
  # inflation must land well above it. A direction test, not a precision test.
  expect_gt(cert$gap_used, 1e-5)

  # NOTE: there is deliberately no upper bound asserted here. `gap_used` is
  # `max_theta Ghat(theta)` inflated, and that maximum is not a consistent
  # estimator of `sup_theta G` at fixed `n_draws`: the oracle can find regions
  # of the face that almost no Q-draw reaches, where `Ghat` is a mean of a few
  # enormous ratios. Certifying this known-optimal mixture 25 times at 40k
  # draws gave a median gap of 0.004 but exceeded 0.1 in 20% of runs (max 2.74),
  # and the inflated runs are exactly those where |theta*| ~ 7 rather than ~ 2.
  # The certificate stays conservative (the safe direction) but can be vacuous.
  # Bounding this needs a variance-aware or uniform-concentration estimator,
  # which is not implemented.

  # KL is recovered to Monte Carlo accuracy (exact value 0.2313753560).
  expect_equal(cert$kl, 0.2313753560, tolerance = 0.05)
})

test_that("`estimate = TRUE` draws a second sample", {
  prob <- cone_problem()
  proj <- cone_optimum()

  count_resamples <- function(estimate) {
    n <- 0L
    real <- ripr:::resample_engine
    local_mocked_bindings(
      resample_engine = function(engine, n_draws) {
        n <<- n + 1L
        real(engine, n_draws)
      },
      .package = "ripr"
    )
    set.seed(5)
    certify(proj, prob, n_draws = 500L, n_seeds = 10, estimate = estimate)
    n
  }

  expect_identical(count_resamples(FALSE), 1L)
  expect_identical(count_resamples(TRUE), 2L)
})

test_that("the certificate reports the effective sample size behind the gap", {
  prob <- cone_problem()
  proj <- cone_optimum()

  set.seed(77)
  cert <- certify(proj, prob, n_draws = 5000L, n_seeds = 30)
  expect_true(is.finite(cert$ess))
  expect_gt(cert$ess, 0)
  expect_lte(cert$ess, 5000)
  # A healthy certificate of the known optimum rests on a decent share of the
  # sample -- it is not carried by a handful of draws.
  expect_gt(cert$ess, 50)

  # Exact engines do not sample, so there is nothing to be effective about.
  fam <- multinomial_family(n_trials = 6, k = 2)
  face <- polytope_face(vertices = cbind(c(0, 1), c(0.5, 0.5)), face_index = 1)
  eprob <- ripr_problem(
    fam,
    null_region(faces = list(face)),
    as_marginal(point_mixing(theta_star = c(0.75, 0.25)), fam)
  )
  ecert <- certify(
    finite_mixing(components = matrix(c(0.5, 0.5), ncol = 1), weights = 1),
    eprob
  )
  expect_true(is.na(ecert$ess))
})

test_that("a collapsed effective sample size warns, and ess_min controls it", {
  prob <- cone_problem()
  proj <- cone_optimum()

  # Drive the threshold rather than fish for a rare seed: the warning path is
  # what needs testing, not the frequency of the event.
  set.seed(78)
  expect_warning(
    certify(proj, prob, n_draws = 2000L, n_seeds = 20, ess_min = Inf),
    "effective sample size"
  )
  set.seed(78)
  expect_no_warning(
    certify(proj, prob, n_draws = 2000L, n_seeds = 20, ess_min = 0)
  )
  # The threshold must not touch the numbers it warns about.
  set.seed(78)
  a <- suppressWarnings(
    certify(proj, prob, n_draws = 2000L, n_seeds = 20, ess_min = Inf)
  )
  set.seed(78)
  b <- certify(proj, prob, n_draws = 2000L, n_seeds = 20, ess_min = 0)
  expect_identical(a$gap_used, b$gap_used)
  expect_identical(a$ess, b$ess)
})

test_that("certify() has no `split` argument", {
  expect_false("split" %in% names(formals(certify)))
  expect_identical(formals(certify)$estimate, FALSE)
})
