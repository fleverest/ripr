# Provenance: "is this gap certified?" must be machine-readable, must be the
# AND across faces, and must never be TRUE by accident. The paths that cannot
# be certified have to degrade quietly to exactly what they did before -- a
# heuristic gap that *looks* certified is worse than the documented status quo.

# The Gaussian cone: a stochastic engine, a family with no Bernstein form, and
# unbounded faces. Nothing about it is certifiable, and none of it may drift.
cone_problem_mc <- function(n_draws = 1000L) {
  fam <- gaussian_family(dim = 2)
  Q <- as_marginal(point_mixing(theta_star = c(2, 1)), fam)
  faces <- list(
    halfspace_face(v = c(1, -1), c = 0, face_index = 1),
    halfspace_face(v = c(1, 1), c = 0, face_index = 2)
  )
  ripr_problem(fam, null_region(faces = faces), Q,
               engine = mc_engine(fam, Q, n_draws = n_draws))
}

cone_mixture <- function() {
  finite_mixing(
    components = cbind(c(1.548311, 1.548311), c(1.267288, -1.267288)),
    weights = c(0.95099, 0.04901)
  )
}

# A certifiable multinomial problem on the plurality null.
mn_problem <- function(n = 5L, K = 3L) {
  fam <- multinomial_family(n_trials = n, k = K)
  alt <- point_mixing(theta_star = c(0.6, 0.25, 0.15))
  ripr_problem(fam, null_region(faces = plurality_faces(K)),
               as_marginal(alt, fam))
}

mn_mixture <- function() {
  finite_mixing(
    components = cbind(c(0.4, 0.4, 0.2), c(0.3, 0.35, 0.35)),
    weights = c(0.6, 0.4)
  )
}

test_that("a Monte Carlo Gaussian problem reports an uncertified gap", {
  prob <- cone_problem_mc()
  proj <- cone_mixture()

  set.seed(404)
  cert <- certify(proj, prob, n_draws = 2000L, n_seeds = 20)

  expect_false(cert$gap_certified)
  expect_identical(cert$gap_method, "multistart_bfgs")
  expect_identical(cert$bnb_iterations, NA_integer_)
})

test_that("`bnb` is inert on the paths it cannot certify", {
  prob <- cone_problem_mc()
  proj <- cone_mixture()

  # Same seed, with and without the new argument: every number must match. The
  # new bound is opt-in, and opting in to something unavailable must change
  # nothing at all.
  set.seed(404)
  before <- certify(proj, prob, n_draws = 2000L, n_seeds = 20)
  set.seed(404)
  after <- certify(proj, prob, n_draws = 2000L, n_seeds = 20, bnb = bnb_control())

  for (field in c("kl", "gap", "gap_se", "gap_used", "kl_lower_bound",
                  "growth_rate", "oracle_theta", "oracle_face", "ess")) {
    expect_identical(before[[field]], after[[field]], info = field)
  }
  expect_false(after$gap_certified)
})

test_that("`require_certified` errors and names what blocked it", {
  prob <- cone_problem_mc()
  proj <- cone_mixture()

  set.seed(404)
  err <- expect_error(
    certify(proj, prob, n_draws = 500L, n_seeds = 10, require_certified = TRUE)
  )
  msg <- conditionMessage(err)
  expect_match(msg, "cannot be certified")
  # The blocker must be named, not merely asserted.
  expect_match(msg, "mc_engine|gaussian_family")

  # And on an exact multinomial problem it is satisfiable, so it must not error.
  mp <- mn_problem()
  expect_no_error(
    cert <- certify(mn_mixture(), mp, require_certified = TRUE)
  )
  expect_true(cert$gap_certified)
})

test_that("the capability lives on the family, not the engine", {
  # Being exact is necessary but not sufficient: the Bernstein identity is a
  # property of the *multinomial* pmf, so an exact engine over some other
  # family must not inherit certifiability from its exactness.
  expect_null(bernstein_form(gaussian_family(dim = 2)))

  bf <- bernstein_form(multinomial_family(n_trials = 5L, k = 3L))
  expect_identical(bf$n, 5L)
  expect_identical(bf$K, 3L)
  # The tally must be the family's own support, row for row: that alignment is
  # what makes `exp(engine@log_q_mass - log_Pw)` the coefficient vector.
  expect_identical(bf$tally, support(multinomial_family(n_trials = 5L, k = 3L)))

  # An exact engine paired with a family that has no Bernstein form is blocked
  # by the family, and says so.
  mn <- multinomial_family(n_trials = 5L, k = 3L)
  eng <- exact_engine(mn, as_marginal(point_mixing(c(0.6, 0.25, 0.15)), mn))
  face <- plurality_faces(3L)[[1L]]
  expect_true(deterministic(eng))
  reasons <- ripr:::certifiable_reason(face, eng, gaussian_family(dim = 2))
  expect_match(paste(reasons, collapse = " "), "gaussian_family")
  expect_match(paste(reasons, collapse = " "), "Bernstein form")
  expect_false(certifiable(face, eng, gaussian_family(dim = 2)))
})

test_that("every applicable blocker is reported, not just the first", {
  # A stochastic engine over a Gaussian family with an unbounded face fails on
  # all three counts. Reporting only the first would send a user round the loop
  # three times.
  prob <- cone_problem_mc(n_draws = 200L)
  reasons <- ripr:::certifiable_reason(
    prob$null[[1L]], prob$engine, prob$family
  )
  expect_length(reasons, 3L)
  joined <- paste(reasons, collapse = " ")
  expect_match(joined, "mc_engine")
  expect_match(joined, "gaussian_family")
  expect_match(joined, "halfspace_face")
})

test_that("one uncertifiable face makes the whole gap uncertified", {
  K <- 3L
  fam <- multinomial_family(n_trials = 5L, k = K)
  alt <- point_mixing(theta_star = c(0.6, 0.25, 0.15))

  simplex_face <- plurality_faces(K)[[1L]]
  # A four-vertex polygon inside the simplex: a perfectly good polytope_face,
  # but not a simplex, so the Bernstein reparametrisation does not apply.
  polygon_face <- polytope_face(
    vertices = cbind(c(0, 1, 0), c(0, 0, 1), c(0.5, 0.5, 0), c(0.5, 0, 0.5)),
    face_index = 2
  )

  eng <- exact_engine(fam, as_marginal(alt, fam))
  ctl <- bnb_control()
  expect_true(certifiable(simplex_face, eng, fam, ctl))
  expect_false(certifiable(polygon_face, eng, fam, ctl))

  mixed <- ripr_problem(
    fam,
    null_region(faces = list(simplex_face, polygon_face)),
    as_marginal(alt, fam)
  )
  set.seed(8)
  cert <- certify(mn_mixture(), mixed, bnb = ctl, n_seeds = 50L)

  # Certification is over a union: partly bounded is not bounded.
  expect_false(cert$gap_certified)
  expect_identical(cert$gap_method, "multistart_bfgs")
  expect_identical(cert$bnb_iterations, NA_integer_)

  # And it must be identical to the untouched path, not merely close.
  set.seed(8)
  plain <- certify(mn_mixture(), mixed, n_seeds = 50L)
  expect_identical(cert$gap_used, plain$gap_used)
})

test_that("exceeding `max_coef` falls back rather than erroring", {
  prob <- mn_problem(n = 5L, K = 3L)
  # 21 coefficients at n = 5, K = 3; a budget of 5 cannot cover it.
  tiny <- bnb_control(max_coef = 5)

  expect_false(certifiable(prob$null[[1L]], prob$engine, prob$family, tiny))

  set.seed(3)
  cert <- expect_no_error(certify(mn_mixture(), prob, bnb = tiny, n_seeds = 50L))
  expect_false(cert$gap_certified)
  expect_identical(cert$gap_method, "multistart_bfgs")

  # Falling back is the planned outcome; being told why is not optional.
  err <- expect_error(
    certify(mn_mixture(), prob, bnb = tiny, n_seeds = 50L,
            require_certified = TRUE)
  )
  expect_match(conditionMessage(err), "max_coef")

  # With a sufficient budget the same call certifies.
  ok <- certify(mn_mixture(), prob, bnb = bnb_control(), n_seeds = 50L)
  expect_true(ok$gap_certified)
})

test_that("defaults are unchanged on a multinomial problem, and bnb only raises the gap", {
  prob <- mn_problem()
  proj <- mn_mixture()

  set.seed(17)
  base <- certify(proj, prob, n_seeds = 200L)
  set.seed(17)
  again <- certify(proj, prob, n_seeds = 200L)
  # The exact path is deterministic given the seed, and the default is the old
  # path verbatim.
  expect_identical(base$gap_used, again$gap_used)
  expect_false(base$gap_certified)
  expect_identical(base$gap_method, "multistart_bfgs")

  cert <- certify(proj, prob, bnb = bnb_control())
  expect_true(cert$gap_certified)
  expect_identical(cert$gap_method, "bernstein_bnb")
  expect_true(is.integer(cert$bnb_iterations))
  # A bound is at least as large as any lower bound on the same quantity.
  expect_gte(cert$gap_used, base$gap_used)
  # And everything derived from the gap still flows through the same formulas.
  expect_equal(cert$gap_used, cert$gap)
  expect_equal(cert$gap_se, 0)
  expect_equal(cert$kl_lower_bound, cert$kl - cert$gap_used)
  expect_equal(cert$growth_rate, cert$kl - log1p(cert$gap_used))
})

test_that("the guarantee level travels with the e_variable", {
  prob <- mn_problem()
  proj <- mn_mixture()

  ev_plain <- e_variable(
    numerator = prob$alternative,
    projection = as_marginal(proj, prob$family),
    gap = 0.1
  )
  expect_false(ev_plain@gap_certified)

  ev_cert <- e_variable(
    numerator = prob$alternative,
    projection = as_marginal(proj, prob$family),
    gap = 0.1,
    gap_certified = TRUE
  )
  expect_true(ev_cert@gap_certified)
  # The flag is provenance, not arithmetic: it must not touch the e-values.
  X <- support(prob$family)
  expect_identical(e_value(ev_plain, X), e_value(ev_cert, X))
})

test_that("a negative certified gap warns, since it can only be a bug", {
  prob <- mn_problem()
  proj <- as_marginal(mn_mixture(), prob$family)

  # sum_j w_j G(theta_j) = 1 forces sup G >= 1, so a certified gap below 0 is
  # not the numerical noise the floor was put there for.
  expect_warning(
    ev <- e_variable(numerator = prob$alternative, projection = proj,
                     gap = -1e-6, gap_certified = TRUE),
    "certified gap is negative"
  )
  expect_identical(ev@gap, 0)

  # Uncertified, the same value is exactly the noise the floor handles.
  expect_no_warning(
    ev2 <- e_variable(numerator = prob$alternative, projection = proj,
                      gap = -1e-6)
  )
  expect_identical(ev2@gap, 0)
})

test_that("run_ripr carries the certificate's guarantee level onto its e_variable", {
  prob <- mn_problem()
  ref <- c(0.6, 0.25, 0.15)
  faces <- prob$null

  fit <- function(...) {
    run_ripr(
      prob,
      init_atoms = init_on_faces(faces, ref),
      init_atom_faces = seq_along(faces),
      fw_iters = 5L,
      em_iters = 5L,
      n_seeds = 50L,
      verbose = FALSE,
      ...
    )
  }

  set.seed(42)
  plain <- fit()
  expect_false(plain$e_variable@gap_certified)
  expect_false(plain$certificate$gap_certified)

  set.seed(42)
  certified <- fit(certify_bnb = bnb_control())
  expect_true(certified$e_variable@gap_certified)
  expect_true(certified$certificate$gap_certified)
  # Certifying can only widen the correction, never narrow it.
  expect_gte(certified$e_variable@gap, plain$e_variable@gap)
})
