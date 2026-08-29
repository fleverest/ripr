# Properties of R/dirichlet.R.
#
# The Dirichlet marginal is a ratio of two integrals of the same shape, so the
# tests come in two groups: that the integral is exact where a closed form
# exists to check it against, and that the induced law is a probability
# distribution even where it is not.

k3_family <- function(n) multinomial_family(n_trials = n, k = 3L)

truncated <- function(alpha, region, ...) {
  truncated_dirichlet(alpha = alpha, region = region, ...)
}

plurality_part <- function(j, k = 3L) {
  vertices <- diag(k)
  vertices[, 1L] <- replace(numeric(k), c(1L, j), 0.5)
  simplex_region(vertices = vertices)
}

# `{theta : theta_1 >= theta_2, theta_3}`, where candidate 1 wins outright.
plurality_complement <- function(k = 3L) {
  family <- multinomial_family(n_trials = 1L, k = k)
  setdiff(
    family@parameter_space,
    union_region(lapply(2:k, plurality_part, k = k))
  )
}

# --- Concentrations -----------------------------------------------------------

test_that("only the truncated case requires integer concentrations", {
  # The restriction is a property of the quadrature, not of the Dirichlet: over
  # the whole simplex the integral is the Beta function, which `lgamma()`
  # evaluates exactly for any positive real.
  region <- plurality_complement()
  expect_error(truncated(c(2, 2.7, 1), region), "positive integer")
  expect_error(truncated(c(2, 2.7, 1), region), "singular")
  expect_identical(truncated(c(4, 3, 2), region)@alpha, c(4L, 3L, 2L))

  expect_true(S7::S7_inherits(dirichlet(c(1.5, 2)), dirichlet))
  expect_equal(dirichlet(c(0.5, 1.5, 2.25))@alpha, c(0.5, 1.5, 2.25))
})

test_that("a non-integer Dirichlet still induces the right mixture", {
  n <- 5L
  family <- k3_family(n)
  alpha <- c(0.5, 1.5, 2.25)
  outcomes <- enumerate_space(family@sample_space)

  hand <- apply(outcomes, 1L, function(x) {
    lgamma(n + 1) -
      sum(lgamma(x + 1)) +
      (sum(lgamma(alpha + x)) - lgamma(sum(alpha + x))) -
      (sum(lgamma(alpha)) - lgamma(sum(alpha)))
  })
  got <- log_density(family(dirichlet(alpha = alpha)), outcomes)

  expect_equal(got, hand, tolerance = 1e-14)
  expect_equal(sum(exp(got)), 1, tolerance = 1e-12)
})

test_that("alpha must be a positive vector of length at least two", {
  expect_error(dirichlet(alpha = 3), "at least 2 entries")
  expect_error(dirichlet(alpha = c(2, 0)), "finite positive")
  expect_error(dirichlet(alpha = c(2, -1)), "finite positive")
})

# --- The untruncated case, against closed forms -------------------------------

test_that("dirichlet induces the Dirichlet-multinomial", {
  n <- 6L
  family <- k3_family(n)
  alpha <- c(4, 3, 2)
  outcomes <- enumerate_space(family@sample_space)

  hand <- apply(outcomes, 1L, function(x) {
    lgamma(n + 1) -
      sum(lgamma(x + 1)) +
      (sum(lgamma(alpha + x)) - lgamma(sum(alpha + x))) -
      (sum(lgamma(alpha)) - lgamma(sum(alpha)))
  })
  got <- log_density(family(dirichlet(alpha = alpha)), outcomes)

  expect_equal(got, hand, tolerance = 1e-14)
  expect_equal(sum(exp(got)), 1, tolerance = 1e-12)
})

test_that("a uniform Dirichlet induces the uniform law over the lattice", {
  # Bose-Einstein: Dir(1, ..., 1) makes every one of the
  # `choose(n + K - 1, K - 1)` count vectors equally likely. It exercises the
  # whole coefficient path, since the multinomial coefficient is exactly what
  # the Beta ratio has to cancel.
  for (k in c(2L, 3L, 4L)) {
    n <- 5L
    family <- multinomial_family(n_trials = n, k = k)
    outcomes <- enumerate_space(family@sample_space)
    flat <- family(dirichlet(alpha = rep(1, k)))
    mass <- exp(log_density(flat, outcomes))

    expect_equal(nrow(outcomes), choose(n + k - 1L, k - 1L))
    expect_equal(
      mass,
      rep(1 / nrow(outcomes), nrow(outcomes)),
      tolerance = 1e-12
    )
  }
})

# --- The quadrature, against closed forms -------------------------------------

test_that("the reference rule integrates monomials on the simplex exactly", {
  # `int_T prod_j lambda_j^m_j / int_T 1 = (K-1)! prod m_j! / (|m|+K-1)!`
  monomial <- function(k, m) {
    rule <- ripr:::reference_simplex_rule(k, degree = sum(m))
    weight <- exp(rule$log_w)
    integral <- sum(weight * apply(rule$lambda, 1L, function(l) prod(l^m)))
    integral / sum(weight)
  }
  exact <- function(k, m) {
    factorial(k - 1L) * prod(factorial(m)) / factorial(sum(m) + k - 1L)
  }

  # The check from first principles: E[lambda_1 lambda_2] under Dir(1, 1, 1).
  expect_equal(monomial(3L, c(1, 1, 0)), 1 / 12, tolerance = 1e-12)

  cases <- list(
    list(k = 2L, m = c(5, 2)),
    list(k = 3L, m = c(4, 3, 2)),
    list(k = 4L, m = c(2, 1, 0, 3)),
    list(k = 4L, m = c(0, 0, 0, 0))
  )
  for (case in cases) {
    expect_equal(
      monomial(case$k, case$m),
      exact(case$k, case$m),
      tolerance = 1e-12
    )
  }
})

test_that("quadrature nodes lie strictly inside the simplex", {
  # This is what keeps `log(theta)` finite: an exponent of zero, which
  # `alpha_j == 1` with `x_j == 0` produces, would otherwise meet `-Inf`.
  for (k in c(2L, 3L, 4L)) {
    rule <- ripr:::reference_simplex_rule(k, degree = 9L)
    expect_true(all(rule$lambda > 0))
    expect_equal(rowSums(rule$lambda), rep(1, nrow(rule$lambda)))
    expect_true(all(is.finite(rule$log_w)))
  }
})

test_that("truncating to the whole simplex reproduces the untruncated law", {
  # Not to Monte Carlo error: with integer alpha the integrand is a polynomial
  # and the rule is exact for it, so anything looser than rounding would mean
  # the degree calculation is wrong.
  family <- k3_family(8L)
  outcomes <- enumerate_space(family@sample_space)
  alpha <- c(4, 3, 2)

  closed <- log_density(family(dirichlet(alpha = alpha)), outcomes)
  quadrature <- log_density(
    family(truncated(alpha, family@parameter_space)),
    outcomes
  )

  expect_equal(quadrature, closed, tolerance = 1e-10)
})

test_that("the cell integrals add up across a decomposition of the simplex", {
  # The one thing a single-cell check cannot see: whether each cell's weights
  # carry the right `abs(det(V))`. Splitting the simplex into pieces that are
  # neither the identity nor of equal volume, then summing, does.
  family <- k3_family(7L)
  alpha <- c(3, 2, 2)
  whole <- truncated(alpha, family@parameter_space)
  pieces <- lapply(
    parts(disjoin(union_region(
      plurality_part(2L),
      plurality_part(3L),
      plurality_complement()
    ))),
    function(cell) truncated(alpha, cell)
  )
  expect_gt(length(pieces), 2L)

  beta <- matrix(c(alpha + c(4, 2, 1), alpha), nrow = 3L)
  total <- rowSums(vapply(
    pieces,
    function(w) exp(ripr:::log_region_integral(w, beta)),
    numeric(ncol(beta))
  ))

  expect_equal(
    log(total),
    ripr:::log_region_integral(whole, beta),
    tolerance = 1e-10
  )
  # And the whole-simplex integral is the multivariate Beta function.
  expect_equal(
    ripr:::log_region_integral(whole, beta),
    colSums(lgamma(beta)) - lgamma(colSums(beta)),
    tolerance = 1e-10
  )
})

# --- Self-normalisation -------------------------------------------------------

test_that("the induced masses sum to one even under an under-degree rule", {
  # The property the whole design rests on. Any node set inside the simplex
  # telescopes, because the counts sum to `n` at every node, so quadrature
  # error distorts the shape of the distribution without making it improper.
  family <- k3_family(10L)
  outcomes <- enumerate_space(family@sample_space)
  region <- plurality_complement()

  exact <- truncated(c(4, 3, 2), region)
  coarse <- truncated(c(4, 3, 2), region, degree_slack = -100L)

  exact_mass <- exp(log_density(family(exact), outcomes))
  coarse_mass <- exp(log_density(family(coarse), outcomes))

  expect_equal(sum(coarse_mass), 1, tolerance = 1e-12)
  expect_true(all(coarse_mass >= 0))
  # The individual masses are wrong, which is the point: it is a distribution
  # regardless of accuracy, not a distribution because it is accurate.
  expect_gt(max(abs(coarse_mass - exact_mass)), 1e-3)
})

test_that("raising the degree beyond the exact one changes nothing", {
  family <- k3_family(6L)
  outcomes <- enumerate_space(family@sample_space)
  region <- plurality_complement()
  at <- function(slack) {
    log_density(
      family(truncated(c(2, 2, 3), region, degree_slack = slack)),
      outcomes
    )
  }
  expect_equal(at(6L), at(0L), tolerance = 1e-10)
})

# --- Additivity under subdivision ---------------------------------------------

# Join the edge midpoints: four congruent children, each a quarter of the area.
quadrisect <- function(v) {
  m <- cbind(
    (v[, 1] + v[, 2]) / 2,
    (v[, 1] + v[, 3]) / 2,
    (v[, 2] + v[, 3]) / 2
  )
  list(
    cbind(v[, 1], m[, 1], m[, 2]),
    cbind(v[, 2], m[, 1], m[, 3]),
    cbind(v[, 3], m[, 2], m[, 3]),
    m
  )
}

subdivide <- function(v, depth) {
  if (depth == 0L) {
    return(list(v))
  }
  unlist(
    lapply(quadrisect(v), subdivide, depth = depth - 1L),
    recursive = FALSE
  )
}

# The unnormalised integral over one cell, at `beta = alpha`, which is the
# measure the cell carries under `Dir(alpha)`.
cell_measure <- function(vertices, alpha) {
  w <- truncated(alpha, simplex_region(vertices = vertices))
  exp(log_region_integral(w, matrix(alpha, ncol = 1L)))
}

test_that("subdividing the simplex leaves the total measure unchanged", {
  # The rule is exact, so a decomposition has to sum back to the closed form to
  # rounding, however fine it is. This is the additive counterpart of the
  # self-normalisation test above: that one fixes the total over outcomes at
  # one `alpha`, this one fixes the total over cells at one outcome.
  for (alpha in list(c(1, 1, 1), c(4, 3, 2))) {
    exact <- exp(sum(lgamma(alpha)) - lgamma(sum(alpha)))
    for (depth in 0:3) {
      cells <- subdivide(diag(3), depth)
      total <- sum(vapply(cells, cell_measure, numeric(1), alpha = alpha))
      expect_equal(length(cells), 4^depth)
      expect_equal(total, exact, tolerance = 1e-13)
    }
  }
})

test_that("subdividing a truncated support leaves its measure unchanged", {
  medial <- cbind(c(.5, .5, 0), c(.5, 0, .5), c(0, .5, .5))
  alpha <- c(4, 3, 2)
  coarse <- cell_measure(medial, alpha)
  for (depth in 1:3) {
    cells <- subdivide(medial, depth)
    total <- sum(vapply(cells, cell_measure, numeric(1), alpha = alpha))
    expect_equal(total, coarse, tolerance = 1e-13)
  }
})

test_that("cells of degenerate aspect ratio do not degrade the total", {
  # Quadrisection keeps every child similar to its parent, so it never tests
  # conditioning. A fan to geometrically spaced points on the opposite edge
  # does: cell `i` has base `2^-i` against a fixed height, so the vertex
  # matrix approaches singularity while the union stays the whole simplex.
  m <- 20L
  s <- c(1 - 2^-(0:m), 1)
  q <- rbind(0, 1 - s, s)
  slivers <- lapply(
    seq_len(m + 1L),
    function(i) cbind(c(1, 0, 0), q[, i], q[, i + 1L])
  )
  expect_gt(max(vapply(slivers, kappa, numeric(1))), 1e6)
  expect_lt(min(vapply(slivers, function(v) abs(det(v)), numeric(1))), 1e-6)

  for (alpha in list(c(1, 1, 1), c(4, 3, 2))) {
    total <- sum(vapply(slivers, cell_measure, numeric(1), alpha = alpha))
    expect_equal(
      total,
      exp(sum(lgamma(alpha)) - lgamma(sum(alpha))),
      tolerance = 1e-13
    )
  }
})

# --- Symmetry -----------------------------------------------------------------

test_that("a symmetric prior over mirrored cells gives permuted densities", {
  family <- k3_family(8L)
  outcomes <- enumerate_space(family@sample_space)
  swapped <- outcomes[, c(2L, 1L, 3L)]

  # `{theta_1 <= theta_2}` and its image under swapping the first two
  # coordinates, `{theta_2 <= theta_1}`.
  lower <- plurality_part(2L)
  upper <- simplex_region(
    vertices = cbind(c(0.5, 0.5, 0), c(1, 0, 0), c(0, 0, 1))
  )

  alpha <- rep(2, 3)
  a <- log_density(family(truncated(alpha, lower)), outcomes)
  b <- log_density(family(truncated(alpha, upper)), swapped)

  expect_equal(a, b, tolerance = 1e-10)
})

# --- Regions the measure refuses ----------------------------------------------

test_that("a region with no full-dimensional cell at all is refused", {
  segment <- simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 0.5, 0.5)))
  point <- point_region(theta = c(0.5, 0.3, 0.2))
  expect_error(truncated(c(2, 2, 2), segment), "measure zero")
  expect_error(truncated(c(2, 2, 2), point), "measure zero")
})

test_that("lower-dimensional cells are dropped with a warning, not refused", {
  # A degenerate piece integrates to zero, so it cannot change an answer --
  # only waste the caller's assumption that it counted for something.
  medial <- simplex_region(
    vertices = cbind(c(.5, .5, 0), c(.5, 0, .5), c(0, .5, .5))
  )
  sliver <- simplex_region(vertices = cbind(c(1, 0, 0), c(0.9, 0.1, 0)))
  mixed <- union_region(medial, sliver)

  expect_warning(
    w <- truncated(c(4, 3, 2), mixed),
    class = "ripr_degenerate_warning"
  )
  expect_warning(truncated(c(4, 3, 2), mixed), "measure zero")
  expect_length(w@cells, 1L)

  # And the surviving cell carries exactly the measure the whole region had.
  clean <- truncated(c(4, 3, 2), medial)
  shape <- matrix(c(4, 3, 2), ncol = 1L)
  expect_equal(
    log_region_integral(w, shape),
    log_region_integral(clean, shape),
    tolerance = 1e-13
  )
})

test_that("a region outside or beyond the simplex is refused", {
  tetrahedron <- simplex_region(
    vertices = cbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
  )
  expect_error(
    truncated(c(2, 2, 2), halfspace_region(normal = c(1, -1, 0))),
    "unbounded"
  )
  # Full-dimensional in R^3 but off the simplex, so a different complaint
  # from the measure-zero one above.
  expect_error(truncated(c(2, 2, 2), tetrahedron), "probability simplex")
  expect_error(
    truncated(c(2, 2), plurality_complement()),
    "3 dimensions but `alpha` has 2"
  )
})

test_that("an oversized quadrature rule is refused, naming what drove it", {
  family <- multinomial_family(n_trials = 60L, k = 5L)
  mixing <- truncated(rep(2, 5), family@parameter_space, max_nodes = 100)
  expect_error(
    log_density(family(mixing), matrix(12L, nrow = 1L, ncol = 5L)),
    "above `max_nodes`"
  )
  expect_error(
    log_density(family(mixing), matrix(12L, nrow = 1L, ncol = 5L)),
    "n_trials = 60"
  )
})

# --- Sampling -----------------------------------------------------------------

test_that("draw lands inside the region", {
  set.seed(1)
  region <- plurality_complement()
  theta <- draw(truncated(c(4, 3, 2), region), 200L)

  expect_equal(dim(theta), c(200L, 3L))
  expect_true(all(apply(theta, 1L, function(t) contains(region, t))))
  expect_equal(rowSums(theta), rep(1, 200L))
})

test_that("draw on an untruncated Dirichlet lands in the simplex", {
  set.seed(1)
  theta <- draw(dirichlet(alpha = c(4, 3, 2)), 100L)
  expect_equal(rowSums(theta), rep(1, 100L))
  expect_true(all(theta > 0))
})

test_that("a region the prior barely reaches errors rather than looping", {
  set.seed(1)
  corner <- simplex_region(vertices = cbind(
    c(1, 0, 0), c(0.999, 0.001, 0), c(0.999, 0, 0.001)
  ))
  pinched <- truncated(c(1, 1, 1), corner)
  expect_error(draw(pinched, 10L), "rejection sampling")
  expect_error(draw(pinched, 10L), "proposals")
})

test_that("draw goes through the kernel one parameter at a time", {
  set.seed(1)
  family <- k3_family(9L)
  Q <- family(truncated(c(4, 3, 2), plurality_complement()))
  x <- draw(Q, 20L)

  expect_equal(dim(x), c(20L, 3L))
  expect_equal(rowSums(x), rep(9, 20L))
})

# --- The mixing-measure interface ---------------------------------------------

test_that("a continuous mixing measure has no support to list", {
  for (mixing in list(
    dirichlet(alpha = c(4, 3, 2)),
    truncated(c(4, 3, 2), plurality_complement())
  )) {
    expect_true(is.na(n_atoms(mixing)))
    name <- attr(S7::S7_class(mixing), "name")
    expect_error(atoms(mixing), name)
    expect_error(weights(mixing), name)
    expect_error(atoms(mixing), "density rather than a support")
  }
})

test_that("reference_point seeds the optimiser from inside the region", {
  # The mode when every concentration exceeds 1, the mean otherwise, projected
  # onto the region when the untruncated point is outside it.
  expect_equal(reference_point(dirichlet(c(4, 3, 2))), c(3, 2, 1) / 6)
  expect_equal(reference_point(dirichlet(c(1, 3, 2))), c(1, 3, 2) / 6)

  region <- plurality_complement()
  inside <- reference_point(truncated(c(4, 3, 2), region))
  expect_true(contains(region, inside))

  # A prior whose mode sits in the null still has to start somewhere legal.
  outside <- reference_point(truncated(c(2, 6, 6), region))
  expect_true(contains(region, outside))
})

test_that("a Dirichlet mixed through the wrong family errors naming both", {
  wrong <- gaussian_family(dim = 3L)(dirichlet(alpha = c(2, 2, 2)))
  expect_error(
    log_density(wrong, c(0, 0, 0)),
    "`dirichlet` over a `gaussian_family`"
  )
  # And it says there is no silent Monte Carlo fallback, and what to reach for
  # instead of one.
  expect_error(
    log_density(wrong, c(0, 0, 0)),
    "not approximated by Monte Carlo"
  )
  expect_error(log_density(wrong, c(0, 0, 0)), "discretise")
})

test_that("the concentration count must match the family's categories", {
  # Caught when the mixture is built, not when it is first evaluated: a
  # `Dir(2, 2)` lives in 2 dimensions and a 3-category family's parameters in 3.
  expect_error(
    k3_family(4L)(dirichlet(c(2, 2))),
    "over 2 dimensions but the family's parameters have 3"
  )
  # The density method keeps its own check, now only reachable by calling it
  # directly, since no mixture over a mismatched pair can be constructed.
  expect_error(
    ripr:::mixture_log_density(dirichlet(c(2, 2)), k3_family(4L), c(2L, 1L, 1L)),
    "2 entries but the family has 3 categories"
  )
})

# --- Plumbing into a fit ------------------------------------------------------

test_that("an exact engine resolves against a truncated Dirichlet", {
  family <- k3_family(8L)
  Q <- family(truncated(c(4, 3, 2), plurality_complement()))
  engine <- resolve_engine(exact_engine(), Q, family)

  expect_true(deterministic(engine))
  expect_equal(sum(exp(engine@log_w)), 1, tolerance = 1e-12)
  expect_equal(n_nodes(engine), nrow(enumerate_space(family@sample_space)))
})

test_that("gauss-hermite refuses a Dirichlet alternative", {
  family <- k3_family(8L)
  Q <- family(truncated(c(4, 3, 2), plurality_complement()))
  expect_error(resolve_engine(gh_engine(5L), Q, family), "Gaussian alternative")
})

test_that("a truncated Dirichlet alternative fits against the plurality null", {
  set.seed(1)
  k <- 3L
  family <- multinomial_family(n_trials = 10L, k = k)
  plurality <- union_region(lapply(2:k, plurality_part, k = k))
  alternative <- setdiff(family@parameter_space, plurality)

  Q <- family(truncated(c(4, 3, 2), alternative))
  state <- ripr_init(Q, null_model(family, plurality), engine = exact_engine())
  state <- fw_step(state, times = 20L)
  state <- em_step(state, times = 60L, until = kl_flat(1e-9))
  fit <- ripr_finish(state, reoptimise = TRUE, identify = TRUE)

  expect_true(is.finite(fit$kl))
  expect_gt(fit$kl, 0)
  # Both verbs are monotone in KL, so a rise would mean the alternative's
  # density and the engine's weights had come apart.
  expect_true(all(diff(state@trace$kl) <= 1e-12))
  expect_true(kl_flat(1e-5)(state))
  expect_lte(fit$kl, state@trace$kl[[1L]])
})
