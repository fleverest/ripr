# Properties of R/certify.R.
#
# The only property that really matters is one-sidedness: `sup_ub` must never be
# below the true supremum.
#
# The reference is a dense barycentric grid over each subnull. This is a *lower*
# bound on the true supremum, so `sup_ub >= grid` is a necessary condition but
# not sufficient. It doesn't prove validity enclosure, but may catch something
# going horribly wrong with certification.

# The plurality null H_0 = union_j {theta : theta_1 <= theta_j}, as simplices.
plurality_subnulls <- function(k) {
  lapply(2:k, function(j) {
    vertices <- diag(k)
    vertices[, 1L] <- replace(numeric(k), c(1L, j), 0.5)
    simplex_null(vertices = vertices)
  })
}

plurality_null <- function(n, k) {
  family <- multinomial_family(n_trials = n, k = k)
  null_model(family, plurality_subnulls(k))
}

# A random variable pinned to given values on the support, so the tests do not
# depend on a fit having converged to anything in particular.
tabulated_rv <- function(family, values) {
  outcomes <- support(family)
  key <- apply(outcomes, 1L, paste, collapse = "-")
  force(values)
  random_variable(
    function(x) {
      values[match(apply(as.matrix(x), 1L, paste, collapse = "-"), key)]
    },
    family = family,
    label = "<tabulated>"
  )
}

# max over a dense barycentric grid on one facet: a lower bound on the truth.
facet_grid_max <- function(family, values, vertices, m = 60L) {
  k <- ncol(vertices)
  weights <- compositions(m, k) / m
  theta <- weights %*% t(vertices)
  outcomes <- support(family)
  max(as.vector(crossprod(
    exp(log_density_batch(family, t(theta), outcomes)),
    values
  )))
}

null_grid_max <- function(null, values, m = 60L) {
  max(vapply(
    null@subnulls,
    function(s) facet_grid_max(null@family, values, s@vertices, m),
    numeric(1L)
  ))
}

# --- One-sidedness ------------------------------------------------------------

test_that("sup_ub is never below a dense grid search", {
  set.seed(101)
  null <- plurality_null(n = 8L, k = 3L)
  for (rep in 1:8) {
    values <- stats::runif(nrow(support(null@family)), 0, 10)
    res <- certify(tabulated_rv(null@family, values), null, tol = 1e-9)
    expect_gte(res$sup_ub, null_grid_max(null, values))
  }
})

test_that("sup_ub stays valid when the node budget is exhausted", {
  # Validity does not depend on convergence. A run cut off after one bisection
  # returns a loose bound, not an invalid one -- this is the property that lets
  # `certify()` be interrupted.
  set.seed(102)
  null <- plurality_null(n = 8L, k = 3L)
  values <- stats::runif(nrow(support(null@family)), 0, 10)
  x <- tabulated_rv(null@family, values)
  truth <- null_grid_max(null, values)

  bounds <- vapply(
    c(1L, 2L, 5L, 20L, 500L),
    function(m) certify(x, null, tol = 0, max_nodes = m)$sup_ub,
    numeric(1L)
  )
  expect_true(all(bounds >= truth))
  expect_false(is.unsorted(rev(bounds))) # non-increasing in the budget
})

test_that("sup_ub brackets sup_lb", {
  set.seed(103)
  null <- plurality_null(n = 10L, k = 3L)
  values <- stats::runif(nrow(support(null@family)), 0, 10)
  x <- tabulated_rv(null@family, values)
  res <- certify(x, null, tol = 1e-9)
  expect_gte(res$sup_ub, res$sup_lb)
  # The searched lower bound is a different algorithm on the same problem, so
  # this is the cross-check between the two halves of the file.
  expect_gte(
    res$sup_ub,
    sup_lb(x, null, n_seeds = 100L, n_restarts = 10L)$sup_lb
  )
})

test_that("dividing by the bound gives an e-variable", {
  # The reason the function exists. `X / sup_ub` must have null expectation at
  # most 1 everywhere on H_0.
  set.seed(104)
  null <- plurality_null(n = 8L, k = 4L)
  family <- null@family
  outcomes <- support(family)
  values <- stats::runif(nrow(outcomes), 0, 5)
  x <- tabulated_rv(family, values)
  res <- certify(x, null, tol = 1e-9)

  # The e-variable:
  e <- x / res$sup_ub
  e_vals <- e(outcomes)

  expect_equal(e_vals, values / res$sup_ub)

  # E_theta[e] for theta supplied as columns.
  expectations <- function(theta) {
    as.vector(crossprod(
      exp(log_density_batch(family, theta, outcomes)),
      e_vals
    ))
  }

  for (s in null@subnulls) {
    weights <- matrix(stats::rgamma(4L * 200L, shape = 1), nrow = 4L)
    theta <- s@vertices %*% div_by_col(weights, colSums(weights))
    expect_lte(max(expectations(theta)), 1)
  }

  tight <- c(
    vapply(
      null@subnulls,
      function(s) max(expectations(s@vertices)),
      numeric(1L)
    ),
    expectations(matrix(
      sup_lb(x, null, n_seeds = 200L, n_restarts = 20L)$theta,
      ncol = 1L
    ))
  )
  expect_lte(max(tight), 1)
  expect_gt(max(tight), 1 - 1e-3)
  expect_gt(res$sup_lb / res$sup_ub, 1 - 1e-6)
  expect_lt(max(tight), 1)
  expect_lt(res$sup_lb / res$sup_ub, 1)
})

test_that("a constant variable certifies to its own value", {
  null <- plurality_null(n = 6L, k = 3L)
  x <- tabulated_rv(null@family, rep(3.5, nrow(support(null@family))))
  res <- certify(x, null, tol = 1e-9)
  expect_equal(res$sup_lb, 3.5)
  expect_lt(res$sup_ub - 3.5, 1e-12)
  expect_true(all(res$nodes == 0L))
})

test_that("a variable maximised at a facet vertex needs no subdivision", {
  # The enclosure is exact at the vertices (PBP 10.2), so if the maximiser is
  # a vertex the starting point is already a tight bound.
  n <- 6L
  null <- plurality_null(n = n, k = 3L)
  family <- null@family
  outcomes <- support(family)
  q <- c(0.1, 0.8, 0.1)
  values <- exp(as.vector(outcomes %*% (log(q) - log(rep(1 / 3, 3)))))
  res <- certify(tabulated_rv(family, values), null, tol = 1e-9)

  expect_true(all(res$nodes == 0L))
  expect_equal(res$sup_ub, res$sup_lb, tolerance = 1e-9)

  expect_equal(res$sup_ub, (3 * max(q))^n, tolerance = 1e-9)
  expect_equal(
    res$bounds,
    rep((3 * max(q))^n, length(null@subnulls)),
    tolerance = 1e-9
  )

  expect_gt(max(values) / min(values), 1e5)
})

test_that("a variable maximised in a facet interior does need subdivision", {
  set.seed(110)
  null <- plurality_null(n = 8L, k = 3L)
  outcomes <- support(null@family)
  values <- stats::runif(nrow(outcomes), 0, 10)
  res <- certify(tabulated_rv(null@family, values), null, tol = 1e-9)

  expect_true(all(res$nodes > 0L))
  # And the work bought something: the bound is below the seed node's, which is
  # the enclosure over the whole facet before any subdivision.
  expect_lt(
    res$sup_ub,
    max(vapply(
      null@subnulls,
      function(s) {
        certify_sup(
          list(list(V = s@vertices, coef = values)),
          bernstein_lattice(8L, 3L),
          tol = 1e-9,
          max_iter = 0L
        )$bound
      },
      numeric(1L)
    ))
  )
})

# --- The pmf is the Bernstein basis -------------------------------------------

test_that("the certified bound is a bound on the expectation itself", {
  # The same correspondence one level up, stated in the terms the caller cares
  # about: `sup_ub` bounds E_theta[X] computed from the family, not merely the
  # polynomial the bounding method happened to be handed.
  set.seed(108)
  null <- plurality_null(n = 6L, k = 3L)
  family <- null@family
  outcomes <- support(family)
  values <- stats::runif(nrow(outcomes), 0, 10)
  res <- certify(tabulated_rv(family, values), null, tol = 1e-9)

  for (s in null@subnulls) {
    weights <- matrix(stats::rgamma(3L * 300L, shape = 1), nrow = 3L)
    theta <- s@vertices %*% div_by_col(weights, colSums(weights))
    expectations <- as.vector(
      crossprod(exp(log_density_batch(family, theta, outcomes)), values)
    )
    expect_lte(max(expectations), res$sup_ub)
  }
})

# --- Return shape -------------------------------------------------------------

test_that("certify() reports one entry per subnull", {
  set.seed(105)
  null <- plurality_null(n = 6L, k = 4L)
  values <- stats::runif(nrow(support(null@family)), 0, 10)
  res <- certify(tabulated_rv(null@family, values), null, tol = 1e-6)

  n_sub <- length(null@subnulls)
  expect_length(res$bounds, n_sub)
  expect_length(res$incumbents, n_sub)
  expect_length(res$nodes, n_sub)
  expect_length(res$converged, n_sub)
  expect_length(res$budget_hit, n_sub)
  expect_type(res$nodes, "integer")
  expect_type(res$converged, "logical")
  expect_type(res$budget_hit, "logical")
  expect_identical(res$method, "bernstein")
  expect_equal(res$sup_ub, max(res$bounds))
  expect_equal(res$sup_lb, max(res$incumbents))
})

test_that("certify() carries back the variable and null it holds for", {
  # A certificate that does not name what it certifies is not a certificate.
  null <- plurality_null(n = 4L, k = 3L)
  x <- tabulated_rv(null@family, rep(1, nrow(support(null@family))))
  res <- certify(x, null)
  expect_identical(res$random_variable, x)
  expect_identical(res$null, null)
})

test_that("converged and budget_hit distinguish the two ways of stopping", {
  set.seed(111)
  null <- plurality_null(n = 8L, k = 3L)
  values <- stats::runif(nrow(support(null@family)), 0, 10)
  x <- tabulated_rv(null@family, values)

  tight <- certify(x, null, tol = 1e-12, max_nodes = 5000L)
  expect_true(all(tight$converged))
  expect_false(any(tight$budget_hit))

  starved <- certify(x, null, tol = 1e-12, max_nodes = 2L)
  expect_true(any(starved$budget_hit))
  # Mutually exclusive per subnull: a search stops one way or the other.
  expect_false(any(starved$converged & starved$budget_hit))
  # Every subnull that ran out of budget used all of it.
  expect_true(all(starved$nodes[starved$budget_hit] == 2L))

  # The starved bound is still valid, just looser -- which is the whole reason
  # the distinction is worth reporting rather than erroring on.
  expect_gte(starved$sup_ub, tight$sup_ub)
})

# --- Registry dispatch --------------------------------------------------------

test_that("certify() sends each subnull to the bound_fn that claimed it", {
  # The reason the registry exists. With one entry the grouping is trivially
  # the identity, so this stubs a second entry to check that results are not
  # merely produced but land in the right slots.
  seen <- list()
  fake_bound_fn <- function(x, family, subnulls, control) {
    seen[[length(seen) + 1L]] <<- vapply(subnulls, class_name, character(1L))
    lapply(seq_along(subnulls), function(i) {
      list(
        bound = 99,
        incumbent = 99,
        theta = NULL,
        iterations = 0L,
        converged = TRUE,
        budget_hit = FALSE
      )
    })
  }
  # Override certify_methods for this test to add our fake one
  local_mocked_bindings(
    certify_methods = function() {
      list(
        list(
          name = "bernstein",
          family = multinomial_family,
          subnull = simplex_null,
          bound_fn = bernstein_bound,
          description = "real"
        ),
        list(
          name = "fake",
          family = multinomial_family,
          subnull = halfspace_null,
          bound_fn = fake_bound_fn,
          description = "stub"
        )
      )
    }
  )

  family <- multinomial_family(n_trials = 4L, k = 3L)
  facet <- diag(3L)
  facet[, 1L] <- c(0.5, 0.5, 0)
  null <- null_model(
    family,
    list(
      halfspace_null(normal = c(1, -1), offset = 0),
      simplex_null(vertices = facet),
      halfspace_null(normal = c(0, 1, -1), offset = 0)
    )
  )
  x <- tabulated_rv(family, rep(1, nrow(support(family))))
  res <- certify(x, null, tol = 1e-9)

  # Subnulls 1 and 3 went to the stub, subnull 2 to the real bound_fn
  expect_equal(res$bounds[c(1L, 3L)], c(99, 99))
  expect_lt(res$bounds[[2L]], 99)
  expect_setequal(res$method, c("bernstein", "fake"))
  # The stub was called once, with both of its subnulls together.
  expect_length(seen, 1L)
  expect_identical(seen[[1L]], c("halfspace_null", "halfspace_null"))
})

test_that("certify() rejects a bound_fn that returns the wrong number of results", {
  local_mocked_bindings(
    certify_methods = function() {
      list(list(
        name = "short",
        family = multinomial_family,
        subnull = simplex_null,
        bound_fn = function(x, family, subnulls, control) list(),
        description = "stub"
      ))
    }
  )
  null <- plurality_null(n = 4L, k = 3L)
  x <- tabulated_rv(null@family, rep(1, nrow(support(null@family))))
  expect_error(certify(x, null), "returned 0 results for 2 subnulls")
})

test_that("check_bound_result() names the field and the bound_fn", {
  ok <- list(
    bound = 1,
    incumbent = 1,
    iterations = 0L,
    converged = TRUE,
    budget_hit = FALSE
  )
  expect_silent(check_bound_result(list(ok), "demo", 1L))

  drop_field <- function(field) {
    bad <- ok
    bad[[field]] <- NULL
    check_bound_result(list(bad), "demo", 1L)
  }
  for (field in c(
    "bound",
    "incumbent",
    "iterations",
    "converged",
    "budget_hit"
  )) {
    expect_error(drop_field(field), field)
    expect_error(drop_field(field), "demo")
  }

  # The state no search can reach, which `isTRUE(NULL)` used to manufacture.
  neither <- ok
  neither$converged <- FALSE
  expect_error(
    check_bound_result(list(neither), "demo", 1L),
    "without recording why"
  )
})

# --- Refusals -----------------------------------------------------------------

test_that("certify() refuses a family and geometry it has no method for", {
  family <- gaussian_family(dim = 2L)
  null <- null_model(
    family,
    list(halfspace_null(normal = c(1, -1), offset = 0))
  )
  x <- random_variable(function(x) rep(1, nrow(as.matrix(x))), family = family)
  expect_error(certify(x, null), "No bounding method is implemented")
  expect_error(certify(x, null), "gaussian_family")
  expect_error(certify(x, null), "halfspace_null")
})

test_that("the refusal names the subnull, not the null model", {
  # The message is meant to say which geometry is missing a bound. Naming the
  # container instead makes it useless.
  family <- gaussian_family(dim = 2L)
  null <- null_model(
    family,
    list(halfspace_null(normal = c(1, -1), offset = 0))
  )
  x <- random_variable(function(x) rep(1, nrow(as.matrix(x))), family = family)
  msg <- tryCatch(certify(x, null), error = conditionMessage)
  expect_false(grepl("null_model", msg, fixed = TRUE))
  expect_false(grepl("FALSE", msg, fixed = TRUE))
})

test_that("the refusal is not repeated once per subnull", {
  family <- gaussian_family(dim = 2L)
  null <- null_model(
    family,
    list(
      halfspace_null(normal = c(1, -1), offset = 0),
      halfspace_null(normal = c(1, 0), offset = 0)
    )
  )
  x <- random_variable(function(x) rep(1, nrow(as.matrix(x))), family = family)
  msg <- tryCatch(certify(x, null), error = conditionMessage)
  expect_identical(
    lengths(regmatches(msg, gregexpr("No bounding method", msg))),
    1L
  )
})

test_that("certify() refuses a variable that is not finite on the support", {
  # An infinite value is legitimate for a likelihood ratio and fatal for a
  # bound: the supremum is then unbounded and no finite certificate exists.
  null <- plurality_null(n = 4L, k = 3L)
  values <- rep(1, nrow(support(null@family)))
  values[3L] <- Inf
  expect_error(
    certify(tabulated_rv(null@family, values), null),
    "not finite everywhere"
  )
})

test_that("certify() refuses a lattice above the coefficient budget", {
  null <- plurality_null(n = 8L, k = 3L)
  x <- tabulated_rv(null@family, rep(1, nrow(support(null@family))))
  expect_error(
    certify(x, null, max_coefficients = 10L),
    "above `max_coefficients`"
  )
  # The refusal must be raised before any of the work is done.
  expect_error(certify(x, null, max_coefficients = 10L), "would fit")
})

test_that("certify() rejects a non-random_variable and bad control values", {
  null <- plurality_null(n = 4L, k = 3L)
  x <- tabulated_rv(null@family, rep(1, nrow(support(null@family))))
  expect_error(certify(function(x) 1, null), "must be a `random_variable`")
  expect_error(certify(x, null, tol = -1))
  expect_error(certify(x, null, max_nodes = 0))
  expect_error(certify(x, null, max_coefficients = 0))
})

# --- sup_lb -------------------------------------------------------------------

test_that("sup_lb() reports a value the objective actually attains", {
  set.seed(106)
  null <- plurality_null(n = 8L, k = 3L)
  family <- null@family
  outcomes <- support(family)
  values <- stats::runif(nrow(outcomes), 0, 10)
  x <- tabulated_rv(family, values)
  found <- sup_lb(x, null, n_seeds = 100L, n_restarts = 10L)

  # E_theta[X] for theta supplied as columns.
  expectations <- function(theta) {
    as.vector(crossprod(
      exp(log_density_batch(family, theta, outcomes)),
      values
    ))
  }

  expect_equal(found$sup_lb, expectations(matrix(found$theta, ncol = 1L)))

  expect_true(contains(null@subnulls[[found$subnull]], found$theta))

  vertex_best <- max(vapply(
    null@subnulls,
    function(s) max(expectations(s@vertices)),
    numeric(1L)
  ))
  expect_gt(found$sup_lb, vertex_best)

  expect_gt(found$sup_lb / certify(x, null, tol = 1e-12)$sup_ub, 0.999)
})

test_that("sup_lb() improves on a single seed given more of them", {
  set.seed(109)
  null <- plurality_null(n = 8L, k = 3L)
  outcomes <- support(null@family)

  ratios <- vapply(
    seq_len(8L),
    function(i) {
      x <- tabulated_rv(null@family, stats::runif(nrow(outcomes), 0, 10))
      one <- sup_lb(x, null, n_seeds = 1L, n_restarts = 1L)$sup_lb
      many <- sup_lb(x, null, n_seeds = 100L, n_restarts = 10L)$sup_lb
      many / one
    },
    numeric(1L)
  )

  expect_gte(mean(ratios), 1)
  expect_true(all(ratios > 0.999))
})

test_that("sup_lb() rejects a non-random_variable", {
  null <- plurality_null(n = 4L, k = 3L)
  expect_error(sup_lb(function(x) 1, null), "must be a `random_variable`")
})
