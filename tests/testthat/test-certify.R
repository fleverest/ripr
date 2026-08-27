# Properties of R/certify.R.
#
# The only property that really matters is one-sidedness: `sup_ub` must never be
# below the true supremum.
#
# The reference is a dense barycentric grid over each part. This is a *lower*
# bound on the true supremum, so `sup_ub >= grid` is a necessary condition but
# not sufficient. It doesn't prove validity enclosure, but may catch something
# going horribly wrong with certification.

# The plurality null H_0 = union_j {theta : theta_1 <= theta_j}, as simplices.
plurality_parts <- function(k) {
  lapply(2:k, function(j) {
    vertices <- diag(k)
    vertices[, 1L] <- replace(numeric(k), c(1L, j), 0.5)
    simplex_region(vertices = vertices)
  })
}

plurality_null <- function(n, k) {
  family <- multinomial_family(n_trials = n, k = k)
  null_model(family, plurality_parts(k))
}

# A null that is convex but not a simplex: the square
# `{theta_1 <= 1/2, theta_2 <= 1/2}` cut out of the 2-simplex. Every vertex is
# an exact double summing to exactly 1, so the hull really does lie in the
# standard simplex rather than an ulp above it, and it triangulates into two
# certifiable cells.
simplex_square <- function() {
  polytope_region(
    vertices = cbind(
      c(0.5, 0.5, 0),
      c(0, 0.5, 0.5),
      c(0, 0, 1),
      c(0.5, 0, 0.5)
    )
  )
}

# A random variable pinned to given values on the support, so the tests do not
# depend on a fit having converged to anything in particular.
tabulated_rv <- function(family, values) {
  outcomes <- enumerate_space(family@sample_space)
  key <- apply(outcomes, 1L, paste, collapse = "-")
  force(values)
  random_variable(
    function(x) {
      values[match(apply(as.matrix(x), 1L, paste, collapse = "-"), key)]
    },
    sample_space = family@sample_space,
    label = "<tabulated>"
  )
}

# max over a dense barycentric grid on one facet: a lower bound on the truth.
facet_grid_max <- function(family, values, vertices, m = 60L) {
  k <- ncol(vertices)
  weights <- compositions(m, k) / m
  theta <- weights %*% t(vertices)
  outcomes <- enumerate_space(family@sample_space)
  max(as.vector(crossprod(
    exp(kernel_loglik_batch(family, t(theta), outcomes)),
    values
  )))
}

null_grid_max <- function(null, values, m = 60L) {
  max(vapply(
    parts(null@region),
    function(s) facet_grid_max(null@family, values, s@vertices, m),
    numeric(1L)
  ))
}

# --- One-sidedness ------------------------------------------------------------

test_that("sup_ub is never below a dense grid search", {
  set.seed(101)
  null <- plurality_null(n = 8L, k = 3L)
  for (rep in 1:8) {
    values <- stats::runif(
      nrow(enumerate_space(null@family@sample_space)),
      0,
      10
    )
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
  values <- stats::runif(nrow(enumerate_space(null@family@sample_space)), 0, 10)
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
  values <- stats::runif(nrow(enumerate_space(null@family@sample_space)), 0, 10)
  x <- tabulated_rv(null@family, values)
  res <- certify(x, null, tol = 1e-9)
  expect_gte(res$sup_ub, res$sup_lb)
  # The searched lower bound is a different algorithm on the same problem, so
  # this is the cross-check between the two halves of the file. SLSQP may attain
  # the supremum, so if the two meet one may round higher than the other, hence
  # the floating-point slack.
  expect_gte(
    res$sup_ub + 1e-8,
    sup_lb(x, null, n_seeds = 100L, n_restarts = 10L)$sup_lb
  )
})

test_that("dividing by the bound gives an e-variable", {
  # The reason the function exists. `X / sup_ub` must have null expectation at
  # most 1 everywhere on H_0.
  set.seed(104)
  null <- plurality_null(n = 8L, k = 4L)
  family <- null@family
  outcomes <- enumerate_space(family@sample_space)
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
      exp(kernel_loglik_batch(family, theta, outcomes)),
      e_vals
    ))
  }

  for (s in parts(null@region)) {
    weights <- matrix(stats::rgamma(4L * 200L, shape = 1), nrow = 4L)
    theta <- s@vertices %*% div_by_col(weights, colSums(weights))
    expect_lte(max(expectations(theta)), 1)
  }

  tight <- c(
    vapply(
      parts(null@region),
      function(s) max(expectations(s@vertices)),
      numeric(1L)
    ),
    expectations(matrix(
      sup_lb(x, null, n_seeds = 200L, n_restarts = 20L)$theta,
      ncol = 1L
    ))
  )
  expect_lte(max(tight), 1 + 1e-9)
  expect_gt(max(tight), 1 - 1e-3)
  expect_gt(res$sup_lb / res$sup_ub, 1 - 1e-6)
  expect_lte(res$sup_lb / res$sup_ub, 1 + 1e-9)
})

test_that("a constant variable certifies to its own value", {
  null <- plurality_null(n = 6L, k = 3L)
  x <- tabulated_rv(
    null@family,
    rep(3.5, nrow(enumerate_space(null@family@sample_space)))
  )
  res <- certify(x, null, tol = 1e-9)
  expect_equal(res$sup_lb, 3.5)
  expect_lt(res$sup_ub - 3.5, 1e-12)
  expect_true(all(res$iterations == 0L))
})

test_that("a variable maximised at a facet vertex needs no subdivision", {
  # The enclosure is exact at the vertices (PBP 10.2), so if the maximiser is
  # a vertex the starting point is already a tight bound.
  n <- 6L
  null <- plurality_null(n = n, k = 3L)
  family <- null@family
  outcomes <- enumerate_space(family@sample_space)
  q <- c(0.1, 0.8, 0.1)
  values <- exp(as.vector(outcomes %*% (log(q) - log(rep(1 / 3, 3)))))
  res <- certify(tabulated_rv(family, values), null, tol = 1e-9)

  expect_true(all(res$iterations == 0L))
  expect_equal(res$sup_ub, res$sup_lb, tolerance = 1e-9)

  expect_equal(res$sup_ub, (3 * max(q))^n, tolerance = 1e-9)
  expect_equal(
    res$bounds,
    rep((3 * max(q))^n, length(parts(null@region))),
    tolerance = 1e-9
  )

  expect_gt(max(values) / min(values), 1e5)
})

test_that("a variable maximised in a facet interior does need subdivision", {
  set.seed(110)
  null <- plurality_null(n = 8L, k = 3L)
  outcomes <- enumerate_space(null@family@sample_space)
  values <- stats::runif(nrow(outcomes), 0, 10)
  res <- certify(tabulated_rv(null@family, values), null, tol = 1e-9)

  expect_true(all(res$iterations > 0L))
  # And the work bought something: the bound is below the seed node's, which is
  # the enclosure over the whole facet before any subdivision.
  expect_lt(
    res$sup_ub,
    max(vapply(
      parts(null@region),
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
  outcomes <- enumerate_space(family@sample_space)
  values <- stats::runif(nrow(outcomes), 0, 10)
  res <- certify(tabulated_rv(family, values), null, tol = 1e-9)

  for (s in parts(null@region)) {
    weights <- matrix(stats::rgamma(3L * 300L, shape = 1), nrow = 3L)
    theta <- s@vertices %*% div_by_col(weights, colSums(weights))
    expectations <- as.vector(
      crossprod(exp(kernel_loglik_batch(family, theta, outcomes)), values)
    )
    expect_lte(max(expectations), res$sup_ub)
  }
})

# --- Return shape -------------------------------------------------------------

test_that("certify() reports one entry per part", {
  set.seed(105)
  null <- plurality_null(n = 6L, k = 4L)
  values <- stats::runif(nrow(enumerate_space(null@family@sample_space)), 0, 10)
  res <- certify(tabulated_rv(null@family, values), null, tol = 1e-6)

  n_sub <- length(parts(null@region))
  expect_length(res$bounds, n_sub)
  expect_length(res$incumbents, n_sub)
  expect_length(res$iterations, n_sub)
  expect_length(res$converged, n_sub)
  expect_length(res$budget_hit, n_sub)
  expect_type(res$iterations, "integer")
  expect_type(res$converged, "logical")
  expect_type(res$budget_hit, "logical")
  expect_identical(res$method, "bernstein")
  expect_equal(res$sup_ub, max(res$bounds))
  expect_equal(res$sup_lb, max(res$incumbents))
})

test_that("certify() carries back the variable and null it holds for", {
  # A certificate that does not name what it certifies is not a certificate.
  null <- plurality_null(n = 4L, k = 3L)
  x <- tabulated_rv(
    null@family,
    rep(1, nrow(enumerate_space(null@family@sample_space)))
  )
  res <- certify(x, null)
  expect_identical(res$random_variable, x)
  expect_identical(res$null, null)
})

test_that("converged and budget_hit distinguish the two ways of stopping", {
  set.seed(111)
  null <- plurality_null(n = 8L, k = 3L)
  values <- stats::runif(nrow(enumerate_space(null@family@sample_space)), 0, 10)
  x <- tabulated_rv(null@family, values)

  tight <- certify(x, null, tol = 1e-12, max_nodes = 5000L)
  expect_true(all(tight$converged))
  expect_false(any(tight$budget_hit))

  starved <- certify(x, null, tol = 1e-12, max_nodes = 2L)
  expect_true(any(starved$budget_hit))
  # Mutually exclusive per part: a search stops one way or the other.
  expect_false(any(starved$converged & starved$budget_hit))
  # Every part that ran out of budget used all of it.
  expect_true(all(starved$iterations[starved$budget_hit] == 2L))

  # The starved bound is still valid, just looser -- which is the whole reason
  # the distinction is worth reporting rather than erroring on.
  expect_gte(starved$sup_ub, tight$sup_ub)
})

# --- Recording ----------------------------------------------------------------

test_that("certify_trace() records every node, and they tile at every step", {
  set.seed(112)
  null <- plurality_null(n = 8L, k = 3L)
  x <- tabulated_rv(
    null@family,
    stats::runif(nrow(enumerate_space(null@family@sample_space)), 0, 10)
  )
  nodes <- certify_trace(x, null, tol = 1e-9)

  expect_s3_class(nodes, "data.frame")
  # One tree per cell, and the reported per-part `iterations` is their total.
  iterations <- attr(nodes, "certificate")$iterations
  for (i in seq_along(null@cells)) {
    rows <- nodes[nodes$cell == i, ]
    it <- sum(rows$fate == "split")
    # `it` splits create two children each, on top of the seed.
    expect_identical(nrow(rows), 1L + 2L * it)
    expect_identical(sum(rows$fate != "split"), it + 1L)
    expect_identical(unique(rows$part), null@cell_part[[i]])

    seed_volume <- abs(det(null@cells[[i]]@vertices))
    for (t in unique(c(0L, seq_len(it)))) {
      drawn <- rows[
        rows$born <= t & !(rows$fate == "split" & rows$retired <= t),
      ]
      expect_equal(sum(drawn$volume), seed_volume)
    }
  }
  for (p in seq_along(parts(null@region))) {
    rows <- nodes[nodes$part == p, ]
    expect_identical(sum(rows$fate == "split"), iterations[[p]])
  }
})

test_that("certify_trace() records the tree and the order it was built in", {
  set.seed(113)
  null <- plurality_null(n = 8L, k = 3L)
  x <- tabulated_rv(
    null@family,
    stats::runif(nrow(enumerate_space(null@family@sample_space)), 0, 10)
  )
  nodes <- certify_trace(x, null, tol = 1e-9)

  for (i in unique(nodes$cell)) {
    rows <- nodes[nodes$cell == i, ]
    # Ids restart at 1 in every cell, so a cell is the scope to read them in.
    it <- sum(rows$fate == "split")

    # Every id issued appears exactly once, so no node goes unrecorded.
    expect_identical(anyDuplicated(rows$id), 0L)
    expect_identical(sort(rows$id), seq_len(1L + 2L * it))

    # One seed, and it is the only node without a parent.
    expect_identical(sum(is.na(rows$parent)), 1L)
    expect_identical(rows$depth[is.na(rows$parent)], 0L)
    expect_identical(rows$born[is.na(rows$parent)], 0L)

    # A node leaves the active set no earlier than it entered, and only a node
    # still live at the end has no retirement.
    expect_true(all(is.na(rows$retired) == (rows$fate == "active")))
    done <- !is.na(rows$retired)
    expect_true(all(rows$retired[done] >= rows$born[done]))
    expect_true(all(rows$born <= it))

    # A child is one level below its parent and born when the parent retired.
    parents <- rows[match(rows$parent, rows$id), ]
    has_parent <- !is.na(rows$parent)
    expect_equal(rows$depth[has_parent], parents$depth[has_parent] + 1L)
    expect_equal(rows$born[has_parent], parents$retired[has_parent])
    expect_true(all(parents$fate[has_parent] == "split"))
  }
})

test_that("certify_trace() agrees with certify() and drops the coefficients", {
  # Same computation, different report. And the record must be cheap: holding
  # `coef` would make it the size of the run rather than the size of the tree.
  set.seed(114)
  null <- plurality_null(n = 8L, k = 3L)
  x <- tabulated_rv(
    null@family,
    stats::runif(nrow(enumerate_space(null@family@sample_space)), 0, 10)
  )

  nodes <- certify_trace(x, null, tol = 1e-9)
  direct <- certify(x, null, tol = 1e-9)
  expect_equal(attr(nodes, "certificate")$sup_ub, direct$sup_ub)
  expect_equal(attr(nodes, "certificate")$sup_lb, direct$sup_lb)

  expect_false("coef" %in% names(nodes))
  expect_true(all(vapply(nodes$vertices, is.matrix, logical(1L))))

  bounds <- attr(nodes, "certificate")$bounds
  for (i in unique(nodes$cell)) {
    rows <- nodes[nodes$cell == i, ]
    # A part's bound is the largest of its cells', so a cell's own nodes sit
    # under it whichever cell of the part they came from.
    expect_lte(
      max(rows$upper[rows$fate != "split"]),
      bounds[[unique(rows$part)]] + 1e-9
    )

    parents <- rows[match(rows$parent, rows$id), ]
    has_parent <- !is.na(rows$parent)
    expect_true(all(
      rows$upper[has_parent] <= parents$upper[has_parent] + 1e-12
    ))
  }
})

test_that("certify() does not record unless asked", {
  null <- plurality_null(n = 6L, k = 3L)
  x <- tabulated_rv(
    null@family,
    rep(1, nrow(enumerate_space(null@family@sample_space)))
  )
  expect_null(certify(x, null)$record)
  expect_null(certify(x, null)$traces)
  expect_type(certify(x, null)$iterations, "integer")
  expect_type(certify(x, null, .record = TRUE)$iterations, "integer")
})

# --- Registry dispatch --------------------------------------------------------

test_that("certify() sends each part to the bound_fn that claimed it", {
  # The reason the registry exists. With one entry the grouping is trivially
  # the identity, so this stubs a second entry to check that results are not
  # merely produced but land in the right slots.
  seen <- list()
  fake_bound_fn <- function(x, family, cells, control) {
    seen[[length(seen) + 1L]] <<- vapply(cells, class_name, character(1L))
    lapply(seq_along(cells), function(i) {
      list(
        bound = 99,
        incumbent = 0,
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
          subject = "The real one",
          fit = function(cell, family) if (bernstein_compatible(cell)) TRUE,
          bound_fn = bernstein_bound
        ),
        list(
          name = "fake",
          subject = "The stub",
          fit = function(cell, family) {
            if (S7_inherits(cell, halfspace_region)) TRUE
          },
          bound_fn = fake_bound_fn
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
      halfspace_region(normal = c(1, -1, 0), offset = 0),
      simplex_region(vertices = facet),
      halfspace_region(normal = c(0, 1, -1), offset = 0)
    )
  )
  x <- tabulated_rv(family, rep(1, nrow(enumerate_space(family@sample_space))))
  res <- certify(x, null, tol = 1e-9)

  # Parts 1 and 3 went to the stub, part 2 to the real bound_fn
  expect_equal(res$bounds[c(1L, 3L)], c(99, 99))
  expect_lt(res$bounds[[2L]], 99)
  expect_setequal(res$method, c("bernstein", "fake"))
  # The stub was called once, with both of its parts together.
  expect_length(seen, 1L)
  expect_identical(seen[[1L]], c("halfspace_region", "halfspace_region"))
})

test_that("certify() rejects a bound_fn that returns the wrong number of results", {
  local_mocked_bindings(
    certify_methods = function() {
      list(list(
        name = "short",
        subject = "The stub",
        fit = function(cell, family) if (bernstein_compatible(cell)) TRUE,
        bound_fn = function(x, family, cells, control) list()
      ))
    }
  )
  null <- plurality_null(n = 4L, k = 3L)
  x <- tabulated_rv(
    null@family,
    rep(1, nrow(enumerate_space(null@family@sample_space)))
  )
  expect_error(certify(x, null), "returned 0 results for 2 cells")
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
    list(halfspace_region(normal = c(1, -1), offset = 0))
  )
  x <- random_variable(
    function(x) rep(1, nrow(as.matrix(x))),
    sample_space = family@sample_space
  )
  expect_error(certify(x, null), "No bounding method is implemented")
  expect_error(certify(x, null), "gaussian_family")
  expect_error(certify(x, null), "halfspace_region")
})

test_that("lower-dimensional nulls fit but do not certify", {
  fam <- multinomial_family(n_trials = 6L, k = 3L)
  # The tie null {theta_1 == theta_2} within the simplex is a segment, and
  # a legit simplex of dimension 1. It is a valid null but currently no
  # certify method is  implemented so we expect it to fail at certification.
  tie <- simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 0, 1)))
  null <- null_model(fam, list(tie))
  Q <- fam(c(0.6, 0.2, 0.2))

  # Fitting works. chart() gives one coordinate, project() and maximise_over()
  # are indifferent to the cell's dimension, and the KL objective is defined.
  set.seed(1)
  state <- ripr_init(Q, null)
  state <- fw_step(
    state,
    times = 20L,
    record_gap = TRUE,
    until = gap_below(1e-10)
  )
  fit <- ripr_finish(state, reoptimise = TRUE, identify = TRUE)
  expect_true(is.finite(fit$kl))
  # Every atom landed on the tie, which is the point: the geometry is honoured.
  expect_equal(atoms(fit$W0)[1L, ], atoms(fit$W0)[2L, ])

  X <- likelihood(Q) / likelihood(fit$P_star)
  expect_true(is.finite(sup_lb(X, null)$sup_lb))

  # Certifying does now.
  # The reason lower-dimensional regions are out of scope is that three future
  # issues will show up:
  #
  #   1. Reparametrisation requires full-dimension targets. `reparametrise_to()`
  #      can probably be fixed to drop dimensions, but that probably is not
  #      going to work without substantial revisions.
  #
  #   2. Every Lebesgue integrals over it will be zero, so truncated mixings
  #      returns -Inf for all outcomes and a rejection sampler will never
  #      terminate. This is because `contains()` on a lower-dimensional cell
  #      is a measure-zero test. So it is entirely tolerance-governed in a way
  #      full-dimensional cells are not. No sampled point ever lands exactly on a
  #      segment.
  #
  #   3. Also we need to decide whether we tolerate taking
  #      `complement(cell, within = simplex)`, which is technically the entire
  #      space when `cell` has lower dimension.
  expect_error(certify(X, null))
})

test_that("the refusal names the failing condition, not the class", {
  # A `simplex_region` the Bernstein enclosure cannot take is not a missing
  # method. The underlying geometry is supported but region is the wrong shape.
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  x <- tabulated_rv(fam, rep(1, nrow(enumerate_space(fam@sample_space))))

  flat <- simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 0, 1)))
  msg <- tryCatch(
    certify(x, null_model(fam, list(flat))),
    error = conditionMessage
  )
  expect_match(msg, "2 vertices in 3 dimensions")
  expect_match(msg, "simplex of dimension 1")
  expect_match(msg, "Only certification is affected")
  expect_false(grepl("No bounding method is implemented", msg))

  # A simplex outside the standard simplex fails on membership, and reports
  # that rather than the vertex count. A tetrahedron in R^3 is affinely
  # independent and not lower-dimensional at all, so a count-based message
  # would be actively wrong.
  tetra <- simplex_region(
    vertices = cbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
  )
  msg <- tryCatch(
    certify(x, null_model(fam, list(tetra))),
    error = conditionMessage
  )
  expect_match(msg, "leave the standard simplex")
  expect_match(msg, "sum to 0 rather than 1")
  expect_false(grepl("simplex of dimension", msg))

  # A negative coordinate is caught before the sum, and named.
  outside <- simplex_region(
    vertices = cbind(c(-0.5, 1.5, 0), c(0, 1, 0), c(0, 0, 1))
  )
  msg <- tryCatch(
    certify(x, null_model(fam, list(outside))),
    error = conditionMessage
  )
  expect_match(msg, "smallest coordinate is -0.5")

  # An unbounded region fails for a reason of its own, and says which. It is
  # not a missing method either: nothing is missing, there is simply no finite
  # simplicial cover to enclose.
  msg <- tryCatch(
    certify(
      x,
      null_model(fam, list(halfspace_region(normal = c(1, -1, 0), offset = 0)))
    ),
    error = conditionMessage
  )
  expect_match(msg, "it is unbounded")
  expect_match(msg, "State the null over a bounded region")
  expect_false(grepl("No bounding method is implemented", msg))

  # Whereas a bounded hull is triangulated rather than refused: the square
  # `{theta_1 <= 1/2, theta_2 <= 1/2}` inside the 2-simplex is two triangles,
  # and each of those the enclosure takes.
  expect_no_error(
    certify(x, null_model(fam, list(simplex_square())))
  )
})

test_that("a region obstruction is not blamed if the family is not implemented", {
  # `bernstein_obstruction()` describes a region, so it must only speak for a
  # family some method actually claims. A `gaussian_family` over a flat
  # `simplex_region` fails because nothing bounds Gaussian expectations at all.
  # A full-dimensional region would fail identically, so raising an error
  # mentioning the region's shape does not give adequate advice.
  fam <- gaussian_family(dim = 2L)
  flat <- simplex_region(vertices = cbind(c(1, 0)))
  x <- random_variable(
    function(x) rep(1, nrow(as.matrix(x))),
    sample_space = fam@sample_space
  )
  msg <- tryCatch(
    certify(x, null_model(fam, list(flat))),
    error = conditionMessage
  )
  expect_match(msg, "No bounding method is implemented")
  expect_match(msg, "gaussian_family")
  expect_false(grepl("Bernstein", msg))
  expect_false(grepl("standard simplex", msg))
})

test_that("bernstein_compatible() accepts when certification is possible", {
  # Full parameter space
  expect_true(bernstein_compatible(simplex_region(vertices = diag(3))))
  # Pairwise plurality
  expect_true(bernstein_compatible(
    simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1)))
  ))

  # A `polytope_region` is refused even when its hull is a simplex, because
  # nothing has established that it is one. That is the predicate's job, not a
  # limitation: `certify()` sends it `cells()`, and the fan hands back
  # `simplex_region`s whatever the part was declared as.
  expect_false(bernstein_compatible(polytope_region(vertices = diag(3))))
  expect_true(all(vapply(
    cells(polytope_region(vertices = diag(3))),
    bernstein_compatible,
    logical(1)
  )))

  expect_false(bernstein_compatible(
    simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 0, 1)))
  ))
  expect_false(bernstein_compatible(simplex_region(
    vertices = cbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
  )))
  expect_false(bernstein_compatible(halfspace_region(
    normal = c(1, -1, 0),
    offset = 0
  )))
  expect_false(bernstein_compatible(unconstrained_region(3L)))

  # The degenerate single vertex has no edges to condition
  expect_true(bernstein_compatible(simplex_region(vertices = matrix(1))))
})


test_that("the refusal names the part, not the null model", {
  # The message is meant to say which geometry is missing a bound. Naming the
  # container instead makes it useless.
  family <- gaussian_family(dim = 2L)
  null <- null_model(
    family,
    list(halfspace_region(normal = c(1, -1), offset = 0))
  )
  x <- random_variable(
    function(x) rep(1, nrow(as.matrix(x))),
    sample_space = family@sample_space
  )
  msg <- tryCatch(certify(x, null), error = conditionMessage)
  expect_false(grepl("null_model", msg, fixed = TRUE))
  expect_false(grepl("FALSE", msg, fixed = TRUE))
})

test_that("the refusal is not repeated once per part", {
  family <- gaussian_family(dim = 2L)
  null <- null_model(
    family,
    list(
      halfspace_region(normal = c(1, -1), offset = 0),
      halfspace_region(normal = c(1, 0), offset = 0)
    )
  )
  x <- random_variable(
    function(x) rep(1, nrow(as.matrix(x))),
    sample_space = family@sample_space
  )
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
  values <- rep(1, nrow(enumerate_space(null@family@sample_space)))
  values[3L] <- Inf
  expect_error(
    certify(tabulated_rv(null@family, values), null),
    "not finite everywhere"
  )
})

test_that("certify() refuses a lattice above the coefficient budget", {
  null <- plurality_null(n = 8L, k = 3L)
  x <- tabulated_rv(
    null@family,
    rep(1, nrow(enumerate_space(null@family@sample_space)))
  )
  expect_error(
    certify(x, null, max_coefficients = 10L),
    "above `max_coefficients`"
  )
  # The refusal must be raised before any of the work is done.
  expect_error(certify(x, null, max_coefficients = 10L), "would fit")
})

test_that("certify() rejects a non-random_variable and bad control values", {
  null <- plurality_null(n = 4L, k = 3L)
  x <- tabulated_rv(
    null@family,
    rep(1, nrow(enumerate_space(null@family@sample_space)))
  )
  expect_error(certify(function(x) 1, null), "must be a `random_variable`")
  expect_error(certify(x, null, tol = -1))
  expect_error(certify(x, null, max_nodes = 0))
  expect_error(certify(x, null, max_coefficients = 0))
})

test_that("an ill-conditioned simplex is refused by the enclosure by name", {
  # The validator's exact independence test lets slivers through; the
  # conditioning heuristic that used to refuse them at construction now lives
  # here, where it can say what is actually wrong.
  sliver <- cbind(c(1, 0, 0), c(0, 1, 0), c(0.5, 0.5 - 1e-12, 1e-12))
  s <- simplex_region(vertices = sliver)
  expect_false(bernstein_compatible(s))
  expect_match(bernstein_obstruction(s)$because, "ill-conditioned")
  expect_match(bernstein_obstruction(s)$remedy, "numerical limit")
  expect_false(grepl(
    "one vertex per coordinate",
    bernstein_obstruction(s)$remedy
  ))

  ok <- simplex_region(
    vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))
  )
  expect_true(bernstein_compatible(ok))
  expect_null(bernstein_obstruction(ok))
})


test_that("an obstruction gets a follow-up that belongs to it", {
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  x <- tabulated_rv(fam, rep(1, nrow(enumerate_space(fam@sample_space))))
  refusal <- function(region) {
    tryCatch(certify(x, null_model(fam, region)), error = conditionMessage)
  }
  # Lower-dimensional
  msg <- refusal(simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 0, 1))))
  expect_match(msg, "simplex of dimension")
  expect_match(msg, "a parametrisation is not yet implemented")
  expect_false(grepl("No bounding method is implemented", msg))
})


test_that("every bounded region triangulates, whatever class states it", {
  # `cells()` keys on boundedness, not on class, so the claim the unbounded
  # refusal makes.
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  x <- tabulated_rv(fam, rep(1, nrow(enumerate_space(fam@sample_space))))
  bare <- polyhedron_region(
    vertices = cbind(c(0.5, 0.5, 0), c(0, 0.5, 0.5), c(0, 0, 1), c(0.5, 0, 0.5))
  )
  expect_false(S7_inherits(bare, polytope_region))
  expect_length(cells(bare), 2L)
  expect_no_error(certify(x, null_model(fam, bare), tol = 1e-9))

  # And an unbounded one is still its own only cell.
  expect_length(cells(halfspace_region(normal = c(1, -1, 0))), 1L)
})


# --- Point nulls --------------------------------------------------------------

test_that("a point null is certified by evaluation, exactly", {
  # The simple null. `sup` over `{theta}` is `E_theta[X]`, so there is nothing
  # to enclose and nothing to subdivide.
  set.seed(9)
  family <- multinomial_family(n_trials = 6L, k = 3L)
  outcomes <- enumerate_space(family@sample_space)
  values <- stats::runif(nrow(outcomes), 0, 10)
  x <- tabulated_rv(family, values)
  theta <- c(0.5, 0.3, 0.2)

  res <- certify(x, null_model(family, point_region(theta = theta)), tol = 0)
  direct <- sum(
    exp(as.vector(compile_loglik(family, outcomes)(matrix(theta, ncol = 1L)))) *
      values
  )

  expect_identical(res$method, "point")
  expect_identical(res$iterations, 0L)
  expect_true(all(res$converged))
  expect_false(any(res$budget_hit))

  # The attained value is the evaluation itself, to the last bit.
  expect_identical(res$sup_lb, direct)
  # A point's certificate is its value
  expect_identical(res$sup_ub, res$sup_lb)
})


test_that("the point method takes any family whose sample space is enumerable", {
  # Nothing here is multinomial-specific, so the registry entry claims
  # `parametric_family` and gates on the sample space instead.
  binomial <- multinomial_family(n_trials = 10L, k = 2L)
  x <- tabulated_rv(
    binomial,
    rep(1, nrow(enumerate_space(binomial@sample_space)))
  )
  res <- certify(x, null_model(binomial, point_region(theta = c(0.5, 0.5))))
  # The expectation of the constant 1 is 1, whatever the parameter.
  expect_equal(res$sup_lb, 1)
  expect_gte(res$sup_ub, 1)

  # A continuous sample space has an integral rather than a sum, and quadrature
  # returns an estimate, which a certificate may not rest on.
  gaussian <- gaussian_family(dim = 2L)
  y <- random_variable(
    function(z) rep(1, nrow(as.matrix(z))),
    sample_space = gaussian@sample_space
  )
  msg <- tryCatch(
    certify(y, null_model(gaussian, point_region(theta = c(0, 0)))),
    error = conditionMessage
  )
  expect_match(msg, "Exact evaluation at a point")
  expect_match(msg, "an integral over a `real_space`")
  expect_false(grepl("No bounding method is implemented", msg))
})


test_that("the wildcard does not make every refusal its business", {
  # The point entry claims `parametric_family`, so every family now has a
  # method that "claims" it. That must not let one method's obstruction speak
  # for all others.
  gaussian <- gaussian_family(dim = 2L)
  y <- random_variable(
    function(z) rep(1, nrow(as.matrix(z))),
    sample_space = gaussian@sample_space
  )
  msg <- tryCatch(
    certify(y, null_model(gaussian, halfspace_region(normal = c(1, -1)))),
    error = conditionMessage
  )
  expect_match(msg, "No bounding method is implemented")
  expect_false(grepl("standard simplex", msg))
  expect_false(grepl("Exact evaluation", msg))
})


test_that("a point outside the parameter space is refused, not evaluated", {
  # An expectation at a point the family has no distribution for is not an
  # expectation under the null, so the method declines it and the family's own
  # method explains why.
  family <- multinomial_family(n_trials = 4L, k = 3L)
  x <- tabulated_rv(family, rep(1, nrow(enumerate_space(family@sample_space))))
  outside <- point_region(theta = c(2, -1, 0))
  expect_null(point_fit(outside, family))

  msg <- tryCatch(
    certify(x, null_model(family, outside)),
    error = conditionMessage
  )
  expect_match(msg, "leave the standard simplex")
  expect_match(msg, "basis polynomials may take negative values")
})


test_that("a null mixing a point with a simplex certifies under both methods", {
  # The reason the registry groups by method rather than refusing a null whose
  # cells differ: each part goes to what can take it, and one certificate comes
  # back.
  set.seed(10)
  family <- multinomial_family(n_trials = 6L, k = 3L)
  outcomes <- enumerate_space(family@sample_space)
  x <- tabulated_rv(family, stats::runif(nrow(outcomes), 0, 10))
  null <- null_model(
    family,
    list(
      point_region(theta = c(0.5, 0.3, 0.2)),
      simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1)))
    )
  )

  res <- certify(x, null, tol = 1e-9)
  expect_setequal(res$method, c("point", "bernstein"))
  expect_length(res$bounds, 2L)
  expect_identical(res$iterations[[1L]], 0L)
  expect_gt(res$iterations[[2L]], 0L)
  expect_equal(res$sup_ub, max(res$bounds))

  # Each part's bound covers its own part, and the whole covers both.
  alone <- vapply(
    parts(null@region),
    function(p) certify(x, null_model(family, p), tol = 1e-9)$sup_ub,
    numeric(1)
  )
  expect_gte(res$sup_ub, max(alone) - 1e-6)
})


# --- sup_lb -------------------------------------------------------------------

test_that("sup_lb() reports a value the objective actually attains", {
  set.seed(106)
  null <- plurality_null(n = 8L, k = 3L)
  family <- null@family
  outcomes <- enumerate_space(family@sample_space)
  values <- stats::runif(nrow(outcomes), 0, 10)
  x <- tabulated_rv(family, values)
  found <- sup_lb(x, null, n_seeds = 100L, n_restarts = 10L)

  # E_theta[X] for theta supplied as columns.
  expectations <- function(theta) {
    as.vector(crossprod(
      exp(kernel_loglik_batch(family, theta, outcomes)),
      values
    ))
  }

  expect_equal(found$sup_lb, expectations(matrix(found$theta, ncol = 1L)))

  expect_true(contains(parts(null@region)[[found$part]], found$theta))

  vertex_best <- max(vapply(
    parts(null@region),
    function(s) max(expectations(s@vertices)),
    numeric(1L)
  ))
  expect_gt(found$sup_lb, vertex_best)

  expect_gt(found$sup_lb / certify(x, null, tol = 1e-12)$sup_ub, 0.999)
})

test_that("sup_lb() improves on a single seed given more of them", {
  # Obviously there is some randomness in test, but should often be the case for
  # non-convex solves. We fix one such seed here as a sanity check.
  set.seed(109)
  null <- plurality_null(n = 8L, k = 3L)
  outcomes <- enumerate_space(null@family@sample_space)

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
  expect_true(all(ratios > 0.99))
})

test_that("sup_lb() rejects a non-random_variable", {
  null <- plurality_null(n = 4L, k = 3L)
  expect_error(sup_lb(function(x) 1, null), "must be a `random_variable`")
})


# --- Convex nulls that are not simplices --------------------------------------

test_that("a polytope null certifies through its triangulation", {
  # The acceptance case for the whole decomposition: a null stated as a convex
  # hull, with no simplex anywhere in what the caller wrote, certifying end to
  # end.
  set.seed(207)
  family <- multinomial_family(n_trials = 6L, k = 3L)
  null <- null_model(family, list(simplex_square()))
  x <- tabulated_rv(
    family,
    stats::runif(nrow(enumerate_space(family@sample_space)), 0, 10)
  )

  res <- certify(x, null, tol = 1e-9)
  expect_length(res$bounds, 1L)
  expect_identical(res$method, "bernstein")
  expect_true(all(res$converged))

  # One part, two cells, and the certificate is reduced back onto the part.
  expect_length(null@cells, 2L)
  expect_identical(null@cell_part, c(1L, 1L))

  # One-sidedness, against a dense grid over the square. The grid is a lower
  # bound on the true supremum, so this is necessary and not sufficient, which
  # is the same standard the simplex cases here are held to.
  chart <- chart(null@cells[[1L]])
  grid <- do.call(
    cbind,
    lapply(null@cells, function(cell) {
      w <- as.matrix(expand.grid(rep(list(seq(0, 1, length.out = 21L)), 2L)))
      w <- cbind(w, 1 - rowSums(w))
      w <- w[rowSums(w >= 0) == 3L, , drop = FALSE]
      cell@vertices %*% t(w)
    })
  )
  outcomes <- enumerate_space(family@sample_space)
  values <- x(outcomes)
  attained <- as.vector(
    crossprod(exp(kernel_loglik_batch(family, grid, outcomes)), values)
  )
  expect_lte(max(attained), res$sup_ub)
  expect_gte(res$sup_lb, max(attained) - 1e-6)
})


test_that("an incumbent found in one cell prunes the searches over the rest", {
  # The point of sharing it: a cell that cannot beat what is already known
  # stops as soon as its active set falls below that value, instead of proving
  # its own supremum to `tol` for an answer the maximum discards anyway.
  set.seed(3)
  family <- multinomial_family(n_trials = 12L, k = 3L)
  outcomes <- enumerate_space(family@sample_space)
  x <- tabulated_rv(family, stats::runif(nrow(outcomes), 0, 10))
  null <- null_model(family, list(simplex_square()))

  # Each cell certified as a null of its own gets no incumbent from the other,
  # which is the comparison: same cells, same tolerance, no sharing.
  alone <- lapply(
    null@cells,
    function(cell) certify(x, null_model(family, list(cell)), tol = 1e-9)
  )
  shared <- certify(x, null, tol = 1e-9)

  expect_lt(
    shared$iterations,
    sum(vapply(alone, \(r) r$iterations, integer(1)))
  )

  # And the certificate is unaffected. Pruning against a value that was
  # actually attained cannot drop the maximiser, and the bound each pruned run
  # returns still accounts for what it dropped.
  separate_ub <- max(vapply(alone, \(r) r$sup_ub, numeric(1)))
  expect_equal(shared$sup_ub, separate_ub, tolerance = 1e-6)
  expect_equal(
    shared$sup_lb,
    max(vapply(alone, \(r) r$sup_lb, numeric(1)))
  )
})


test_that("a dominated part reports its own bound, not the incumbent's", {
  family <- multinomial_family(n_trials = 10L, k = 3L)
  outcomes <- enumerate_space(family@sample_space)
  x <- tabulated_rv(family, outcomes[, 1L]^2)
  large <- simplex_region(
    vertices = cbind(c(1, 0, 0), c(0.5, 0.5, 0), c(0.5, 0, 0.5))
  )
  small <- simplex_region(
    vertices = cbind(c(0, 0, 1), c(0, 0.5, 0.5), c(0.5, 0, 0.5))
  )
  res <- certify(x, null_model(family, list(large, small)), tol = 1e-9)
  alone <- certify(x, null_model(family, list(small)), tol = 1e-9)

  # The large part's attained value dominates everything the poor part has.
  expect_gt(res$incumbents[1L], alone$sup_ub)
  # The small part still reports a valid bound on itself, far below the
  # incumbent it was pruned against.
  expect_gte(res$bounds[2L], alone$sup_lb)
  expect_lt(res$bounds[2L], res$incumbents[1L])
})


test_that("the incumbent carries from one bound_fn group to the next", {
  # Nothing about the reduction stops a value attained under one method from
  # pruning a search under another: `sup_lb` is a maximum over every cell of
  # every group, so any group's incumbent bounds the supremum from below. With
  # one entry in the registry there is only ever one group, so this stubs a
  # second to check that the value is actually handed on.
  seen <- numeric(0)
  fake_bound_fn <- function(x, family, cells, control) {
    seen <<- c(seen, control$incumbent)
    lapply(seq_along(cells), function(i) {
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
  local_mocked_bindings(
    certify_methods = function() {
      list(
        list(
          name = "fake",
          subject = "The stub",
          fit = function(cell, family) {
            if (S7_inherits(cell, halfspace_region)) TRUE
          },
          bound_fn = fake_bound_fn
        ),
        list(
          name = "bernstein",
          subject = "The real one",
          fit = function(cell, family) if (bernstein_compatible(cell)) TRUE,
          bound_fn = bernstein_bound
        )
      )
    }
  )

  family <- multinomial_family(n_trials = 4L, k = 3L)
  facet <- diag(3L)
  facet[, 1L] <- c(0.5, 0.5, 0)
  x <- tabulated_rv(family, rep(1, nrow(enumerate_space(family@sample_space))))

  # The stub runs first and attains 99, which the real method then has to be
  # handed. Reversing the order shows the first group starts from nothing.
  first <- null_model(
    family,
    list(
      halfspace_region(normal = c(1, -1, 0), offset = 0),
      simplex_region(vertices = facet)
    )
  )
  certify(x, first, tol = 1e-9)
  expect_identical(seen, -Inf)

  seen <- numeric(0)
  second <- null_model(
    family,
    list(
      simplex_region(vertices = facet),
      halfspace_region(normal = c(1, -1, 0), offset = 0)
    )
  )
  certify(x, second, tol = 1e-9)
  # The Bernstein group ran first, so the stub was handed what it attained --
  # a real value, not the `-Inf` a first group starts from.
  expect_length(seen, 1L)
  expect_true(is.finite(seen))
  expect_gt(seen, 0)
})
