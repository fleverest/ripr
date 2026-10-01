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
    vertices[1L, ] <- replace(numeric(k), c(1L, j), 0.5)
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
    vertices = rbind(
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
  k <- nrow(vertices)
  weights <- enumerate_space(count_space(n_trials = m, k = k)) / m
  theta <- weights %*% vertices
  outcomes <- enumerate_space(family@sample_space)
  max(as.vector(crossprod(
    exp(compile_loglik(family, outcomes)(theta)),
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
    # The bound and the grid evaluate the same expectations in different
    # orders, so the two can disagree in rounding.
    expect_gte(
      res@sup_ub + rounding_tol(res@sup_ub),
      null_grid_max(null, values)
    )
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
    function(m) certify(x, null, tol = 0, max_splits = m)@sup_ub,
    numeric(1L)
  )
  expect_true(all(bounds + rounding_tol(truth) >= truth))
  expect_false(is.unsorted(rev(bounds))) # non-increasing in the budget
})

test_that("sup_ub brackets sup_lb", {
  set.seed(103)
  null <- plurality_null(n = 10L, k = 3L)
  values <- stats::runif(nrow(enumerate_space(null@family@sample_space)), 0, 10)
  x <- tabulated_rv(null@family, values)
  res <- certify(x, null, tol = 1e-9)
  expect_gte(res@sup_ub, res@sup_lb)
  # The searched lower bound is a different algorithm on the same problem, so
  # this is the cross-check between the two halves of the file. SLSQP may attain
  # the supremum, so if the two meet one may round higher than the other, hence
  # the floating-point slack.
  expect_gte(
    res@sup_ub + rounding_tol(res@sup_ub),
    sup_lb(x, null, n_seeds = 100L, n_restarts = 10L)@sup_lb
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

  # The e-variable, which `e_variable()` builds by dividing by the bound:
  expect_gt(res@sup_ub, 1)
  e <- e_variable(res)
  e_vals <- e(outcomes)

  expect_equal(e_vals, values / res@sup_ub)

  # E_theta[e] for theta supplied as rows.
  expectations <- function(theta) {
    as.vector(crossprod(
      exp(compile_loglik(family, outcomes)(theta)),
      e_vals
    ))
  }

  for (s in parts(null@region)) {
    weights <- matrix(stats::rgamma(4L * 200L, shape = 1), ncol = 4L)
    theta <- (weights / rowSums(weights)) %*% s@vertices
    expect_lte(max(expectations(theta)), 1 + rounding_tol(1))
  }

  attained <- c(
    vapply(
      parts(null@region),
      function(s) max(expectations(s@vertices)),
      numeric(1L)
    ),
    expectations(matrix(
      sup_lb(x, null, n_seeds = 200L, n_restarts = 20L)@theta,
      nrow = 1L
    ))
  )
  expect_lte(max(attained), 1 + rounding_tol(1))
  expect_gt(max(attained), 1 - 1e-3)
  expect_gt(res@sup_lb / res@sup_ub, 1 - 1e-6)
  expect_lte(res@sup_lb / res@sup_ub, 1)
})

test_that("a constant variable certifies to its own value", {
  null <- plurality_null(n = 6L, k = 3L)
  x <- tabulated_rv(
    null@family,
    rep(3.5, nrow(enumerate_space(null@family@sample_space)))
  )
  res <- certify(x, null, tol = 1e-9)
  expect_equal(res@sup_lb, 3.5)
  expect_lt(res@sup_ub - 3.5, rounding_tol(3.5))
  expect_true(all(res@iterations == 0L))
})

test_that("a variable maximised at a facet vertex needs no subdivision", {
  # The enclosure is exact at the vertices (PBP 10.2), so if the maximiser is
  # a vertex the starting bound is already attained.
  n <- 6L
  null <- plurality_null(n = n, k = 3L)
  family <- null@family
  outcomes <- enumerate_space(family@sample_space)
  q <- c(0.1, 0.8, 0.1)
  values <- exp(as.vector(outcomes %*% (log(q) - log(rep(1 / 3, 3)))))
  res <- certify(tabulated_rv(family, values), null, tol = 1e-9)

  expect_true(all(res@iterations == 0L))
  expect_equal(res@sup_ub, res@sup_lb, tolerance = 1e-9)

  expect_equal(res@sup_ub, (3 * max(q))^n, tolerance = 1e-9)
  expect_equal(
    res@bounds,
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

  expect_true(all(res@iterations > 0L))
  # And the work bought something: the bound is below every seed node's, which
  # is the enclosure over a whole facet before any subdivision.
  nodes <- certify_trace(tabulated_rv(null@family, values), null, tol = 1e-9)
  seeds <- nodes[is.na(nodes$parent), ]
  expect_identical(nrow(seeds), length(null@cells))
  expect_lt(res@sup_ub, max(seeds$upper))
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
    weights <- matrix(stats::rgamma(3L * 300L, shape = 1), ncol = 3L)
    theta <- (weights / rowSums(weights)) %*% s@vertices
    expectations <- as.vector(
      crossprod(exp(compile_loglik(family, outcomes)(theta)), values)
    )
    expect_lte(max(expectations), res@sup_ub)
  }
})

# --- Return shape -------------------------------------------------------------

test_that("converged and budget_hit distinguish the two ways of stopping", {
  set.seed(111)
  null <- plurality_null(n = 8L, k = 3L)
  values <- stats::runif(nrow(enumerate_space(null@family@sample_space)), 0, 10)
  x <- tabulated_rv(null@family, values)

  full <- certify(x, null, tol = 1e-12, max_splits = 5000L)
  expect_true(all(full@converged))
  expect_false(any(full@budget_hit))

  starved <- certify(x, null, tol = 1e-12, max_splits = 2L)
  expect_true(any(starved@budget_hit))
  # Mutually exclusive per part: a search stops one way or the other.
  expect_false(any(starved@converged & starved@budget_hit))
  # Every part that ran out of budget used all of it.
  expect_true(all(starved@iterations[starved@budget_hit] == 2L))

  # The starved bound is still valid, just looser -- which is the whole reason
  # the distinction is worth reporting rather than erroring on.
  expect_gte(starved@sup_ub, full@sup_ub)
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
  iterations <- attr(nodes, "certificate")@iterations
  for (i in seq_along(null@cells)) {
    rows <- nodes[nodes$cell == i, ]
    it <- sum(rows$fate == "split")
    # `it` splits create two children each, on top of the seed.
    expect_identical(nrow(rows), 1L + 2L * it)
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
  expect_true(S7::S7_inherits(attr(nodes, "certificate"), ripr_certificate))
  expect_equal(attr(nodes, "certificate")@sup_ub, direct@sup_ub)
  expect_equal(attr(nodes, "certificate")@sup_lb, direct@sup_lb)

  expect_false("coef" %in% names(nodes))
  expect_true(all(vapply(nodes$vertices, is.matrix, logical(1L))))

  bounds <- attr(nodes, "certificate")@bounds
  for (i in unique(nodes$cell)) {
    rows <- nodes[nodes$cell == i, ]
    # A part's bound is the largest of its cells', so a cell's own nodes sit
    # under it whichever cell of the part they came from.
    expect_lte(
      max(rows$upper[rows$fate != "split"]),
      bounds[[unique(rows$part)]]
    )

    parents <- rows[match(rows$parent, rows$id), ]
    has_parent <- !is.na(rows$parent)
    expect_true(all(
      rows$upper[has_parent] <=
        parents$upper[has_parent] + rounding_tol(parents$upper[has_parent])
    ))
  }
})

# --- Refusals -----------------------------------------------------------------

test_that("certify() refuses a family and geometry it has no method for", {
  family <- gaussian_family(d = 2L)
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

test_that("lower-dimensional nulls fit and certify", {
  fam <- multinomial_family(n_trials = 6L, k = 3L)
  # The tie null {theta_1 == theta_2} within the simplex is a segment, a
  # simplex of dimension 1.
  tie <- simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 0, 1)))
  null <- null_model(fam, list(tie))
  Q <- fam(c(0.6, 0.2, 0.2))

  # chart() gives one coordinate, project() and the oracle's search are
  # indifferent to the cell's dimension, and the KL objective is defined.
  set.seed(1)
  state <- ripr_init(Q, null)
  state <- fw_step(state, times = 20L, until = gap_below(1e-10))
  fit <- ripr_finish(state, reoptimise = TRUE, identify = TRUE)
  expect_true(is.finite(fit@kl))
  # Every atom landed on the tie, which is the point: the geometry is honoured.
  expect_equal(atoms(fit@W0)[, 1L], atoms(fit@W0)[, 2L])

  X <- likelihood(Q) / likelihood(fit@P_star)
  cert <- certify(X, null, tol = 1e-9)
  values <- evaluate_on_space(X, fam)
  outcomes <- as.matrix(enumerate_space(fam@sample_space))
  along <- vapply(seq(0, 1, length.out = 2001L), function(s) {
    theta <- s * c(0.5, 0.5, 0) + (1 - s) * c(0, 0, 1)
    sum(values * apply(outcomes, 1L, stats::dmultinom, prob = theta))
  }, numeric(1))
  expect_gte(cert@sup_ub, max(along) - 1e-9)
  expect_lte(cert@sup_ub, max(along) + 1e-6)
})

test_that("the refusal names the failing condition, not the class", {
  # A `simplex_region` the Bernstein enclosure cannot take is not a missing
  # method. The underlying geometry is supported but region is the wrong shape.
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  x <- tabulated_rv(fam, rep(1, nrow(enumerate_space(fam@sample_space))))

  # A simplex outside the standard simplex fails on membership.
  tetra <- simplex_region(
    vertices = rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
  )
  msg <- tryCatch(
    certify(x, null_model(fam, list(tetra))),
    error = conditionMessage
  )
  expect_match(msg, "leave the standard simplex")
  expect_match(msg, "sum to 0 rather than 1")
  expect_false(grepl("No bounding method is implemented", msg))

  # A negative coordinate is caught before the sum, and named.
  outside <- simplex_region(
    vertices = rbind(c(-0.5, 1.5, 0), c(0, 1, 0), c(0, 0, 1))
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
})

test_that("a region obstruction is not blamed if the family is not implemented", {
  # The Bernstein obstruction describes a region, so it must only speak for a
  # family the enclosure actually claims. A `gaussian_family` over a flat
  # `simplex_region` fails because nothing bounds Gaussian expectations at all.
  # A full-dimensional region would fail identically, so raising an error
  # mentioning the region's shape does not give adequate advice.
  fam <- gaussian_family(d = 2L)
  flat <- simplex_region(vertices = rbind(c(1, 0)))
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

test_that("certification goes ahead exactly where the enclosure applies", {
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  x <- tabulated_rv(fam, rep(1, nrow(enumerate_space(fam@sample_space))))
  certifies <- function(region) {
    !inherits(
      try(certify(x, null_model(fam, list(region)), tol = 1e-9), silent = TRUE),
      "try-error"
    )
  }
  # Full parameter space
  expect_true(certifies(simplex_region(vertices = diag(3))))
  # Pairwise plurality
  expect_true(certifies(
    simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1)))
  ))
  # A `polytope_region` whose hull is a simplex certifies too: `certify()`
  # sends it `cells()`, and the fan hands back `simplex_region`s whatever the
  # part was declared as.
  expect_true(certifies(polytope_region(vertices = diag(3))))
  expect_true(all(vapply(
    cells(polytope_region(vertices = diag(3))),
    \(cell) S7_inherits(cell, simplex_region),
    logical(1)
  )))

  # Lower-dimensional
  expect_true(certifies(
    simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 0, 1)))
  ))
  expect_false(certifies(simplex_region(
    vertices = rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
  )))
  expect_false(certifies(halfspace_region(
    normal = c(1, -1, 0),
    offset = 0
  )))
  expect_false(certifies(real_region(3L)))
})


test_that("the refusal names the part, not the null model", {
  # The message is meant to say which geometry is missing a bound. Naming the
  # container instead makes it useless.
  family <- gaussian_family(d = 2L)
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
  family <- gaussian_family(d = 2L)
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
  expect_error(
    certify(x, null, max_coefficients = 10L),
    "reduce `n_trials`"
  )
})

test_that("certify() rejects a non-random_variable and bad control values", {
  null <- plurality_null(n = 4L, k = 3L)
  x <- tabulated_rv(
    null@family,
    rep(1, nrow(enumerate_space(null@family@sample_space)))
  )
  expect_error(certify(function(x) 1, null), "must be a `random_variable`")
  expect_error(certify(x, null, tol = -1))
  expect_error(certify(x, null, max_splits = 0))
  expect_error(certify(x, null, max_coefficients = 0))
})

test_that("an ill-conditioned simplex certifies", {
  sliver <- rbind(c(1, 0, 0), c(0, 1, 0), c(0.5, 0.5 - 1e-12, 1e-12))
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  x <- tabulated_rv(fam, rep(1, nrow(enumerate_space(fam@sample_space))))
  cert <- certify(x, null_model(fam, simplex_region(vertices = sliver)))
  expect_equal(cert@sup_ub, 1, tolerance = 1e-12)
})


test_that("every bounded region triangulates, whatever class states it", {
  # `cells()` keys on boundedness, not on class, so the claim the unbounded
  # refusal makes.
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  x <- tabulated_rv(fam, rep(1, nrow(enumerate_space(fam@sample_space))))
  bare <- polyhedron_region(
    vertices = rbind(c(0.5, 0.5, 0), c(0, 0.5, 0.5), c(0, 0, 1), c(0.5, 0, 0.5))
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
    exp(as.vector(compile_loglik(family, outcomes)(matrix(theta, nrow = 1L)))) *
      values
  )

  expect_identical(res@method, "point")
  expect_identical(res@iterations, 0L)
  expect_true(all(res@converged))
  expect_false(any(res@budget_hit))

  # The attained value is the evaluation itself, to the last bit.
  expect_identical(res@sup_lb, direct)
  # A point's certificate is its value
  expect_identical(res@sup_ub, res@sup_lb)
})


test_that("the point method takes any family whose sample space is enumerable", {
  # Nothing here is multinomial-specific, so the point method takes any family
  # and gates on the sample space instead.
  binomial <- multinomial_family(n_trials = 10L, k = 2L)
  x <- tabulated_rv(
    binomial,
    rep(1, nrow(enumerate_space(binomial@sample_space)))
  )
  res <- certify(x, null_model(binomial, point_region(theta = c(0.5, 0.5))))
  # The expectation of the constant 1 is 1, whatever the parameter.
  expect_equal(res@sup_lb, 1)
  expect_gte(res@sup_ub, 1)

  # A continuous sample space has an integral rather than a sum, and quadrature
  # returns an estimate, which a certificate may not rest on.
  gaussian <- gaussian_family(d = 2L)
  y <- random_variable(
    function(z) rep(1, nrow(as.matrix(z))),
    sample_space = gaussian@sample_space
  )
  msg <- tryCatch(
    certify(y, null_model(gaussian, point_region(theta = c(0, 0)))),
    error = conditionMessage
  )
  expect_match(msg, "Evaluation at a point")
  expect_match(msg, "an integral over a `real_region`")
  expect_false(grepl("No bounding method is implemented", msg))
})


test_that("the wildcard does not make every refusal its business", {
  # The point entry claims `parametric_family`, so every family has a
  # method that "claims" it. That must not let one method's obstruction speak
  # for all others.
  gaussian <- gaussian_family(d = 2L)
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

  msg <- tryCatch(
    certify(x, null_model(family, outside)),
    error = conditionMessage
  )
  expect_match(msg, "leave the standard simplex")
  expect_match(msg, "basis polynomials may take negative values")
})


test_that("a null mixing a point with a simplex certifies under both methods", {
  # A null whose cells differ is split across methods rather than refused:
  # each part goes to what can take it, and one certificate comes back.
  set.seed(10)
  family <- multinomial_family(n_trials = 6L, k = 3L)
  outcomes <- enumerate_space(family@sample_space)
  x <- tabulated_rv(family, stats::runif(nrow(outcomes), 0, 10))
  null <- null_model(
    family,
    list(
      point_region(theta = c(0.5, 0.3, 0.2)),
      simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1)))
    )
  )

  res <- certify(x, null, tol = 1e-9)
  expect_setequal(res@method, c("point", "bernstein"))
  expect_length(res@bounds, 2L)
  expect_identical(res@iterations[[1L]], 0L)
  expect_gt(res@iterations[[2L]], 0L)
  expect_equal(res@sup_ub, max(res@bounds))

  # Each part's bound covers its own part, and the whole covers both.
  alone <- vapply(
    parts(null@region),
    function(p) certify(x, null_model(family, p), tol = 1e-9)@sup_ub,
    numeric(1)
  )
  expect_gte(res@sup_ub, max(alone) - 1e-6)
})


# --- sup_lb -------------------------------------------------------------------

test_that("sup_lb() reports a value the objective actually attains", {
  set.seed(106)
  null <- plurality_null(n = 8L, k = 3L)
  family <- null@family
  outcomes <- enumerate_space(family@sample_space)
  values <- stats::runif(nrow(outcomes), 0, 10)
  x <- tabulated_rv(family, values)
  found <- sup_lb(x, null, n_seeds = 200L, n_restarts = 25L)

  # E_theta[X] for theta supplied as rows.
  expectations <- function(theta) {
    as.vector(crossprod(
      exp(compile_loglik(family, outcomes)(theta)),
      values
    ))
  }

  expect_equal(found@sup_lb, expectations(matrix(found@theta, nrow = 1L)))

  expect_true(contains(parts(null@region)[[found@part]], found@theta))

  vertex_best <- max(vapply(
    parts(null@region),
    function(s) max(expectations(s@vertices)),
    numeric(1L)
  ))
  expect_gt(found@sup_lb, vertex_best)

  expect_gt(found@sup_lb / certify(x, null, tol = 1e-12)@sup_ub, 0.999)
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
      one <- sup_lb(x, null, n_seeds = 1L, n_restarts = 1L)@sup_lb
      many <- sup_lb(x, null, n_seeds = 100L, n_restarts = 10L)@sup_lb
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

test_that("sup_lb() rejects an `engine` that is not an engine spec", {
  null <- plurality_null(n = 4L, k = 3L)
  x <- likelihood(null@family(c(0.4, 0.35, 0.25)))
  expect_error(sup_lb(x, null, engine = gh_engine), "must be an engine spec")
})

test_that("sup_lb() searches a sample space it cannot enumerate", {
  # `E_theta[p_0]` is the integral of one `N(0, I)` density against another,
  # i.e. a `N(0, 2I)` density evaluated at `theta`. Over the halfspace
  # `{theta_1 <= theta_2}` that is largest at the origin, on the boundary,
  # where it is `1 / (4 pi)`.
  set.seed(11)
  family <- gaussian_family(d = 2L)
  null <- null_model(
    family,
    halfspace_region(normal = c(1, -1), offset = 0)
  )
  x <- likelihood(family(c(0, 0)))

  found <- sup_lb(
    x,
    null,
    engine = gh_engine(n_nodes = 9L),
    n_seeds = 20L,
    n_restarts = 2L
  )

  expect_equal(found@sup_lb, 1 / (4 * pi), tolerance = 1e-3)
  expect_equal(found@theta, c(0, 0), tolerance = 1e-3)
  expect_equal(found@log_sup_lb, log(found@sup_lb))
})

test_that("sup_lb() reads a ratio at nodes where the ratio itself is NaN", {
  # `E_theta[q / p]` for two Gaussian densities is `exp(theta'a + |a|^2 / 2 +
  # (|m_p|^2 - |m_q|^2) / 2)` with `a = m_q - m_p`: an exponential of a linear
  # function of theta, whatever the distance involved.
  family <- gaussian_family(d = 2L)
  m_q <- c(0.5, -0.25)
  m_p <- c(0, 0.25)
  x <- likelihood(family(m_q)) / likelihood(family(m_p))
  far <- c(40, 40)
  null <- null_model(family, point_region(theta = far))
  a <- m_q - m_p
  truth <- sum(far * a) + sum(a^2) / 2 + (sum(m_p^2) - sum(m_q^2)) / 2

  engine <- gh_engine(n_nodes = 9L)
  found <- sup_lb(x, null, engine = engine)
  expect_equal(found@log_sup_lb, truth)
  expect_equal(found@sup_lb, exp(truth))

  # The same variable evaluated directly, which is what the log form is for:
  # both densities have underflowed at nodes this far out, and every node
  # comes back `0 / 0`.
  nodes <- ripr_init(family(far), null, engine = engine)@engine@nodes
  expect_true(all(is.nan(x(nodes))))

  # And with the log form stripped off there is nothing left to report.
  direct <- random_variable(
    function(y) x(y),
    sample_space = family@sample_space
  )
  expect_error(sup_lb(direct, null, engine = engine), "is `NaN` at every part")
})

test_that("log space changes the range of sup_lb(), not its answer", {
  set.seed(12)
  null <- plurality_null(n = 8L, k = 3L)
  family <- null@family
  x <- likelihood(family(c(0.5, 0.3, 0.2)))
  outcomes <- enumerate_space(family@sample_space)
  # The same variable, tabulated, so it carries no log form and is integrated
  # directly rather than in log space.
  direct <- tabulated_rv(family, x(outcomes))

  set.seed(13)
  in_log <- sup_lb(x, null, n_seeds = 50L, n_restarts = 5L)
  set.seed(13)
  linear <- sup_lb(direct, null, n_seeds = 50L, n_restarts = 5L)

  expect_equal(in_log@sup_lb, linear@sup_lb)
  expect_equal(in_log@theta, linear@theta)
  expect_equal(in_log@log_sup_lb, linear@log_sup_lb)
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
  expect_length(res@bounds, 1L)
  expect_identical(res@method, "bernstein")
  expect_true(all(res@converged))

  # One part, two cells, and the certificate is reduced back onto the part.
  expect_length(null@cells, 2L)
  expect_identical(null@cell_part, c(1L, 1L))

  # One-sidedness, against a dense grid over the square. The grid is a lower
  # bound on the true supremum, so this is necessary and not sufficient, which
  # is the same standard the simplex cases here are held to.
  chart <- chart(null@cells[[1L]])
  grid <- do.call(
    rbind,
    lapply(null@cells, function(cell) {
      w <- as.matrix(expand.grid(rep(list(seq(0, 1, length.out = 21L)), 2L)))
      w <- cbind(w, 1 - rowSums(w))
      w <- w[rowSums(w >= 0) == 3L, , drop = FALSE]
      w %*% cell@vertices
    })
  )
  outcomes <- enumerate_space(family@sample_space)
  values <- x(outcomes)
  attained <- as.vector(
    crossprod(exp(compile_loglik(family, outcomes)(grid)), values)
  )
  expect_lte(max(attained), res@sup_ub)
  expect_gte(res@sup_lb, max(attained) - 1e-6)
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
    shared@iterations,
    sum(vapply(alone, \(r) r@iterations, integer(1)))
  )

  # And the certificate is unaffected. Pruning against a value that was
  # actually attained cannot drop the maximiser, and the bound each pruned run
  # returns still accounts for what it dropped.
  separate_ub <- max(vapply(alone, \(r) r@sup_ub, numeric(1)))
  expect_equal(shared@sup_ub, separate_ub, tolerance = 1e-6)
  expect_equal(
    shared@sup_lb,
    max(vapply(alone, \(r) r@sup_lb, numeric(1)))
  )
})


test_that("a dominated part reports its own bound, not the incumbent's", {
  family <- multinomial_family(n_trials = 10L, k = 3L)
  outcomes <- enumerate_space(family@sample_space)
  x <- tabulated_rv(family, outcomes[, 1L]^2)
  large <- simplex_region(
    vertices = rbind(c(1, 0, 0), c(0.5, 0.5, 0), c(0.5, 0, 0.5))
  )
  small <- simplex_region(
    vertices = rbind(c(0, 0, 1), c(0, 0.5, 0.5), c(0.5, 0, 0.5))
  )
  res <- certify(x, null_model(family, list(large, small)), tol = 1e-9)
  alone <- certify(x, null_model(family, list(small)), tol = 1e-9)

  # The large part's attained value dominates everything the poor part has.
  expect_gt(res@incumbents[1L], alone@sup_ub)
  # The small part still reports a valid bound on itself, far below the
  # incumbent it was pruned against.
  expect_gte(res@bounds[2L], alone@sup_lb)
  expect_lte(res@bounds[2L], res@incumbents[1L])
})


test_that("the incumbent carries from the point cells to the Bernstein runs", {
  # Nothing about the reduction stops a value attained under one method from
  # pruning a search under another: `sup_lb` is a maximum over every cell, so
  # any cell's incumbent bounds the supremum from below. Points are evaluated
  # first, whatever order the parts were declared in, so the Bernstein runs are
  # handed what they attained.
  set.seed(110)
  family <- multinomial_family(n_trials = 8L, k = 3L)
  outcomes <- enumerate_space(family@sample_space)
  # Large only where theta_1 is large, so the point attains far more than the
  # facet's seed enclosure, while the facet's own supremum is interior.
  values <- stats::runif(nrow(outcomes), 0, 10) + 100 * (outcomes[, 1L] == 8L)
  x <- tabulated_rv(family, values)
  facet <- diag(3L)
  facet[1L, ] <- c(0.5, 0.5, 0)
  point <- point_region(theta = c(0.9, 0.05, 0.05))

  # A null of simplices alone starts from nothing, so it has to subdivide.
  alone <- certify(
    x,
    null_model(family, list(simplex_region(vertices = facet))),
    tol = 1e-9
  )
  expect_gt(alone@iterations, 0L)

  # Declared after the facet, the point is still evaluated first, and what it
  # attains already beats the facet's seed bound: nothing is left to split.
  res <- certify(
    x,
    null_model(family, list(simplex_region(vertices = facet), point)),
    tol = 1e-9
  )
  expect_identical(res@method, c("bernstein", "point"))
  expect_gt(res@incumbents[[2L]], res@bounds[[1L]])
  expect_identical(res@iterations[[1L]], 0L)
  expect_true(all(res@converged))
  # The pruned facet's bound is looser for it, but still a bound.
  expect_gte(res@bounds[[1L]], alone@sup_ub)
  expect_equal(res@sup_ub, res@incumbents[[2L]])
})


# --- Result objects and e_variable() ------------------------------------------

test_that("certify() and sup_lb() return classed results that print", {
  set.seed(120)
  null <- plurality_null(n = 4L, k = 3L)
  values <- stats::runif(nrow(enumerate_space(null@family@sample_space)), 0, 3)
  x <- tabulated_rv(null@family, values)

  cert <- certify(x, null, tol = 1e-9)
  expect_true(S7::S7_inherits(cert, ripr_certificate))
  expect_identical(cert@random_variable, x)
  expect_identical(cert@null, null)
  expect_length(cert@bounds, 2L)
  expect_match(
    format(cert),
    "ripr_certificate: sup E[<tabulated>] <= ",
    fixed = TRUE
  )
  out <- paste(capture.output(print(cert)), collapse = "\n")
  expect_match(out, "<ripr_certificate>", fixed = TRUE)
  expect_match(out, "(certified, by bernstein)", fixed = TRUE)
  expect_match(out, "part 2", fixed = TRUE)
  expect_match(out, "e_variable(): X / ", fixed = TRUE)
  expect_no_match(out, "@")

  found <- sup_lb(x, null, n_seeds = 10L, n_restarts = 2L)
  expect_true(S7::S7_inherits(found, ripr_search))
  expect_match(
    format(found),
    "ripr_search: sup E[<tabulated>] >= ",
    fixed = TRUE
  )
  out <- paste(capture.output(print(found)), collapse = "\n")
  expect_match(out, "searched, not certified", fixed = TRUE)
  expect_match(out, paste0("in part ", found@part), fixed = TRUE)
})

test_that("a certificate says in print when a search ran out of nodes", {
  set.seed(121)
  null <- plurality_null(n = 8L, k = 3L)
  x <- tabulated_rv(
    null@family,
    stats::runif(nrow(enumerate_space(null@family@sample_space)), 0, 10)
  )
  starved <- certify(x, null, tol = 0, max_splits = 1L)
  expect_true(any(starved@budget_hit))
  expect_match(
    paste(capture.output(print(starved)), collapse = "\n"),
    "node budget reached",
    fixed = TRUE
  )
})

test_that("e_variable() leaves an e-variable alone and rescales anything else", {
  set.seed(122)
  null <- plurality_null(n = 4L, k = 3L)
  family <- null@family
  outcomes <- enumerate_space(family@sample_space)

  # Values in [0, 1] cannot have expectation above 1 anywhere.
  small <- tabulated_rv(family, stats::runif(nrow(outcomes), 0, 1))
  cert <- certify(small, null, tol = 1e-9)
  expect_lte(cert@sup_ub, 1)
  expect_identical(e_variable(cert), small)
  expect_match(
    paste(capture.output(print(cert)), collapse = "\n"),
    "already an e-variable",
    fixed = TRUE
  )

  # A likelihood ratio against a null point that is not the RIPr overshoots,
  # and is divided by exactly the bound, keeping its log form.
  x <- likelihood(family(c(0.4, 0.35, 0.25))) /
    likelihood(family(c(1, 1, 1) / 3))
  cert <- certify(x, null, tol = 1e-9)
  expect_gt(cert@sup_ub, 1)
  e <- e_variable(cert)
  expect_true(S7::S7_inherits(e, random_variable))
  expect_equal(e(outcomes), x(outcomes) / cert@sup_ub)
  expect_false(is.null(e@log_f))
  expect_lte(certify(e, null, tol = 1e-9)@sup_ub, 1 + 1e-9)
})

test_that("e_variable() refuses a searched lower bound", {
  null <- plurality_null(n = 4L, k = 3L)
  x <- likelihood(null@family(c(0.4, 0.35, 0.25)))
  found <- sup_lb(x, null, n_seeds = 5L, n_restarts = 1L)
  expect_error(e_variable(found), "lower bound")
  expect_error(e_variable(1))
})


test_that("a one-category null is refused with a reason, not an assertion", {
  fam <- multinomial_family(n_trials = 3L, k = 1L)
  x <- random_variable(function(x) rep(2, nrow(x)), fam@sample_space)
  null <- null_model(fam, simplex_region(vertices = matrix(1)))
  expect_error(certify(x, null), "single coordinate")
})
