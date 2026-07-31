# The property the whole exercise exists for: with a certified gap, the
# rescaled e-variable's expectation under *every* null theta is at most 1 --
# not at the thetas an oracle happened to visit, but everywhere.
#
# The comparison against the heuristic gap is the point of the first test. It is
# not that multi-start BFGS is bad; it is that "no counterexample found" and
# "no counterexample exists" are different claims, and only one of them is what
# a certificate should assert.

# A random small multinomial problem on the plurality null.
random_problem <- function(n, K = 3L) {
  comp <- matrix(runif(K * 2L, 0.5, 3), nrow = K)
  comp <- sweep(comp, 2L, colSums(comp), "/")
  # Keep candidate 1 leading, so Q sits inside the plurality alternative.
  comp[1L, ] <- comp[1L, ] + 0.6
  comp <- sweep(comp, 2L, colSums(comp), "/")
  fam <- multinomial_family(n_trials = n, k = K)
  alt <- finite_mixing(components = comp, weights = c(0.5, 0.5))
  faces <- plurality_faces(K)
  ripr_problem(fam, null_region(faces = faces), as_marginal(alt, fam))
}

# E_theta[(Q / P*) / (1 + gap)] by direct summation over the enumerated support.
null_expectation <- function(ev, fam, theta) {
  X <- support(fam)
  sum(exp(log_density(fam, theta) + e_value(ev, X, log = TRUE)))
}

# Uniform draws from a face: Dirichlet weights on its vertices.
sample_face <- function(face, n_pts) {
  V <- face@vertices
  a <- matrix(rgamma(n_pts * ncol(V), shape = 1), nrow = n_pts)
  a <- a / rowSums(a)
  V %*% t(a)
}

# A dense grid over the K = 3 probability simplex, restricted to the null.
null_grid <- function(faces, step = 0.01) {
  g <- as.matrix(expand.grid(
    a = seq(0, 1, by = step),
    b = seq(0, 1, by = step)
  ))
  g <- g[rowSums(g) <= 1 + 1e-12, , drop = FALSE]
  theta <- cbind(g, pmax(1 - rowSums(g), 0))
  keep <- apply(theta, 1L, function(th) {
    any(vapply(faces, function(f) contains(f, th, tol = 1e-9), logical(1L)))
  })
  theta[keep, , drop = FALSE]
}

fit_and_certify <- function(prob, bnb = bnb_control(tol = 1e-6)) {
  ref <- alt_mean(prob$alternative@mixing)
  faces <- prob$null
  res <- run_ripr(
    prob,
    init_atoms = init_on_faces(faces, ref),
    init_atom_faces = seq_along(faces),
    fw_iters = 12L,
    em_iters = 10L,
    n_seeds = 100L,
    verbose = FALSE
  )
  list(
    res = res,
    heuristic = res$certificate,
    certified = certify(res$projection, prob, bnb = bnb, n_seeds = 100L)
  )
}

test_that("the certified gap makes the null expectation at most 1 everywhere", {
  set.seed(1234)
  for (trial in seq_len(4L)) {
    n <- sample(3:6, 1L)
    prob <- random_problem(n, K = 3L)
    fit <- fit_and_certify(prob)
    expect_true(fit$certified$gap_certified)

    ev <- e_variable(
      numerator = prob$alternative,
      projection = fit$res$projection,
      gap = fit$certified$gap_used,
      gap_certified = TRUE
    )

    thetas <- do.call(cbind, lapply(prob$null, sample_face, n_pts = 150L))
    vals <- apply(thetas, 2L, function(th) {
      null_expectation(ev, prob$family, th)
    })
    expect_true(
      all(vals <= 1 + 1e-10),
      info = sprintf("trial %d (n = %d): max null expectation %.12f",
                     trial, n, max(vals))
    )
  }
})

test_that("the heuristic gap carries no such guarantee (recorded counterexample)", {
  # The case the accuracy section of ?certify warns about, made concrete and
  # deterministic rather than left as folklore.
  #
  # Construction. With P_W a single atom at theta_W and Q a two-component
  # multinomial mixture, the ratio coefficients collapse and
  #
  #   G(theta) = 0.5 <theta, r_a>^n + 0.5 <theta, r_b>^n,   r = theta_alt / theta_W
  #
  # which is *convex*, so every local maximum sits at a vertex of the face. The
  # two components are chosen so that two different vertices are local maxima at
  # different heights: e_3 globally (r_a is third-heavy) and e_2 locally (r_b is
  # second-heavy). BFGS ascends within whichever basin its seed lands in, so an
  # under-seeded oracle returns the *lower* peak and reports a gap that is 23%
  # too small. Nothing about the search is broken -- a lower bound is all it
  # ever promised.
  n <- 20L
  K <- 3L
  fam <- multinomial_family(n_trials = n, k = K)
  alt <- finite_mixing(
    components = cbind(c(0.01, 0.01, 0.98), c(0.01, 0.97, 0.02)),
    weights = c(0.5, 0.5)
  )
  face <- plurality_faces(K)[[1L]]
  prob <- ripr_problem(fam, null_region(faces = list(face)), as_marginal(alt, fam))
  uniform <- finite_mixing(
    components = matrix(rep(1 / K, K), ncol = 1L),
    weights = 1
  )

  cert <- certify(uniform, prob, bnb = bnb_control(tol = 1e-6), n_seeds = 1L)
  set.seed(2)
  heur <- certify(uniform, prob, n_seeds = 1L)

  expect_true(cert$gap_certified)
  expect_false(heur$gap_certified)
  # The oracle settled on the lower peak; the bound found the higher one.
  expect_gt(cert$gap_used / heur$gap_used, 1.2)

  ev_of <- function(g, certified) {
    e_variable(
      numerator = prob$alternative,
      projection = as_marginal(uniform, fam),
      gap = g,
      gap_certified = certified
    )
  }
  # The three face vertices plus a spread of interior points. The vertices are
  # where a convex G attains its maxima, so they are exactly the points a
  # vertex-blind search can under-report.
  set.seed(5)
  thetas <- cbind(face@vertices, sample_face(face, 200L))
  worst <- function(ev) {
    max(apply(thetas, 2L, function(th) null_expectation(ev, fam, th)))
  }

  # Certified: at most 1 under every null theta, as a proven inequality.
  expect_lte(worst(ev_of(cert$gap_used, TRUE)), 1 + 1e-10)
  # Heuristic: not merely untight -- actually invalid, by 23%.
  expect_gt(worst(ev_of(heur$gap_used, FALSE)), 1.2)
})

test_that("the bound dominates the heuristic oracle and a dense grid", {
  set.seed(202)
  prob <- random_problem(5L, K = 3L)
  fit <- fit_and_certify(prob)

  state <- ripr:::build_mixture_state(
    prob$engine,
    lapply(
      seq_len(ncol(fit$res$projection@mixing@components)),
      function(k) fit$res$projection@mixing@components[, k]
    ),
    rep(1L, ncol(fit$res$projection@mixing@components)),
    fit$res$projection@mixing@weights
  )
  log_Pw <- ripr:::state_log_p_mixture(state)
  orr <- ripr:::oracle_step(state, prob, n_seeds = 400L)

  bounds <- lapply(prob$null, function(f) {
    oracle_bound(f, prob$engine, log_Pw, prob$family,
                 control = bnb_control(tol = 1e-8))
  })
  b <- vapply(bounds, `[[`, "bound", FUN.VALUE = numeric(1L))
  inc <- vapply(bounds, `[[`, "incumbent", FUN.VALUE = numeric(1L))

  # The bound is an upper bound; the oracle's E_star is a lower one.
  expect_gte(max(b), orr$E_star)
  # incumbent <= bound, and the incumbent is attained at a real point.
  expect_true(all(inc <= b))

  grid <- null_grid(prob$null, step = 0.01)
  gvals <- apply(grid, 1L, function(th) {
    ripr:::expect_ratio(prob$engine, log_density(prob$family, th), log_Pw)
  })
  # The bound dominates the grid maximum, to within floating point. The slack
  # is relative and deliberate: `gvals` come from the engine's log-space sum
  # while `b` comes from de Casteljau, and the two agree only to about 1e-14
  # relative (see the floating point section of ?oracle_bound). Asserting exact
  # domination would encode a claim that is false in the ill-conditioned tail.
  expect_gte(max(b), max(gvals) * (1 - 1e-12))
  # The incumbent is *not* asserted to sit below the grid maximum: it is an
  # exact value of G at a point the grid need not contain, so it may
  # legitimately be the larger of the two.
  expect_gte(max(inc), max(gvals) - 1e-3)

  # Every returned theta lies on its own face, and its incumbent is the exact
  # value of G there -- not an approximation to it.
  for (i in seq_along(prob$null)) {
    expect_true(contains(prob$null[[i]], bounds[[i]]$theta, tol = 1e-9))
    expect_equal(
      ripr:::expect_ratio(
        prob$engine,
        log_density(prob$family, bounds[[i]]$theta),
        log_Pw
      ),
      inc[[i]],
      tolerance = 1e-9
    )
  }
  expect_identical(bounds[[1L]]$method, "bernstein_bnb")
})

test_that("keep_argmax returns an enclosure of every maximiser", {
  set.seed(31)
  prob <- random_problem(6L, K = 3L)
  fit <- fit_and_certify(prob)
  mix <- fit$res$projection@mixing
  state <- ripr:::build_mixture_state(
    prob$engine,
    lapply(seq_len(ncol(mix@components)), function(k) mix@components[, k]),
    rep(1L, ncol(mix@components)),
    mix@weights
  )
  log_Pw <- ripr:::state_log_p_mixture(state)

  # Without keep_argmax the pruning rule discards boxes that may attain the
  # maximum, so no enclosure is claimed -- and none is returned.
  plain <- oracle_bound(prob$null[[1L]], prob$engine, log_Pw, prob$family,
                        control = bnb_control(tol = 1e-6))
  expect_null(plain$active)

  ob <- oracle_bound(prob$null[[1L]], prob$engine, log_Pw, prob$family,
                     control = bnb_control(tol = 1e-6, keep_argmax = TRUE))
  expect_true(length(ob$active) > 0L)

  # Barycentric coordinates of a point in a box; all non-negative means inside.
  inside <- function(box, theta) {
    a <- tryCatch(qr.solve(box$V, theta), error = function(e) NULL)
    !is.null(a) &&
      max(abs(as.vector(box$V %*% a) - theta)) < 1e-9 &&
      all(a >= -1e-9)
  }

  # Find the grid argmax over this face and check the enclosure covers it --
  # unless the grid point is not a maximiser at all, in which case excluding it
  # is the correct behaviour, not a failure.
  face <- prob$null[[1L]]
  grid <- null_grid(list(face), step = 0.01)
  gvals <- apply(grid, 1L, function(th) {
    ripr:::expect_ratio(prob$engine, log_density(prob$family, th), log_Pw)
  })
  arg <- grid[which.max(gvals), ]
  covered <- any(vapply(ob$active, inside, logical(1L), theta = arg))
  expect_true(covered || max(gvals) < ob$incumbent)
})

test_that("oracle_bound() refuses a mixture with no support where Q has mass", {
  set.seed(6)
  fam <- multinomial_family(n_trials = 4L, k = 3L)
  alt <- point_mixing(theta_star = c(0.5, 0.3, 0.2))
  faces <- plurality_faces(3L)
  prob <- ripr_problem(fam, null_region(faces = faces), as_marginal(alt, fam))

  # A single atom on a simplex vertex: P_W gives zero mass to most outcomes
  # that Q reaches, so G is unbounded and no rescaling saves it.
  degenerate <- ripr:::build_mixture_state(
    prob$engine, list(c(0, 1, 0)), 1L, 1
  )
  log_Pw <- ripr:::state_log_p_mixture(degenerate)
  expect_error(
    oracle_bound(faces[[1L]], prob$engine, log_Pw, prob$family),
    "no rescaling"
  )
})

test_that("halfspace faces have no bound method, permanently", {
  hs <- halfspace_face(v = c(1, -1), c = 0, face_index = 1)
  fam <- multinomial_family(n_trials = 3L, k = 2L)
  alt <- point_mixing(theta_star = c(0.6, 0.4))
  eng <- exact_engine(fam, as_marginal(alt, fam))
  expect_false(certifiable(hs, eng, fam))
  expect_error(oracle_bound(hs, eng, rep(log(0.5), eng@M), fam))
})
