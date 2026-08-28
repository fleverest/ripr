# Properties of R/steps.R.
#
# The identity `sum_c w_c G(theta_c) = 1` is the workhorse. It is algebraic, not
# a convergence condition, so it holds at every iterate and fails only if the
# quadrature weights are unnormalised, the atoms and weights are misaligned, or
# `log_p` was computed aGst the wrong mixture.
#
# Most of what is here is exercised again through the verbs in test-fit.R. These
# tests aim at the pieces directly, so that a failure says which piece.

plurality <- function(n = 12, k = 4, q = c(0.42, 0.31, 0.16, 0.11), ...) {
  fam <- multinomial_family(n_trials = n, k = k)
  Q <- induced_distribution(fam, point_mixing(theta_star = q))
  parts <- lapply(2:k, function(j) {
    basis <- lapply(setdiff(seq_len(k), 1L), function(i) {
      replace(numeric(k), i, 1)
    })
    tie <- replace(numeric(k), c(1L, j), 0.5)
    simplex_region(vertices = do.call(cbind, c(basis, list(tie))))
  })
  ripr_init(
    Q,
    null_model(fam, parts),
    exact_engine(),
    control = ripr_control(n_seeds = 30L, n_restarts = 4L, ...)
  )
}

# Everything below `plan_step` works on `(ld_all, w)`, so most tests need the
# candidate column spliced in at a known index.
with_candidate <- function(state, theta, at = NULL) {
  ld <- compile_engine(state@engine)
  w <- flat_weights(state)
  at <- if (is.null(at)) length(w) + 1L else at
  list(
    ld_all = insert_col(
      ld(flat_atoms(state)),
      as.vector(ld(matrix(theta, ncol = 1L))),
      at
    ),
    w = append(w, 0, after = at - 1L),
    new_idx = at,
    log_p = log_p_at_nodes(state, ld),
    engine = state@engine
  )
}

Gs <- function(ld_all, w, engine) {
  exp(col_logsumexp(ld_all - mixture_log_p(ld_all, w) + engine@log_w))
}

# --- The algebraic identity ---------------------------------------------------

test_that("weighted Gs sum to 1 at any weights", {
  # sum_c w_c E_Q[p_c / P_W] = E_Q[P_W / P_W] = 1. Algebra, so it does not care
  # whether the weights are any good.
  st <- plurality()
  cand <- with_candidate(st, c(0.30, 0.40, 0.20, 0.10))
  for (w in list(cand$w, rev(cand$w), c(0.7, 0.1, 0.1, 0.1))) {
    expect_equal(
      sum(w * Gs(cand$ld_all, w, cand$engine)),
      1,
      tolerance = rounding_tol(1)
    )
  }
})


# --- Step paths ---------------------------------------------------------------

test_that("every path stays on the simplex", {
  st <- plurality()
  cand <- with_candidate(st, c(0.30, 0.40, 0.20, 0.10))
  paths <- step_paths(
    c("forward", "pairwise", "away"),
    cand$ld_all,
    cand$w,
    cand$new_idx,
    cand$log_p,
    cand$engine
  )
  for (p in paths) {
    for (gamma in seq(0, p$gamma_max, length.out = 5L)) {
      w <- p$w_of(gamma)
      expect_equal(sum(w), 1, tolerance = rounding_tol(1))
      expect_true(all(w >= -rounding_tol(1)))
    }
  }
})

test_that("each path's log density agrees with rebuilding the mixture", {
  # `path_forward` takes a two-column shortcut rather than rebuilding, which is
  # only sound if it agrees with the long way round.
  st <- plurality()
  cand <- with_candidate(st, c(0.30, 0.40, 0.20, 0.10))
  paths <- step_paths(
    c("forward", "pairwise", "away"),
    cand$ld_all,
    cand$w,
    cand$new_idx,
    cand$log_p,
    cand$engine
  )
  for (p in paths) {
    gamma <- p$gamma_max / 3
    expect_equal(p$log_p_at(gamma), mixture_log_p(cand$ld_all, p$w_of(gamma)))
  }
})

test_that("pairwise coincides with forward when one atom is active", {
  # `worst` is then that atom, and both paths are `[1 - gamma, gamma]` with a
  # cap of 1. So the one-atom case needs no special handling.
  st <- plurality()
  ld <- compile_engine(st@engine)
  w <- c(1, rep(0, length(flat_weights(st)) - 1L), 0)
  ld_all <- cbind(
    ld(flat_atoms(st)),
    as.vector(ld(matrix(c(0.3, 0.4, 0.2, 0.1), ncol = 1L)))
  )
  new_idx <- length(w)
  log_p <- mixture_log_p(ld_all, w)

  fwd <- step_paths("forward", ld_all, w, new_idx, log_p, st@engine)[[1L]]
  pw <- step_paths("pairwise", ld_all, w, new_idx, log_p, st@engine)[[1L]]
  expect_equal(fwd$gamma_max, pw$gamma_max)
  for (gamma in c(0, 0.25, 0.5, 1)) {
    expect_equal(fwd$w_of(gamma), pw$w_of(gamma))
  }
})

test_that("away alone errors below two active atoms", {
  # There is nothing to redistribute the mass to. A silent fallback is what
  # previously made `"away"` behave as pairwise while its docs claimed forward.
  st <- plurality()
  ld <- compile_engine(st@engine)
  w <- c(1, rep(0, length(flat_weights(st)) - 1L), 0)
  ld_all <- cbind(
    ld(flat_atoms(st)),
    as.vector(ld(matrix(c(0.3, 0.4, 0.2, 0.1), ncol = 1L)))
  )
  expect_error(
    step_paths(
      "away",
      ld_all,
      w,
      length(w),
      mixture_log_p(ld_all, w),
      st@engine
    ),
    "second active atom"
  )
  # Offered alongside forward it is simply dropped.
  expect_length(
    step_paths(
      c("forward", "away"),
      ld_all,
      w,
      length(w),
      mixture_log_p(ld_all, w),
      st@engine
    ),
    1L
  )
})

test_that("the away step never uses the candidate", {
  st <- plurality()
  cand <- with_candidate(st, c(0.30, 0.40, 0.20, 0.10))
  away <- step_paths(
    "away",
    cand$ld_all,
    cand$w,
    cand$new_idx,
    cand$log_p,
    cand$engine
  )[[1L]]
  expect_equal(away$w_of(away$gamma_max / 2)[cand$new_idx], 0)
})

test_that("pairwise and away empty the worst atom at their cap", {
  st <- plurality()
  cand <- with_candidate(st, c(0.30, 0.40, 0.20, 0.10))
  worst <- worst_atom(cand$ld_all, cand$w, cand$log_p, cand$engine)
  for (d in c("pairwise", "away")) {
    p <- step_paths(
      d,
      cand$ld_all,
      cand$w,
      cand$new_idx,
      cand$log_p,
      cand$engine
    )[[1L]]
    expect_equal(p$w_of(p$gamma_max)[worst], 0, tolerance = rounding_tol(1))
  }
})

# --- Step sizes ---------------------------------------------------------------

test_that("a line search cannot increase KL", {
  # `gamma = 0` is in every interval, so the step can always decline to move.
  st <- plurality()
  cand <- with_candidate(st, c(0.30, 0.40, 0.20, 0.10))
  before <- expect_q(cand$engine, cand$engine@log_q - cand$log_p)
  for (d in list("forward", "pairwise", c("forward", "away"))) {
    res <- apply_step(
      cand$ld_all,
      cand$w,
      cand$new_idx,
      cand$log_p,
      cand$engine,
      directions = d
    )
    expect_true(res$kl <= before + 1e-12)
  }
})

test_that("the fixed schedule is capped at the path's own maximum", {
  # Pairwise and away cap below 1; an uncapped schedule value would take the
  # weights off the simplex.
  st <- plurality()
  cand <- with_candidate(st, c(0.30, 0.40, 0.20, 0.10))
  res <- apply_step(
    cand$ld_all,
    cand$w,
    cand$new_idx,
    cand$log_p,
    cand$engine,
    directions = "pairwise",
    size = "fixed",
    gamma_fixed = 1
  )
  expect_equal(sum(res$weights), 1, tolerance = rounding_tol(1))
  expect_true(all(res$weights >= 0))
})

test_that("apply_step takes whichever offered direction reaches the lowest KL", {
  st <- plurality()
  cand <- with_candidate(st, c(0.30, 0.40, 0.20, 0.10))
  single <- vapply(
    c("forward", "pairwise", "away"),
    function(d) {
      apply_step(
        cand$ld_all,
        cand$w,
        cand$new_idx,
        cand$log_p,
        cand$engine,
        directions = d
      )$kl
    },
    numeric(1)
  )
  both <- apply_step(
    cand$ld_all,
    cand$w,
    cand$new_idx,
    cand$log_p,
    cand$engine,
    directions = c("forward", "pairwise", "away")
  )
  expect_equal(both$kl, min(single))
  expect_identical(both$direction, names(single)[which.min(single)])
})

test_that("uses_candidate is derived from the weight it ends with", {
  st <- plurality()
  cand <- with_candidate(st, c(0.30, 0.40, 0.20, 0.10))
  fwd <- apply_step(cand$ld_all, cand$w, cand$new_idx, cand$log_p, cand$engine)
  expect_identical(fwd$uses_candidate, fwd$weights[cand$new_idx] > 0)
  away <- apply_step(
    cand$ld_all,
    cand$w,
    cand$new_idx,
    cand$log_p,
    cand$engine,
    directions = "away"
  )
  expect_false(away$uses_candidate)
})

# --- The corrective solve -----------------------------------------------------

test_that("a weight sweep is monotone and preserves the identity", {
  # `w <- w * G` is the exact M-step for the weights, so `sum_c w_c G_c = 1`
  # makes it self-normalising.
  st <- plurality()
  ld <- compile_engine(st@engine)
  ld_all <- ld(flat_atoms(st))
  w <- flat_weights(st)
  kl <- expect_q(st@engine, st@engine@log_q - mixture_log_p(ld_all, w))
  for (i in 1:5) {
    w <- weight_sweep(ld_all, w, mixture_log_p(ld_all, w), st@engine)$weights
    expect_equal(sum(w), 1, tolerance = rounding_tol(1))
    kl_new <- expect_q(st@engine, st@engine@log_q - mixture_log_p(ld_all, w))
    expect_true(kl_new <= kl + rounding_tol(kl))
    kl <- kl_new
  }
})

test_that("the sweep residual is the restricted Frank-Wolfe gap", {
  st <- plurality()
  ld <- compile_engine(st@engine)
  ld_all <- ld(flat_atoms(st))
  w <- flat_weights(st)
  sweep <- weight_sweep(ld_all, w, mixture_log_p(ld_all, w), st@engine)
  expect_equal(sweep$residual, max(Gs(ld_all, w, st@engine)) - 1)
})

test_that("solve_weights respects both its budget and its tolerance", {
  st <- plurality()
  ld <- compile_engine(st@engine)
  ld_all <- ld(flat_atoms(st))
  w0 <- flat_weights(st)
  # A tolerance nothing can meet leaves the cap as the only stopping rule.
  tight <- solve_weights(ld_all, w0, st@engine, tol = 0, max_iter = 20L)
  loose <- solve_weights(ld_all, w0, st@engine, tol = 1e10, max_iter = 20L)
  expect_equal(loose, w0)
  expect_equal(sum(tight), 1, tolerance = rounding_tol(1))
})

# --- Oracles ------------------------------------------------------------------

test_that("the linear oracle returns G, and its batch form agrees", {
  st <- plurality()
  ld <- compile_engine(st@engine)
  log_p <- log_p_at_nodes(st, ld)
  obj <- linear_oracle(st, log_p, ld)
  atoms <- flat_atoms(st)
  expect_equal(
    vapply(seq_len(ncol(atoms)), function(i) obj$value(atoms[, i]), numeric(1)),
    obj$value_batch(atoms)
  )
  expect_equal(
    obj$value_batch(atoms),
    Gs(ld(atoms), flat_weights(st), st@engine)
  )
})

test_that("the linear oracle's gradient matches finite differences", {
  st <- plurality()
  ld <- compile_engine(st@engine)
  obj <- linear_oracle(st, log_p_at_nodes(st, ld), ld)
  theta <- c(0.30, 0.40, 0.20, 0.10)
  # A feasible direction on the simplex: components must sum to zero.
  d <- c(1, -1, 0, 0) * 1e-6
  fd <- (obj$value(theta + d) - obj$value(theta - d)) / 2
  expect_equal(sum(obj$grad(theta) * d), fd, tolerance = 1e-6)
})

test_that("the Li-Barron oracle scores by KL after the step", {
  st <- plurality()
  ld <- compile_engine(st@engine)
  log_p <- log_p_at_nodes(st, ld)
  obj <- nonlinear_oracle(st, log_p, ld)
  theta <- c(0.30, 0.40, 0.20, 0.10)
  planned <- plan_step(st, log_p, ld)(theta)
  expect_equal(obj$value(theta), -planned$kl)
})

test_that("away is stripped from the Li-Barron inner directions", {
  # An away path never touches the candidate, so it scores identically for every
  # one and its gradient is exactly zero -- a plateau that strands any restart
  # landing on it. Dropping it costs nothing, since the outer step still weighs
  # away aGst whatever the oracle returns.
  expect_identical(inner_directions(c("forward", "away")), "forward")
  expect_identical(inner_directions("away"), "forward")
  expect_identical(
    inner_directions(c("forward", "pairwise")),
    c("forward", "pairwise")
  )

  st <- plurality()
  ld <- compile_engine(st@engine)
  log_p <- log_p_at_nodes(st, ld)
  theta <- c(0.30, 0.40, 0.20, 0.10)
  expect_equal(
    nonlinear_oracle(st, log_p, ld, directions = c("forward", "away"))$value(
      theta
    ),
    nonlinear_oracle(st, log_p, ld, directions = "forward")$value(theta)
  )
})

test_that("the Li-Barron gradient matches finite differences", {
  # Under `size = "fixed"` the step size does not depend on theta, so the
  # objective is smooth. Under a line search it is flat wherever G <= 1, since
  # the search then declines to move at all.
  st <- plurality()
  ld <- compile_engine(st@engine)
  obj <- nonlinear_oracle(
    st,
    log_p_at_nodes(st, ld),
    ld,
    size = "fixed",
    gamma_fixed = 0.4
  )
  theta <- c(0.30, 0.40, 0.20, 0.10)
  d <- c(1, -1, 0, 0) * 1e-6
  fd <- (obj$value(theta + d) - obj$value(theta - d)) / 2
  expect_equal(sum(obj$grad(theta) * d), fd, tolerance = 1e-6)
})

test_that("the Li-Barron objective is flat wherever G is at most 1", {
  # `value(theta) = -min_gamma KL(...)`, and when G(theta) <= 1 the minimising
  # gamma is 0, so the value is the current KL whatever theta is. Not a defect,
  # but it is why the search depends on some seed landing where G > 1.
  st <- plurality()
  ld <- compile_engine(st@engine)
  log_p <- log_p_at_nodes(st, ld)
  G <- linear_oracle(st, log_p, ld)$value
  obj <- nonlinear_oracle(st, log_p, ld)

  poor <- Filter(
    function(th) G(th) < 1,
    lapply(c(0.05, 0.10, 0.15), function(t1) c(t1, (1 - t1) * c(0.5, 0.3, 0.2)))
  )
  expect_true(length(poor) >= 2L)
  values <- vapply(poor, obj$value, numeric(1))
  expect_equal(values, rep(values[1L], length(values)))
})

test_that("plan_step memoises its last call", {
  # `optim()` asks for the value and the gradient at the same point through
  # separate slots, and each evaluation costs a line search.
  st <- plurality()
  calls <- 0L
  ld <- compile_engine(st@engine)
  counted <- function(m) {
    calls <<- calls + 1L
    ld(m)
  }
  obj <- nonlinear_oracle(st, log_p_at_nodes(st, ld), counted)
  theta <- c(0.30, 0.40, 0.20, 0.10)
  invisible(obj$value(theta))
  after_value <- calls
  invisible(obj$grad(theta))
  expect_identical(calls, after_value)
})

test_that("plan_step puts the candidate where it is told", {
  st <- plurality()
  ld <- compile_engine(st@engine)
  at <- insert_index(st, 2L)
  planned <- plan_step(st, log_p_at_nodes(st, ld), ld, at = at)(c(
    0.30,
    0.40,
    0.20,
    0.10
  ))
  expect_identical(planned$new_idx, at)
  expect_length(planned$weights, length(flat_weights(st)) + 1L)
})

# --- Support identification ---------------------------------------------------

test_that("removing an atom adjusts the mixture in place", {
  # `identify_support` tests every candidate atom, so rebuilding the mixture
  # each time would make a pass O(MC^2). `log_p_without` adjusts the existing
  # one instead, in O(M):
  #
  #   log P + log(1 - w_j p_j / P) - log(1 - w_j)  =  log((P - w_j p_j)/(1 - w_j))
  #
  # which is not hard to get wrong so we check against the definition here.
  st <- plurality()
  ld <- compile_engine(st@engine)
  ld_all <- ld(flat_atoms(st))
  w <- flat_weights(st)
  log_p <- mixture_log_p(ld_all, w)
  dropped <- w
  dropped[2L] <- 0
  expect_equal(
    log_p_without(log_p, ld_all[, 2L], w[2L]),
    mixture_log_p(ld_all, dropped / sum(dropped))
  )
})

test_that("the downdate survives an atom carrying nearly all the mass", {
  st <- plurality()
  ld <- compile_engine(st@engine)
  ld_all <- ld(flat_atoms(st))
  n <- ncol(ld_all)

  heavy <- c(1 - 1e-15, rep(1e-15 / (n - 1), n - 1))
  log_p <- mixture_log_p(ld_all, heavy)

  # Removing the dominant atom: `share` is 1 to machine precision throughout,
  # and the renormalising term divides by 1e-15.
  expect_true(all(!is.nan(log_p_without(log_p, ld_all[, 1L], heavy[1L]))))

  # Removing a negligible one must still agree with rebuilding.
  dropped <- heavy
  dropped[2L] <- 0
  expect_equal(
    log_p_without(log_p, ld_all[, 2L], heavy[2L]),
    mixture_log_p(ld_all, dropped / sum(dropped))
  )

  # And the whole pass must come through with usable weights.
  identified <- identify_support(ld_all, heavy, st@engine)
  expect_true(all(!is.na(identified)))
  expect_equal(sum(identified), 1, tolerance = rounding_tol(1))
})

test_that("the downdate is finite at a weight of exactly 1", {
  st <- plurality()
  ld <- compile_engine(st@engine)
  ld_all <- ld(flat_atoms(st))
  log_p <- ld_all[, 1L]
  out <- log_p_without(log_p, ld_all[, 1L], 1)
  expect_true(all(!is.nan(out)))
  expect_true(all(out == -Inf))
})


test_that("identify_support does not increase KL", {
  # It compares two feasible points and keeps the better
  st <- plurality()
  ld <- compile_engine(st@engine)
  ld_all <- ld(flat_atoms(st))
  w <- flat_weights(st)
  before <- expect_q(st@engine, st@engine@log_q - mixture_log_p(ld_all, w))
  after_w <- identify_support(ld_all, w, st@engine)
  expect_equal(sum(after_w), 1, tolerance = rounding_tol(1))
  expect_true(
    expect_q(st@engine, st@engine@log_q - mixture_log_p(ld_all, after_w)) <=
      before + rounding_tol(before)
  )
})

test_that("identify_support leaves at least one atom", {
  st <- plurality()
  ld <- compile_engine(st@engine)
  w <- identify_support(ld(flat_atoms(st)), flat_weights(st), st@engine)
  expect_true(sum(w > 0) >= 1L)
})

# --- EM -----------------------------------------------------------------------

test_that("an EM sweep keeps every atom in its own part", {
  # The M-step optimises each component over its own chart, so an atom cannot
  # migrate between parts however the responsibilities fall.
  st <- fw_step(plurality(), 2L)
  moved <- em_sweep(st, compile_engine(st@engine))
  for (i in seq_along(moved@atoms)) {
    for (j in seq_len(ncol(moved@atoms[[i]]))) {
      expect_true(contains(
        parts(st@null@region)[[i]],
        moved@atoms[[i]][, j],
        tol = 1e-5
      ))
    }
  }
})

test_that("an EM sweep neither grows nor shrinks the support", {
  st <- fw_step(plurality(), 2L)
  moved <- em_sweep(st, compile_engine(st@engine))
  expect_identical(block_sizes(moved), block_sizes(st))
})

test_that("the atom M-step is deterministic", {
  # No random seeding: the same state must give the same moved atoms whatever
  # the RNG is doing, which is what makes a fit reproducible without a seed.
  st <- plurality()
  ld <- compile_engine(st@engine)
  wt <- exp(st@engine@log_w) * em_responsibilities(st, ld)
  set.seed(41)
  a <- em_atom_step(st, ld, wt)
  set.seed(42)
  b <- em_atom_step(st, ld, wt)
  expect_equal(flat_atoms(a), flat_atoms(b))
})

test_that("an EM sweep preserves the identity", {
  st <- fw_step(plurality(), 2L)
  moved <- em_sweep(st, compile_engine(st@engine))
  ld <- compile_engine(moved@engine)
  w <- flat_weights(moved)
  expect_equal(
    sum(w * Gs(ld(flat_atoms(moved)), w, moved@engine)),
    1,
    tolerance = rounding_tol(1)
  )
})

# --- An atom on the boundary of the support -----------------------------------

test_that("the EM M-step survives a node its atom gives zero probability", {
  # `0 * -Inf` is `NaN`, not the 0 the term is worth under `0 log 0 = 0`: a node
  # an atom cannot produce has both zero responsibility and `-Inf` log density,
  # for the same reason. Left unhandled it reaches `maximise_over`, which
  # compares `NaN < Inf` and errors with "missing value where TRUE/FALSE
  # needed", naming nothing that points back here. The plain multinomial
  # reaches this state because the chart covers the closed region: an atom
  # with an exact zero coordinate gives zero probability to every outcome
  # using that category, and EM genuinely produces such atoms whenever a
  # responsibility-weighted optimum lies on a facet.
  fam <- multinomial_family(n_trials = 10L, k = 3L)
  null <- null_model(
    fam,
    lapply(2:3, function(j) {
      v <- diag(3L)
      v[, 1L] <- replace(numeric(3L), c(1L, j), 0.5)
      simplex_region(vertices = v)
    })
  )
  # One interior atom keeps KL finite; the second sits on the facet
  # `theta_3 = 0` of its part, so it gives zero probability to every outcome
  # using the third category.
  st <- ripr_init(
    fam(c(0.40, 0.35, 0.25)),
    null,
    atoms = list(
      cbind(c(1 / 3, 1 / 3, 1 / 3), c(0.5, 0.5, 0)),
      cbind(c(1 / 3, 1 / 3, 1 / 3))
    )
  )

  # The state really does contain such an atom, or this guards nothing.
  ld <- compile_engine(st@engine)
  wt <- exp(st@engine@log_w) * em_responsibilities(st, ld)
  lda <- ld(flat_atoms(st))
  hazard <- which(apply(wt == 0 & !is.finite(lda), 2L, any))
  expect_gte(length(hazard), 1L)
  # And the unguarded objective at that atom is exactly the NaN in question.
  expect_true(is.nan(sum(wt[, hazard[1L]] * lda[, hazard[1L]])))

  after <- em_step(st, times = 5L)
  expect_false(is.na(kl_divergence(after)))
  expect_lte(kl_divergence(after), kl_divergence(st))
})


# --- Cells under the steps ----------------------------------------------------

# A null whose one part is a convex hull rather than a simplex: the square
# `{theta_1 <= 1/2, theta_2 <= 1/2}` in the 2-simplex, which triangulates into
# two cells. Every part of every other null in this file is its own only cell,
# so this is where the cell wiring has anything to do.
square_null <- function(n = 12, q = c(0.15, 0.35, 0.5), ...) {
  fam <- multinomial_family(n_trials = n, k = 3L)
  Q <- induced_distribution(fam, point_mixing(theta_star = q))
  square <- polytope_region(
    vertices = cbind(
      c(0.5, 0.5, 0),
      c(0, 0.5, 0.5),
      c(0, 0, 1),
      c(0.5, 0, 0.5)
    )
  )
  ripr_init(
    Q,
    null_model(fam, list(square)),
    exact_engine(),
    control = ripr_control(n_seeds = 30L, n_restarts = 4L, ...)
  )
}

test_that("search_null searches cells and reports the part they belong to", {
  st <- square_null()
  expect_length(st@null@cells, 2L)

  log_p <- mixture_log_p(
    compile_engine(st@engine)(flat_atoms(st)),
    flat_weights(st)
  )
  found <- search_null(st, linear_oracle(st, log_p, compile_engine(st@engine)))

  # The index is a part's, because that is what an atom is filed under, and
  # `state@atoms` has exactly that many blocks.
  expect_identical(found$part, 1L)
  expect_true(found$part <= n_parts(st@null@region))
  expect_true(contains(parts(st@null@region)[[found$part]], found$theta))

  # And the maximiser is in one of the cells searched, not merely in the hull.
  expect_true(any(vapply(
    st@null@cells,
    \(cell) contains(cell, found$theta, tol = 1e-6),
    logical(1)
  )))
})

test_that("the atom M-step crosses a fan diagonal when the optimum does", {
  # Responsibilities concentrated on the outcome (0, 6, 6) put the M-step
  # optimum at theta = (0, 0.5, 0.5), a vertex of exactly one of the square's
  # two cells. An atom started deep inside the other cell must cross the fan's
  # diagonal to reach it -- which confinement to a containing cell forbade.
  st <- square_null()
  target <- c(0, 0.5, 0.5)
  holds <- vapply(
    st@null@cells,
    \(cell) contains(cell, target, tol = rounding_tol(1)),
    logical(1)
  )
  expect_identical(sum(holds), 1L)
  far_cell <- st@null@cells[[which(!holds)[1L]]]
  st@atoms[[1L]][, 1L] <- rowMeans(far_cell@vertices)

  ld <- compile_engine(st@engine)
  wt <- exp(st@engine@log_w) * em_responsibilities(st, ld)
  hit <- st@engine@nodes[, 1L] == 0 & st@engine@nodes[, 2L] == 6
  wt[, 1L] <- 0
  wt[hit, 1L] <- 1

  moved <- em_atom_step(st, ld, wt)
  expect_equal(moved@atoms[[1L]][, 1L], target, tolerance = 1e-6)
})

test_that("an atom whose step finds nothing finite stays where it is", {
  # Unreachable through the multinomial family, whose log-likelihood is finite
  # on the interior; pinned with a synthetic log-density so the contract is a
  # decision rather than an accident of seed ordering.
  st <- square_null()
  wt <- exp(st@engine@log_w) *
    em_responsibilities(st, compile_engine(st@engine))
  bottomless <- function(theta_mat) {
    matrix(-Inf, nrow = nrow(st@engine@nodes), ncol = ncol(theta_mat))
  }
  moved <- em_atom_step(st, bottomless, wt)
  expect_identical(flat_atoms(moved), flat_atoms(st))
})

test_that("an EM sweep over a triangulated part keeps its atoms in the part", {
  # The M-step optimises over the whole part, so the only boundary an atom
  # respects is the one `state@atoms` rests on: it may cross the fan's
  # diagonals freely but may not leave its part.
  st <- fw_step(square_null(), 3L)
  moved <- em_sweep(st, compile_engine(st@engine))
  for (j in seq_len(ncol(moved@atoms[[1L]]))) {
    theta <- moved@atoms[[1L]][, j]
    expect_true(contains(parts(st@null@region)[[1L]], theta, tol = 1e-5))
  }
})
