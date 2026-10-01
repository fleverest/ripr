# Properties of R/steps.R, observed through the verbs.
#
# The identity `sum_c w_c G(theta_c) = 1` is the workhorse. It is algebraic, not
# a convergence condition, so it holds at every iterate and fails only if the
# quadrature weights are unnormalised, the atoms and weights are misaligned, or
# `log_p` was computed against the wrong mixture.
#
# Every quantity is checked against the brute-force references in
# helper-fit.R, which enumerate the sample space through `log_density()` and
# share nothing with the step layer. So these tests say what a step must do,
# not how, and a refactor of the internals that keeps the answers leaves them
# alone. test-fit.R covers the verbs' bookkeeping; this file covers their
# arithmetic. The few numerical kernels no verb exposes are pinned directly at
# the foot of the file.

plurality <- function(
  n = 12,
  k = 4,
  q = c(0.42, 0.31, 0.16, 0.11),
  atoms = NULL,
  weights = NULL,
  record_gap = FALSE,
  ...
) {
  fam <- multinomial_family(n_trials = n, k = k)
  Q <- mixture(fam, dirac(theta = q))
  parts <- lapply(2:k, function(j) {
    basis <- lapply(setdiff(seq_len(k), 1L), function(i) {
      replace(numeric(k), i, 1)
    })
    tie <- replace(numeric(k), c(1L, j), 0.5)
    simplex_region(vertices = do.call(rbind, c(basis, list(tie))))
  })
  ripr_init(
    Q,
    null_model(fam, parts),
    exact_engine(),
    atoms = atoms,
    weights = weights,
    record_gap = record_gap,
    control = ripr_control(n_seeds = 30L, n_restarts = 4L, ...)
  )
}

# One atom, in the first part, so nothing is left to move mass away from.
lone_atom <- function(...) {
  empty <- matrix(numeric(0), nrow = 0L, ncol = 4L)
  plurality(
    atoms = list(rbind(c(0.2, 0.4, 0.2, 0.2)), empty, empty),
    ...
  )
}

# --- G, as the weight verb sees it --------------------------------------------

test_that("a weight step multiplies each weight by its atom's G", {
  # `w_c <- w_c G(theta_c)` is the exact M-step for the weights, and the
  # identity makes it self-normalising. Checked at weights that are nowhere
  # near optimal, since the identity is algebra and does not care.
  st <- fw_step(plurality(), 3L)
  a <- atoms(st@mixing)
  n <- nrow(a)
  by_part <- split(seq_len(n), factor(st@part, levels = 1:3))
  ord <- unlist(by_part, use.names = FALSE)
  for (w in list(
    weights(st@mixing),
    rev(weights(st@mixing)),
    c(0.7, rep(0.3 / (n - 1), n - 1))
  )) {
    # Rebuilt part by part, as `ripr_init()` takes them.
    fresh <- plurality(
      atoms = lapply(by_part, \(i) a[i, , drop = FALSE]),
      weights = lapply(by_part, \(i) w[i])
    )
    G <- exact_g(fresh, atoms(fresh@mixing))
    expect_equal(sum(w[ord] * G), 1, tolerance = rounding_tol(1))
    swept <- weight_step(fresh, 1L)
    expect_equal(weights(swept@mixing), w[ord] * G, tolerance = rounding_tol(1))
  }
})

test_that("a weight row's oracle value is the restricted Frank-Wolfe gap", {
  # `max_c G(theta_c)`, measured before the sweep: what the sweep had left to
  # gain over the current support.
  st <- fw_step(plurality(), 3L)
  G <- exact_g(st, atoms(st@mixing))
  swept <- weight_step(st, 1L)
  expect_equal(utils::tail(swept@trace$oracle_value, 1L), max(G))
  # And the predicate reading the same number agrees on either side of it.
  expect_true(support_gap_below(max(G) - 1 + 1e-9)(st))
  expect_false(support_gap_below(max(G) - 1 - 1e-9)(st))
})

test_that("weight steps are monotone and record the KL they reach", {
  st <- fw_step(plurality(), 3L)
  swept <- weight_step(st, 5L)
  kl <- swept@trace$kl
  expect_true(all(diff(kl) <= rounding_tol(kl[-1L])))
  expect_equal(
    utils::tail(kl, 1L),
    exact_kl(swept),
    tolerance = rounding_tol(1)
  )
})

# --- Step paths ---------------------------------------------------------------

test_that("every variant stays on the simplex and records its KL exactly", {
  # The forward path takes a two-column shortcut for the mixture rather than
  # rebuilding it, and every line search reads KL off that path. So each
  # snapshot's KL, recomputed from scratch, must match the row that made it.
  # And a line search cannot increase KL: `gamma = 0` is in every interval,
  # so the step can always decline to move.
  for (v in c("standard", "away-step", "pairwise")) {
    set.seed(3)
    st <- fw_step(plurality(snapshot = "all"), 6L, variant = v)
    kl <- st@trace$kl
    expect_true(all(diff(kl) <= rounding_tol(kl[-1L])))
    for (snap in st@snapshots) {
      w <- weights(snap$mixing)
      expect_equal(sum(w), 1, tolerance = rounding_tol(1))
      expect_true(all(w > 0))
      expect_equal(
        st@trace$kl[snap$step + 1L],
        exact_kl(st, snap$mixing),
        tolerance = rounding_tol(1)
      )
    }
  }
})

test_that("pairwise and away-step reduce to forward with one active atom", {
  # Pairwise then has only that atom to take mass from, so both paths are
  # `[1 - gamma, gamma]` with a cap of 1; away has nothing to move mass to and
  # is dropped from the choice. So the one-atom case needs no special handling.
  st <- lone_atom()
  set.seed(5)
  fwd <- fw_step(st, 1L)
  for (v in c("pairwise", "away-step")) {
    set.seed(5)
    other <- fw_step(st, 1L, variant = v)
    expect_equal(weights(other@mixing), weights(fwd@mixing))
    expect_equal(atoms(other@mixing), atoms(fwd@mixing))
    # The two line searches minimise the same KL along differently computed
    # paths; KL is flat at its minimum, so its argmin is only determined to
    # about sqrt(machine epsilon).
    expect_equal(
      utils::tail(other@trace$step_size, 1L),
      utils::tail(fwd@trace$step_size, 1L),
      tolerance = 1e-6
    )
  }
})

test_that("a drop step removes its atom from the support", {
  # Only a removal can make the support fall: every other step either adds an
  # atom or leaves the count alone. Under this seed the away cap's residual
  # would otherwise strand an atom at ~1e-18, and a line search stopping a
  # whisker short of the cap would leave one at ~1e-13.
  set.seed(2)
  st <- fw_step(plurality(snapshot = "all"), 18L, variant = "away-step")
  fw <- st@trace[st@trace$phase == "fw", ]
  expect_true(any(diff(fw$support_size) < 0))
  w <- weights(st@mixing)
  expect_true(all(w > 0))
  expect_false(any(w < 1e-9))
  expect_identical(nrow(atoms(st@mixing)), length(w))

  # The atom goes with its part: what survives a drop is a subset of what was
  # there, each atom still filed where it was.
  drops <- which(diff(st@trace$support_size) < 0)
  for (i in drops) {
    before <- st@snapshots[[i]]
    after <- st@snapshots[[i + 1L]]
    kept <- match(
      apply(atoms(after$mixing), 1L, paste, collapse = ","),
      apply(atoms(before$mixing), 1L, paste, collapse = ",")
    )
    expect_false(anyNA(kept))
    expect_identical(after$part, before$part[kept])
  }
})

test_that("repeated away steps do not let the weights drift off the simplex", {
  # Two atoms on top of each other at the tie, with the projection already
  # reached: the away step then shuffles mass between them and nothing else.
  # `(1 + gamma) w - gamma e_v` sums to one exactly only in exact arithmetic,
  # so without renormalising the total drifts an ulp per step.
  fam <- multinomial_family(n_trials = 20, k = 2)
  Q <- mixture(fam, dirac(theta = c(0.6, 0.4)))
  v <- diag(2)
  v[1L, ] <- c(0.5, 0.5)
  tie <- c(0.5, 0.5)
  st <- ripr_init(
    Q,
    null_model(fam, list(simplex_region(vertices = v))),
    atoms = list(rbind(tie, tie, deparse.level = 0)),
    weights = list(c(0.5, 0.5))
  )
  expect_no_error(stepped <- fw_step(st, 60L, variant = "away-step"))
  expect_equal(sum(weights(stepped@mixing)), 1, tolerance = rounding_tol(1))
})

test_that("the fixed schedule is capped at the path's own maximum", {
  # Pairwise caps at the worst atom's weight and away below 1; an uncapped
  # schedule value would take the weights off the simplex. Under this seed
  # the cap binds, so a step shorter than the schedule is recorded, and at the
  # cap the worst atom must leave exactly, not as a sliver of weight.
  for (v in c("pairwise", "away-step")) {
    set.seed(1)
    st <- plurality()
    # Long enough for the cap to bind under each.
    stepped <- fw_step(
      st,
      if (v == "pairwise") 4L else 7L,
      variant = v,
      size = "fixed"
    )
    fw <- stepped@trace[stepped@trace$phase == "fw", ]
    k <- st@trace$support_size + seq_len(nrow(fw)) - 1L
    expect_true(all(fw$step_size <= 2 / (k + 2) + rounding_tol(1)))
    expect_true(any(fw$step_size < 2 / (k + 2) - 1e-6))
    w <- weights(stepped@mixing)
    expect_equal(sum(w), 1, tolerance = rounding_tol(1))
    expect_false(any(w < 1e-9))
  }
})

test_that("a step's part says whether it used the candidate", {
  # Derived from the weight the candidate ends with: an atom that entered is
  # the oracle's point, filed under the recorded part; a row with no part
  # added nothing, so its support cannot have grown.
  set.seed(6)
  st <- fw_step(plurality(snapshot = "all"), 10L, variant = "away-step")
  tr <- st@trace
  for (i in which(tr$phase == "fw")) {
    if (is.na(tr$part[i])) {
      expect_lte(tr$support_size[i], tr$support_size[i - 1L])
    } else {
      after <- st@snapshots[[i]]
      hit <- apply(atoms(after$mixing), 1L, identical, tr$oracle_theta[[i]])
      expect_true(tr$part[i] %in% after$part[hit])
    }
  }
  expect_true(any(is.na(tr$part[tr$phase == "fw"])))
  expect_true(all(is.na(tr$part[tr$direction %in% "away"])))
})

# --- The corrective solve -----------------------------------------------------

test_that("the finishing weight solve respects both its budget and tolerance", {
  # `fc_tol` and `fc_max_iter` are read only by a corrective solve, so they
  # change nothing about the plain steps taken to get here; a looser tolerance
  # or a smaller budget can only stop the solve earlier, never beat it.
  full <- ripr_finish(fw_step(plurality(), 4L), reoptimise = TRUE)
  loose <- ripr_finish(fw_step(plurality(fc_tol = 1e10), 4L), reoptimise = TRUE)
  capped <- ripr_finish(fw_step(plurality(fc_max_iter = 2L), 4L), reoptimise = TRUE)
  expect_gte(loose@kl, full@kl - 1e-12)
  expect_gte(capped@kl, full@kl - 1e-12)
  expect_gt(capped@kl, full@kl)
})

# --- Oracles ------------------------------------------------------------------

test_that("the linear oracle's value is G at the point it reports", {
  # And at least G at every atom, since the atoms seed the search and the
  # identity makes their maximum at least 1.
  st <- fw_step(plurality(), 2L, record_gap = TRUE)
  expect_equal(st@oracle$value, exact_g(st, st@oracle$theta))
  expect_gte(
    st@oracle$value,
    max(exact_g(st, atoms(st@mixing))) - rounding_tol(1)
  )
  expect_true(contains(
    parts(st@null@region)[[st@oracle$part]],
    st@oracle$theta,
    tol = 1e-6
  ))
})

test_that("the linear oracle finds the maximum of G over a small null", {
  # A wrong gradient or batch evaluator does not break the search outright; it
  # leaves it short of the maximum. On a two-part K = 3 null a grid is dense
  # enough to know the maximum, so the search must reach it.
  fam <- multinomial_family(n_trials = 6L, k = 3L)
  parts <- lapply(2:3, function(j) {
    simplex_region(
      vertices = rbind(
        c(0, 1, 0),
        c(0, 0, 1),
        replace(numeric(3), c(1L, j), 0.5)
      )
    )
  })
  set.seed(7)
  st <- ripr_init(
    mixture(fam, dirac(c(0.5, 0.3, 0.2))),
    null_model(fam, parts),
    record_gap = TRUE,
    control = ripr_control(n_seeds = 30L, n_restarts = 4L)
  )
  grid <- as.matrix(expand.grid(
    a = seq(0, 0.5, by = 0.025),
    b = seq(0, 1, by = 0.025)
  ))
  grid <- cbind(grid, 1 - rowSums(grid))
  in_null <- apply(grid, 1L, \(t) all(t >= 0) && contains(st@null@region, t))
  grid_max <- max(exact_g(st, grid[in_null, , drop = FALSE]))
  expect_gte(st@oracle$value, grid_max - 1e-9)
})

test_that("a Li-Barron row records the KL its step leaves as the objective", {
  # The Li-Barron oracle scores a candidate by `-KL` after its step, so the
  # value it attains is exactly the negated KL of the row it produced; and
  # since the forward step towards the linear oracle's point is one of the
  # candidates it considers, its step does at least as well as that one.
  st <- fw_step(plurality(), 2L)
  set.seed(8)
  lb <- lb_step(st, 1L)
  row <- utils::tail(lb@trace, 1L)
  expect_equal(row$oracle_value, -row$kl)
  expect_equal(row$kl, exact_kl(lb), tolerance = rounding_tol(1))
  set.seed(8)
  fw <- fw_step(st, 1L)
  expect_lte(row$kl, utils::tail(fw@trace$kl, 1L) + rounding_tol(1))
})

# --- Support identification ---------------------------------------------------

test_that("identification does not increase KL", {
  # It compares two feasible points and keeps the better.
  st <- fw_step(plurality(), 6L)
  plain <- ripr_finish(st)
  identified <- ripr_finish(st, identify = TRUE)
  expect_lte(identified@kl, plain@kl + rounding_tol(plain@kl))
  expect_equal(
    identified@kl,
    exact_kl(st, identified@W0),
    tolerance = rounding_tol(1)
  )
})

test_that("identification survives an atom carrying nearly all the mass", {
  # Removing the dominant atom divides by 1e-15 and puts the removal share at
  # 1 to machine precision; the pass must still come through with usable
  # weights.
  a <- atoms(plurality()@mixing)
  heavy <- c(1 - 1e-15, 5e-16, 5e-16)
  st <- plurality(
    atoms = lapply(1:3, \(i) a[i, , drop = FALSE]),
    weights = as.list(heavy)
  )
  fit <- ripr_finish(st, identify = TRUE)
  expect_false(anyNA(weights(fit@W0)))
  expect_equal(sum(weights(fit@W0)), 1, tolerance = rounding_tol(1))
  expect_false(is.nan(fit@kl))
})

# --- EM -----------------------------------------------------------------------

test_that("an EM sweep keeps every atom in its own part", {
  # The M-step optimises each component over its own chart, so an atom cannot
  # migrate between parts however the responsibilities fall.
  st <- fw_step(plurality(), 2L)
  expect_true(in_own_part(em_step(st, 1L)))
})

test_that("an EM sweep neither grows nor shrinks the support", {
  st <- fw_step(plurality(), 2L)
  expect_identical(em_step(st, 1L)@part, st@part)
})

test_that("an EM sweep's weights are the weight step's", {
  # The M-step for the weights, `E_Q[r_c]`, is the multiplicative update
  # `w_c G(theta_c)` reached from the other direction, and the atom M-step
  # that follows leaves the weights alone.
  st <- fw_step(plurality(), 2L)
  expect_equal(
    weights(em_step(st, 1L)@mixing),
    weights(weight_step(st, 1L)@mixing),
    tolerance = rounding_tol(1)
  )
})

test_that("an EM sweep records the KL of the mixture it made", {
  st <- em_step(fw_step(plurality(), 2L), 1L)
  expect_equal(
    utils::tail(st@trace$kl, 1L),
    exact_kl(st),
    tolerance = rounding_tol(1)
  )
})

test_that("the atom M-step is deterministic", {
  # No random seeding: the same state must give the same moved atoms whatever
  # the RNG is doing, which is what makes a fit reproducible without a seed.
  st <- plurality()
  set.seed(41)
  a <- em_step(st, 1L)
  set.seed(42)
  b <- em_step(st, 1L)
  expect_identical(atoms(a@mixing), atoms(b@mixing))
})

# --- An atom on the boundary of the support -----------------------------------

test_that("the EM M-step survives a node its atom gives zero probability", {
  # `0 * -Inf` is `NaN`, not the 0 the term is worth under `0 log 0 = 0`: a node
  # an atom cannot produce has both zero responsibility and `-Inf` log density,
  # for the same reason. Left unhandled it reaches the optimiser, which
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
      v[1L, ] <- replace(numeric(3L), c(1L, j), 0.5)
      simplex_region(vertices = v)
    })
  )
  # One interior atom keeps KL finite; the second sits on the facet
  # `theta_3 = 0` of its part, so it gives zero probability to every outcome
  # using the third category.
  facet <- c(0.5, 0.5, 0)
  st <- ripr_init(
    fam(c(0.40, 0.35, 0.25)),
    null,
    atoms = list(
      rbind(c(1 / 3, 1 / 3, 1 / 3), facet),
      rbind(c(1 / 3, 1 / 3, 1 / 3))
    )
  )

  # The state really does contain such an atom, at outcomes Q can produce,
  # or this guards nothing.
  x <- enumerate_space(fam@sample_space)
  expect_true(any(
    log_density(st@alternative, x) > -Inf & log_density(fam(facet), x) == -Inf
  ))

  after <- em_step(st, times = 5L)
  kl <- after@trace$kl
  expect_false(anyNA(kl))
  expect_lte(utils::tail(kl, 1L), kl[1L])
})

# --- Cells under the steps ----------------------------------------------------

# A null whose one part is a convex hull rather than a simplex: the square
# `{theta_1 <= 1/2, theta_2 <= 1/2}` in the 2-simplex, which triangulates into
# two cells. Every part of every other null in this file is its own only cell,
# so this is where the cell wiring has anything to do.
square <- function() {
  polytope_region(
    vertices = rbind(
      c(0.5, 0.5, 0),
      c(0, 0.5, 0.5),
      c(0, 0, 1),
      c(0.5, 0, 0.5)
    )
  )
}

square_null <- function(
  n = 12,
  q = c(0.15, 0.35, 0.5),
  atoms = NULL,
  record_gap = FALSE,
  ...
) {
  fam <- multinomial_family(n_trials = n, k = 3L)
  ripr_init(
    mixture(fam, dirac(theta = q)),
    null_model(fam, list(square())),
    exact_engine(),
    atoms = atoms,
    record_gap = record_gap,
    control = ripr_control(n_seeds = 30L, n_restarts = 4L, ...)
  )
}

test_that("the oracle searches cells and reports the part they belong to", {
  st <- square_null(record_gap = TRUE)
  expect_length(st@null@cells, 2L)
  found <- st@oracle

  # The index is a part's, because that is what an atom is filed under in
  # `state@part`, and the trace files the gap the same way.
  expect_identical(found$part, 1L)
  expect_identical(st@trace$gap_after_part, 1L)
  expect_true(contains(parts(st@null@region)[[found$part]], found$theta))

  # And the maximiser is in one of the cells searched, not merely in the hull.
  expect_true(any(vapply(
    st@null@cells,
    \(cell) contains(cell, found$theta, tol = 1e-6),
    logical(1)
  )))
})

test_that("the atom M-step crosses a fan diagonal when the optimum does", {
  # With a lone atom its responsibility is 1 everywhere, so the M-step is the
  # maximum likelihood projection of Q onto the part: for Q at (0, 0.5, 0.5),
  # Q itself, a vertex of exactly one of the square's two cells. An atom
  # started deep inside the other cell must cross the fan's diagonal to reach
  # it -- which confinement to a containing cell forbade.
  target <- c(0, 0.5, 0.5)
  fam <- multinomial_family(n_trials = 12, k = 3L)
  cells <- null_model(fam, list(square()))@cells
  holds <- vapply(
    cells,
    \(cell) contains(cell, target, tol = rounding_tol(1)),
    logical(1)
  )
  expect_identical(sum(holds), 1L)
  far <- colMeans(cells[[which(!holds)[1L]]]@vertices)

  st <- square_null(q = target, atoms = list(matrix(far, nrow = 1L)))
  moved <- em_step(st, 1L)
  expect_equal(atoms(moved@mixing)[1L, ], target, tolerance = 1e-6)
})

test_that("an EM sweep over a triangulated part keeps its atoms in the part", {
  # The M-step optimises over the whole part, so the only boundary an atom
  # respects is the one `state@part` rests on: it may cross the fan's
  # diagonals freely but may not leave its part.
  st <- fw_step(square_null(), 3L)
  moved <- em_step(st, 1L)
  for (j in seq_len(n_atoms(moved@mixing))) {
    expect_true(contains(square(), atoms(moved@mixing)[j, ], tol = 1e-5))
  }
})

# --- Internal kernels ---------------------------------------------------------
#
# Numerical pieces whose correctness no verb can show on its own. A wrong
# gradient or a sliver of weight left at a cap does not fail a fit; it makes
# the search weaker or the support untidy by amounts too small or too
# seed-dependent to assert on through `fw_step()`. These call the internals by
# name, so a refactor that renames them must update this section and nothing
# above it.

test_that("pairwise and away empty the worst atom exactly at their cap", {
  # Exactly zero, not merely small: the step drops atoms whose weight is zero,
  # and a residual would strand one. For away, `(1 + gamma) w - gamma` at
  # `gamma = w / (1 - w)` leaves a positive residual for `w = 1/3`. No verb
  # can be steered to a given weight at a cap, so the maps are asked directly;
  # neither reads the densities to compute its weights.
  w <- c(1 / 3, 1 / 3, 1 / 3, 0)
  away <- path_away(NULL, w, NULL, 1L, 0.5)
  expect_identical(away$w_of(away$gamma_max)[1L], 0)
  w <- c(0.1, 0.3, 0.6, 0)
  pairwise <- path_pairwise(NULL, w, NULL, 1L, 0.5)
  expect_identical(pairwise$w_of(pairwise$gamma_max)[1L], 0)
})

test_that("removing an atom adjusts the mixture in place", {
  # Identification tests every candidate atom, so rebuilding the mixture each
  # time would make a pass O(MC^2). `log_p_without` adjusts the existing one
  # instead, in O(M):
  #
  #   log P + log(1 - w_j p_j / P) - log(1 - w_j)
  #     =  log((P - w_j p_j) / (1 - w_j))
  #
  # which is not hard to get wrong. A wrong downdate only mis-ranks removals,
  # which `ripr_finish(identify = TRUE)` would hide behind its own KL check,
  # so it is checked against the definition here, including where one atom
  # carries all or nearly all the mass.
  set.seed(9)
  p <- matrix(stats::rexp(20), nrow = 5L)
  p <- sweep(p, 2L, colSums(p), "/")
  mix <- function(w) log(as.vector(p %*% w))
  for (w in list(c(0.4, 0.3, 0.2, 0.1), c(1 - 1e-15, rep(1e-15 / 3, 3)))) {
    dropped <- replace(w, 2L, 0)
    expect_equal(
      log_p_without(mix(w), log(p[, 2L]), w[2L]),
      mix(dropped / sum(dropped))
    )
  }
  heavy <- c(1 - 1e-15, rep(1e-15 / 3, 3))
  expect_false(anyNA(log_p_without(mix(heavy), log(p[, 1L]), heavy[1L])))
  # At a weight of exactly 1 nothing is left: finite arithmetic, all `-Inf`.
  out <- log_p_without(log(p[, 1L]), log(p[, 1L]), 1)
  expect_false(anyNA(out))
  expect_true(all(out == -Inf))
})

test_that("the oracles' gradients and batch forms agree with their values", {
  # The search uses all three, but a wrong gradient only slows or misdirects
  # it, and the multistart can paper over that on any one problem.
  st <- plurality()
  ld <- compile_engine(st@engine)
  log_p <- log_p_at_nodes(st, ld)
  theta <- c(0.30, 0.40, 0.20, 0.10)
  # A feasible direction on the simplex: components must sum to zero.
  d <- c(1, -1, 0, 0) * 1e-6

  lin <- linear_oracle(st, log_p, ld)
  a <- atoms(st@mixing)
  expect_equal(lin$value_batch(a), exact_g(st, a))
  expect_equal(
    vapply(seq_len(nrow(a)), \(i) lin$value(a[i, ]), numeric(1)),
    lin$value_batch(a)
  )
  fd <- (lin$value(theta + d) - lin$value(theta - d)) / 2
  expect_equal(sum(lin$grad(theta) * d), fd, tolerance = 1e-6)

  # Under `size = "fixed"` the Li-Barron step size does not depend on theta,
  # so the objective is smooth. Under a line search it is flat wherever
  # G <= 1, since the search then declines to move at all: not a defect, but
  # why the search depends on some seed landing where G > 1.
  lb <- nonlinear_oracle(st, log_p, ld, size = "fixed", gamma_fixed = 0.4)
  fd <- (lb$value(theta + d) - lb$value(theta - d)) / 2
  expect_equal(sum(lb$grad(theta) * d), fd, tolerance = 1e-6)

  flat <- nonlinear_oracle(st, log_p, ld)
  poor <- Filter(
    function(th) exact_g(st, th) < 1,
    lapply(c(0.05, 0.10, 0.15), function(t1) c(t1, (1 - t1) * c(0.5, 0.3, 0.2)))
  )
  expect_true(length(poor) >= 2L)
  values <- vapply(poor, flat$value, numeric(1))
  expect_equal(values, rep(-utils::tail(st@trace$kl, 1L), length(values)))
})

test_that("an atom whose M-step finds nothing finite stays where it is", {
  # Unreachable through the multinomial family, whose log-likelihood is finite
  # on the interior; pinned with a synthetic log-density so the contract is a
  # decision rather than an accident of seed ordering.
  st <- square_null()
  wt <- matrix(1, nrow = nrow(st@engine@nodes), ncol = n_atoms(st@mixing))
  bottomless <- function(theta_mat) {
    matrix(-Inf, nrow = nrow(st@engine@nodes), ncol = nrow(theta_mat))
  }
  moved <- em_atom_step(st, bottomless, wt)
  expect_identical(atoms(moved@mixing), atoms(st@mixing))
})
