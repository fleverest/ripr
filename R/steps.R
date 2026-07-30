# Algorithm steps. Every step maps (state, problem) -> state, where `problem`
# is the list built by ripr_problem() and `state` is a mutable mixture_state
# (steps mutate it in place and return auxiliary information for the
# scheduler's records). The steps only ever touch the problem through the
# family / face / engine generics, so new sample spaces and null geometries
# need no changes here. All gradients are analytic (family score chained
# through the face parametrisation); no numerical differentiation is used.

# =============================================================================
# Objective builders
# =============================================================================

#' Frank-Wolfe objective `G(theta) = E_Q[p_theta / P_w]` as a face-oracle spec
#'
#' Builds the `value` / `grad_theta` / `value_batch` closures the face oracle
#' consumes, written purely against the engine and family generics. `log_Pw` is
#' the candidate mixture's log-density over the engine's outcome set.
#' @param engine An `expectation_engine`.
#' @param family A `family`.
#' @param log_Pw Length-M mixture log-density.
#' @return List of oracle objective closures.
#' @keywords internal
fw_objective <- function(engine, family, log_Pw) {
  list(
    value = function(theta) {
      expect_ratio(engine, eval_log_density(engine, theta), log_Pw)
    },
    grad_theta = function(theta) {
      expect_ratio_grad(
        engine,
        eval_log_density(engine, theta),
        log_Pw,
        eval_score(engine, theta)
      )
    },
    value_batch = function(theta_mat) {
      expect_ratio_batch(
        engine,
        eval_log_density_batch(engine, theta_mat),
        log_Pw
      )
    }
  )
}

#' EM M-step objective `E_Q[r_k(X) log p_theta(X)]` as a face-oracle spec
#'
#' As [fw_objective()], for the responsibility-weighted M-step: `log_r_k` is
#' atom k's log-responsibility column.
#' @param engine An `expectation_engine`.
#' @param family A `family`.
#' @param log_r_k Length-M log-responsibility column of atom k.
#' @return List of oracle objective closures.
#' @keywords internal
em_objective <- function(engine, family, log_r_k) {
  list(
    value = function(theta) {
      expect_weighted(engine, log_r_k, eval_log_density(engine, theta))
    },
    grad_theta = function(theta) {
      expect_weighted(engine, log_r_k, eval_em_score(engine, theta))
    },
    value_batch = function(theta_mat) {
      expect_weighted(engine, log_r_k, eval_log_density_batch(engine, theta_mat))
    }
  )
}

# =============================================================================
# Line searches and weight reoptimisation
# =============================================================================

#' Optimal mixing weight eps in `[1e-10, 1-1e-10]` for adding a new atom to P_w
#' @param engine An `expectation_engine`.
#' @param log_Pw Length-M mixture log-density.
#' @param log_tm_new Length-M log-density of the new atom.
#' @param tol Line-search tolerance.
#' @return Numeric mixing weight.
#' @keywords internal
fw_line_search <- function(engine, log_Pw, log_tm_new, tol = 1e-12) {
  m <- pmax(log_Pw, log_tm_new)
  exp_a <- exp(log_Pw - m)
  exp_b <- exp(log_tm_new - m)
  g <- function(eps) {
    -expect(engine, m + log((1 - eps) * exp_a + eps * exp_b))
  }
  eps_star <- optimize(g, interval = c(1e-10, 1 - 1e-10), tol = tol)$minimum
  if (g(eps_star) >= g(1e-10)) {
    eps_star <- 1e-10
  }
  eps_star
}

# Pairwise (vertex-exchange) line search: transfer alpha in [0, w_worst] from
# the worst atom to the new FW atom. Works in linear probability space since
# the update direction can be negative.
pairwise_line_search <- function(
  engine,
  log_Pw,
  log_tm_new,
  log_tm_worst,
  w_worst,
  tol = 1e-12
) {
  if (w_worst < 1e-10) {
    return(0)
  }
  Pw <- exp(log_Pw)
  p_new <- exp(log_tm_new)
  p_worst <- exp(log_tm_worst)
  delta <- p_new - p_worst
  g <- function(alpha) {
    -expect(engine, log(Pw + alpha * delta))
  }
  alpha_star <- optimize(g, interval = c(0, w_worst), tol = tol)$minimum
  if (g(alpha_star) >= g(0)) {
    alpha_star <- 0
  }
  alpha_star
}

#' Fully-corrective weight reoptimisation via entropic mirror descent
#'
#' Minimises the KL loss over the weight simplex with the atoms held fixed,
#' using exponentiated-gradient updates `w <- w * exp(-eta * grad)` (normalised)
#' with backtracking on the step size eta. `loss_and_grad(w)` must return
#' `list(loss = , grad = )`; convergence is declared when the per-iteration
#' decrease falls below `atol + rtol * abs(loss)`.
#' @param loss_and_grad Function of `w` returning `list(loss = , grad = )`.
#' @param weights Initial weight vector.
#' @param max_iters Iteration cap.
#' @param atol,rtol Convergence tolerances.
#' @return `list(weights = , loss = )`.
#' @keywords internal
mirror_descent_weights <- function(
  loss_and_grad,
  weights,
  max_iters,
  atol = 1e-12,
  rtol = 1e-8
) {
  cur <- loss_and_grad(weights)
  eta <- 1
  for (it in seq_len(max_iters)) {
    accepted <- FALSE
    while (eta >= 1e-12) {
      log_w <- log(pmax(weights, 1e-300)) - eta * cur$grad
      log_w <- log_w - max(log_w)
      w_cand <- exp(log_w)
      w_cand <- w_cand / sum(w_cand)
      cand <- loss_and_grad(w_cand)
      if (cand$loss <= cur$loss) {
        accepted <- TRUE
        break
      }
      eta <- eta / 2
    }
    if (!accepted) {
      break
    }
    decrease <- cur$loss - cand$loss
    weights <- w_cand
    cur <- cand
    eta <- eta * 2
    if (decrease < atol + rtol * abs(cur$loss)) {
      break
    }
  }
  list(weights = weights, loss = cur$loss)
}

# =============================================================================
# Steps
# =============================================================================

# Seed centres for a face's global search: the current mixture's atoms, as a
# `d x m` matrix. Each face projects them onto itself, so atoms belonging to
# other faces still contribute a (projected) centre rather than nothing. Faces
# holding no atoms of their own therefore search around the projection of the
# support -- `init_point()`'s map -- instead of an arbitrary anchor.
state_seed_centres <- function(state) {
  atoms <- state_atoms(state)[seq_len(state_n_atoms(state))]
  if (length(atoms) == 0L) {
    return(NULL)
  }
  do.call(cbind, atoms)
}

#' Frank-Wolfe oracle sweep over every null face
#'
#' Finds `argmax_theta G(theta) = E_Q[p_theta / P_W]` over the union of faces at
#' the current state. Does not modify the state. The returned `log_Pw` is the
#' mixture log-density the sweep was run against; downstream weight-update steps
#' reuse it so the whole iteration sees one consistent P_W.
#'
#' @param state A `mixture_state`.
#' @param problem Problem list from [ripr_problem()].
#' @param n_seeds Random Dirichlet seeds per face.
#' @return List with `E_star`, `best_theta`, `best_fi`, `se`, `face_results`,
#'   `log_Pw`.
#' @keywords internal
oracle_step <- function(state, problem, n_seeds) {
  log_Pw <- state_log_p_mixture(state)
  objective <- fw_objective(problem$engine, problem$family, log_Pw)
  centres <- state_seed_centres(state)
  face_results <- lapply(problem$null, function(f) {
    oracle(f, objective, n_seeds = n_seeds, seed_centres = centres)
  })
  E_ratios <- vapply(face_results, `[[`, "value", FUN.VALUE = numeric(1L))
  best_fi <- which.max(E_ratios)
  best_theta <- face_results[[best_fi]]$theta
  se <- 0
  ess <- NA_real_
  if (!deterministic(problem$engine)) {
    log_tm_best <- eval_log_density(problem$engine, best_theta)
    ratio_best <- exp(nan_to_neginf(log_tm_best - log_Pw))
    se <- expect_se(problem$engine, ratio_best)
    # Kish effective sample size of the importance ratio at the selected
    # theta*. `E_star` is a mean of `ratio_best`, so this is the number of
    # draws the estimate is *effectively* built from. It collapses when the
    # oracle picks a theta almost no draw supports -- the regime where the
    # maximum reflects sampling noise rather than the true face maximum.
    s1 <- sum(ratio_best)
    s2 <- sum(ratio_best^2)
    ess <- if (is.finite(s1) && is.finite(s2) && s2 > 0) s1^2 / s2 else 0
  }
  list(
    E_star = E_ratios[best_fi],
    best_theta = best_theta,
    best_fi = best_fi,
    se = se,
    ess = ess,
    face_results = face_results,
    log_Pw = log_Pw
  )
}

#' Vanilla Frank-Wolfe step: fixed schedule gamma = 2 / (step_index + 2)
#'
#' Adds the oracle atom with the classical step size; existing weights shrink by
#' the complementary factor.
#'
#' @param state A `mixture_state`; mutated in place.
#' @param problem Problem list.
#' @param oracle_result Result of [oracle_step()] at this state.
#' @param step_index 1-based outer iteration index driving the schedule.
#' @return List with `eps_star` and `prop_star` (always `NA` here).
#' @keywords internal
fw_step <- function(state, problem, oracle_result, step_index) {
  eps_star <- 2 / (step_index + 2)
  weights <- c(state_weights(state) * (1 - eps_star), eps_star)
  state_add_atom(state, oracle_result$best_theta, oracle_result$best_fi)
  state_set_weights(state, weights)
  list(eps_star = eps_star, prop_star = NA_real_)
}

#' Frank-Wolfe step with exact 1-D line search on the new atom's weight
#'
#' @inheritParams fw_step
#' @param ls_tol Line-search tolerance.
#' @return List with `eps_star` and `prop_star` (always `NA` here).
#' @keywords internal
fw_linesearch_step <- function(state, problem, oracle_result, ls_tol = 1e-12) {
  log_tm_new <- eval_log_density(problem$engine, oracle_result$best_theta)
  eps_star <- fw_line_search(
    problem$engine,
    oracle_result$log_Pw,
    log_tm_new,
    tol = ls_tol
  )
  weights <- c(state_weights(state) * (1 - eps_star), eps_star)
  state_add_atom(state, oracle_result$best_theta, oracle_result$best_fi)
  state_set_weights(state, weights)
  list(eps_star = eps_star, prop_star = NA_real_)
}

#' Pairwise (vertex-exchange) Frank-Wolfe step
#'
#' Moves mass from the atom with the lowest `E_Q[p_k / P_W]` to the oracle atom
#' via a 1-D line search over the transferred mass; when the worst atom's
#' remaining weight falls below `removal_thresh` it is replaced outright, so the
#' support size need not grow.
#'
#' @inheritParams fw_linesearch_step
#' @param removal_thresh Weight threshold for in-place replacement.
#' @param step_index Outer iteration index (used only in the degenerate-step
#'   warning).
#' @return List with `eps_star` (transferred mass) and `prop_star` (fraction of
#'   the worst atom's weight that moved).
#' @keywords internal
fw_pairwise_step <- function(
  state,
  problem,
  oracle_result,
  ls_tol = 1e-12,
  removal_thresh = 1e-8,
  step_index = NA_integer_
) {
  log_Pw <- oracle_result$log_Pw
  log_tm_new <- eval_log_density(problem$engine, oracle_result$best_theta)
  e_ratios <- state_atom_ratios(state, log_Pw)
  weights <- state_weights(state)
  active_idx <- which(weights > removal_thresh)
  k_worst <- active_idx[which.min(e_ratios[active_idx])]
  w_worst <- weights[k_worst]
  alpha_star <- pairwise_line_search(
    problem$engine,
    log_Pw,
    log_tm_new,
    state_atom_log_density(state, k_worst),
    w_worst,
    tol = ls_tol
  )
  if (alpha_star == 0) {
    warning(sprintf(
      "pairwise_line_search returned alpha=0 at fw_iter %d; oracle atom and worst atom may coincide.",
      step_index
    ))
  }
  if (w_worst - alpha_star < removal_thresh) {
    state_replace_atom(
      state,
      k_worst,
      oracle_result$best_theta,
      oracle_result$best_fi
    )
    weights[k_worst] <- w_worst
    state_set_weights(state, weights)
    prop_star <- 1
  } else {
    weights[k_worst] <- w_worst - alpha_star
    weights <- c(weights, alpha_star)
    state_add_atom(state, oracle_result$best_theta, oracle_result$best_fi)
    state_set_weights(state, weights)
    prop_star <- alpha_star / w_worst
  }
  list(eps_star = alpha_star, prop_star = prop_star)
}

#' Fully-corrective Frank-Wolfe step
#'
#' Seeds the new atom's weight with the 1-D line search, then reoptimises every
#' weight over the enlarged support by entropic mirror descent on the KL
#' objective.
#'
#' @inheritParams fw_linesearch_step
#' @param mirror_iters Mirror-descent iteration cap.
#' @param kl_atol,kl_rtol Convergence tolerances for the reoptimisation.
#' @return List with `eps_star` (the line-search seed) and `prop_star` (`NA`).
#' @keywords internal
fully_corrective_step <- function(
  state,
  problem,
  oracle_result,
  ls_tol = 1e-12,
  mirror_iters = 500L,
  kl_atol = 1e-12,
  kl_rtol = 1e-8
) {
  info <- fw_linesearch_step(state, problem, oracle_result, ls_tol = ls_tol)
  md <- mirror_descent_weights(
    function(w) state_objective(state, w),
    state_weights(state),
    max_iters = mirror_iters,
    atol = kl_atol,
    rtol = kl_rtol
  )
  state_set_weights(state, md$weights)
  list(eps_star = info$eps_star, prop_star = NA_real_)
}

#' Csiszar-Tusnady EM refinement on the current support
#'
#' Alternates the exact responsibility-weight update with per-atom location
#' M-steps: each atom is re-optimised on its own face (seeded from its current
#' position via [face_coordinates()]) and the move is accepted only if the
#' weighted log-likelihood improves. Runs at most `em_iters` iterations,
#' stopping early when the KL decrease falls below `kl_atol + kl_rtol * |KL|`.
#'
#' @param state A `mixture_state`; mutated in place.
#' @param problem Problem list.
#' @param em_iters Maximum EM iterations.
#' @param kl_atol,kl_rtol Convergence tolerances.
#' @param kl_init KL value at entry (used for the first convergence check).
#' @return List with `kl` and `kl_trace` (KL after each EM iteration run).
#' @keywords internal
em_step <- function(
  state,
  problem,
  em_iters,
  kl_atol,
  kl_rtol,
  kl_init
) {
  engine <- problem$engine
  faces <- problem$null
  kl <- kl_init
  n_live <- state_n_atoms(state)
  kl_trace <- numeric(0L)

  for (em_idx in seq_len(em_iters)) {
    log_Pw <- state_log_p_mixture(state)
    log_r <- state_responsibilities(state, log_Pw)

    new_weights <- colSums(nan_to_zero(exp(log_r) * engine@q_mass))
    new_weights <- pmax(new_weights, 1e-300)
    new_weights <- new_weights / sum(new_weights)

    atoms <- state_atoms(state)
    face_idx <- state_face_idx(state)
    for (k in seq_len(n_live)) {
      f_k <- faces[[face_idx[k]]]
      # Recover current atom's alpha via pseudo-inverse, clip + renormalise for
      # FP safety.
      seed_alpha_k <- face_coordinates(f_k, atoms[[k]])
      result <- oracle(
        f_k,
        em_objective(engine, problem$family, log_r[, k]),
        seed_alpha = seed_alpha_k
      )

      old_ll <- expect_weighted(
        engine,
        log_r[, k],
        state_atom_log_density(state, k)
      )
      if (result$value >= old_ll) {
        state_replace_atom(state, k, result$theta)
      }
    }

    state_set_weights(state, new_weights)
    kl_new <- state_objective(state)$loss
    kl_trace <- c(kl_trace, kl_new)
    if (kl - kl_new < kl_atol + kl_rtol * abs(kl)) {
      break
    }
    kl <- kl_new
  }
  list(kl = kl, kl_trace = kl_trace)
}

#' Li-Barron greedy step: exact-objective oracle plus line search
#'
#' Where the Frank-Wolfe oracle maximises the linearisation
#' `E_Q[p_theta / P_W]`, the Li-Barron step selects the atom maximising the
#' exact post-update objective `max_eps E_Q[log((1 - eps) P_W + eps p_theta)]`,
#' with the inner `eps` solved by the 1-D line search and the outer `theta` by
#' the face oracle. The oracle's gradient uses the envelope theorem (the inner
#' maximiser is held fixed), and its seed grid is ranked by the cheap linearised
#' objective before the exact objective drives BFGS refinement.
#'
#' Runs its own oracle sweep, so schedulers pairing it with a duality-gap record
#' still need a separate [oracle_step()].
#'
#' @param state A `mixture_state`; mutated in place.
#' @param problem Problem list.
#' @param n_seeds Random Dirichlet seeds per face.
#' @param ls_tol Line-search tolerance.
#' @return List with `eps_star`, `prop_star` (`NA`), `best_theta`, `best_fi`,
#'   and `value`.
#' @keywords internal
li_barron_step <- function(state, problem, n_seeds, ls_tol = 1e-12) {
  engine <- problem$engine
  log_Pw <- state_log_p_mixture(state)

  log_mix_at <- function(log_tm, eps) {
    m <- pmax(log_Pw, log_tm)
    m + log(exp(log_Pw - m) * (1 - eps) + exp(log_tm - m) * eps)
  }

  lb_line_search <- function(log_tm) {
    g <- function(eps) -expect(engine, log_mix_at(log_tm, eps))
    eps_star <- optimize(g, interval = c(1e-10, 1 - 1e-10), tol = ls_tol)$minimum
    if (g(eps_star) >= g(1e-10)) {
      eps_star <- 1e-10
    }
    eps_star
  }

  # Per-theta cache of the inner line search, shared between the value and
  # gradient callbacks the BFGS wrapper calls back-to-back.
  last <- new.env(parent = emptyenv())
  solve_inner <- function(theta) {
    if (!is.null(last$theta) && identical(theta, last$theta)) {
      return(last$sol)
    }
    log_tm <- eval_log_density(engine, theta)
    eps_star <- lb_line_search(log_tm)
    log_mix <- log_mix_at(log_tm, eps_star)
    sol <- list(log_tm = log_tm, eps_star = eps_star, log_mix = log_mix)
    last$theta <- theta
    last$sol <- sol
    sol
  }

  objective <- list(
    value = function(theta) {
      sol <- solve_inner(theta)
      expect(engine, sol$log_mix)
    },
    grad_theta = function(theta) {
      sol <- solve_inner(theta)
      # Envelope theorem: eps* is optimal, so only the explicit theta dependence
      # contributes.
      sol$eps_star *
        expect_ratio_grad(
          engine,
          sol$log_tm,
          sol$log_mix,
          eval_score(engine, theta)
        )
    },
    value_batch = function(theta_mat) {
      # Seed ranking by the linearised objective: a cheap proxy that orders
      # candidate atoms almost identically for seed selection.
      expect_ratio_batch(
        engine,
        eval_log_density_batch(engine, theta_mat),
        log_Pw
      )
    }
  )

  centres <- state_seed_centres(state)
  face_results <- lapply(problem$null, function(f) {
    oracle(f, objective, n_seeds = n_seeds, seed_centres = centres)
  })
  values <- vapply(face_results, `[[`, "value", FUN.VALUE = numeric(1L))
  best_fi <- which.max(values)
  best_theta <- face_results[[best_fi]]$theta

  log_tm_new <- eval_log_density(engine, best_theta)
  eps_star <- lb_line_search(log_tm_new)
  weights <- c(state_weights(state) * (1 - eps_star), eps_star)
  state_add_atom(state, best_theta, best_fi)
  state_set_weights(state, weights)
  list(
    eps_star = eps_star,
    prop_star = NA_real_,
    best_theta = best_theta,
    best_fi = best_fi,
    value = values[best_fi]
  )
}
