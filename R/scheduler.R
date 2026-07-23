#' Compute the reverse information projection (RIPr)
#'
#' The sole public entry point. Runs the reference hybrid schedule on a
#' [ripr_problem()]: optional initial EM, then `fw_iters` outer iterations of
#' (weight-update step, EM refinement, oracle sweep), with per-iteration
#' recording. All variant behaviour comes from which weight-update step
#' `fw_variant` selects; the loop itself is variant-agnostic.
#'
#' The result carries the atoms/weights of the final mixture, the KL and gap
#' traces, a `metrics` data.frame (one row per outer iteration), and the
#' terminal [certify()] result at the final state.
#'
#' @param problem Problem list from [ripr_problem()].
#' @param init_atoms `d x l` matrix of initial atom locations on the null
#'   faces. Each column is one atom.
#' @param init_atom_faces Integer vector of face indices, one per column.
#' @param fw_iters Number of Frank-Wolfe iterations. Set to `0L` for pure EM.
#' @param em_iters Max EM iterations per outer step. Set to `0L` to skip EM.
#' @param n_seeds Random Dirichlet/Gaussian seeds per face for the oracle.
#' @param kl_atol,kl_rtol Absolute/relative tolerances for EM convergence.
#' @param gap_tol Convergence tolerance on the FW gap (`E_star - 1`).
#' @param ls_tol Line-search tolerance.
#' @param removal_thresh Weight threshold below which the worst atom is replaced
#'   in place (pairwise variant only).
#' @param mirror_iters Max entropic-mirror-descent iterations per weight
#'   reoptimisation (fully-corrective variant only).
#' @param prune_threshold Drop projection atoms with weight `<= prune_threshold`
#'   (renormalising the rest) before the final certification, so the returned
#'   e-variable and certificate describe one clean mixture. Default `0` keeps
#'   every atom.
#' @param fw_variant Weight-update variant: `"line-search"` (default),
#'   `"pairwise"`, `"vanilla"`, `"fully-corrective"`, or `"li-barron"` (the
#'   exact-objective greedy step).
#' @param checkpoint_iters Integer vector of outer-iteration indices (0 =
#'   post-init) at which to snapshot the mixture. Default `NULL` records nothing.
#' @param verbose Print per-iteration progress. Default TRUE.
#' @return A list. Its principal element is `e_variable`, the fitted RIPr
#'   [e_variable()] (numerator `Q`, projection `P*`, and the certified gap); call
#'   [e_value()] on it to score data. Also returned: the projection's `atoms` /
#'   `weights` / `atom_face_idx`, the `certificate`, the `kl_trace` / `history` /
#'   `metrics` records, the terminal `gap`, `kl`, `converged` flag, and the raw
#'   optimisation `state`.
#' @export
run_ripr <- function(
  problem,
  init_atoms,
  init_atom_faces,
  fw_iters,
  em_iters,
  n_seeds = 200L,
  kl_atol = 1e-12,
  kl_rtol = 1e-8,
  gap_tol = 1e-8,
  ls_tol = 1e-12,
  removal_thresh = 1e-8,
  mirror_iters = 500L,
  prune_threshold = 0,
  fw_variant = c(
    "line-search",
    "pairwise",
    "vanilla",
    "fully-corrective",
    "li-barron"
  ),
  checkpoint_iters = NULL,
  verbose = TRUE
) {
  n_init <- ncol(init_atoms)
  if (length(init_atom_faces) != n_init) {
    stop(sprintf(
      "length(init_atom_faces) (%d) must equal ncol(init_atoms) (%d).",
      length(init_atom_faces),
      n_init
    ))
  }

  state <- mixture_state(problem$engine, n_init + fw_iters)
  t_start <- proc.time()[["elapsed"]]

  # --- History accumulators ---
  trace_rows <- list()
  outer_rows <- list()
  support_sizes <- integer(0L)
  gap_ses <- numeric(0L)
  checkpoints <- list()
  kl_ulb <- -Inf
  best_gr <- -Inf

  record_trace <- function(iter, type, n_atoms, kl_val) {
    trace_rows[[length(trace_rows) + 1L]] <<- data.frame(
      iter = iter,
      step_type = type,
      n_atoms = n_atoms,
      kl = kl_val
    )
  }

  record_checkpoint <- function(iter, oracle_theta_cp = NULL) {
    if (is.null(checkpoint_iters) || !(iter %in% checkpoint_iters)) {
      return(invisible(NULL))
    }
    snap <- state_snapshot(state)
    checkpoints[[length(checkpoints) + 1L]] <<- list(
      iter = iter,
      atoms = snap$atoms,
      weights = snap$weights,
      atom_face_idx = snap$atom_face_idx,
      oracle_theta = oracle_theta_cp
    )
  }

  fw_variant <- match.arg(fw_variant)
  fw_mode <- if (fw_iters == 0L) "em-only" else fw_variant

  if (verbose) {
    message(sprintf(
      "run_ripr: K=%d, M=%d outcomes, n_init=%d, fw_iters=%d, mode=%s, em_iters=%d, kl_atol=%g, kl_rtol=%g, gap_tol=%g, n_seeds=%d",
      param_dim(problem$family),
      n_outcomes(problem$engine),
      n_init,
      fw_iters,
      fw_mode,
      em_iters,
      kl_atol,
      kl_rtol,
      gap_tol,
      n_seeds
    ))
  }

  # --- Initialisation ---
  for (j in seq_len(n_init)) {
    state_add_atom(state, init_atoms[, j], init_atom_faces[j])
  }
  state_set_weights(state, rep(1 / n_init, n_init))
  kl <- state_objective(state)$loss
  record_trace(0L, "init", state_n_atoms(state), kl)

  # --- Initial EM ---
  if (em_iters > 0L) {
    em_result <- em_step(state, problem, em_iters, kl_atol, kl_rtol, kl)
    for (kl_em in em_result$kl_trace) {
      record_trace(0L, "em", state_n_atoms(state), kl_em)
    }
    kl <- em_result$kl
  }
  # --- Initial oracle ---
  converged <- FALSE
  fw <- oracle_step(state, problem, n_seeds = n_seeds)
  gap <- fw$E_star - 1
  oracle_theta <- fw$best_theta
  kl_ulb <- kl - gap
  gr <- kl - log1p(gap)
  best_gr <- max(best_gr, gr)
  outer_rows[[1L]] <- data.frame(
    iter = 0L,
    face_idx = fw$best_fi,
    gap = gap,
    eps_star = NA_real_,
    prop_star = NA_real_,
    kl_after_fw = NA_real_,
    kl_after_em = kl,
    kl_ulb = kl_ulb,
    gr = gr,
    best_gr = best_gr,
    elapsed_s = proc.time()[["elapsed"]] - t_start
  )
  support_sizes <- c(support_sizes, state_n_atoms(state))
  gap_ses <- c(gap_ses, fw$se %||% 0)
  if (verbose) {
    message(sprintf(
      "Init [%d atoms]: Gap %e, KL %e, ULB %e, GR %e",
      state_n_atoms(state),
      gap,
      kl,
      kl_ulb,
      gr
    ))
  }
  record_checkpoint(0L, oracle_theta_cp = fw$best_theta)
  if (gap < gap_tol) {
    converged <- TRUE
  }

  # --- FW loop ---
  kl_prev <- kl
  for (fw_idx in seq_len(fw_iters)) {
    if (converged) {
      break
    }

    # Weight-update step -- uses fw from the oracle at end of previous phase.
    step_info <- switch(
      fw_variant,
      "pairwise" = fw_pairwise_step(
        state,
        problem,
        fw,
        ls_tol = ls_tol,
        removal_thresh = removal_thresh,
        step_index = fw_idx
      ),
      "line-search" = fw_linesearch_step(
        state,
        problem,
        fw,
        ls_tol = ls_tol
      ),
      "fully-corrective" = fully_corrective_step(
        state,
        problem,
        fw,
        ls_tol = ls_tol,
        mirror_iters = mirror_iters,
        kl_atol = kl_atol,
        kl_rtol = kl_rtol
      ),
      "vanilla" = fw_step(state, problem, fw, step_index = fw_idx),
      "li-barron" = li_barron_step(
        state,
        problem,
        n_seeds = n_seeds,
        ls_tol = ls_tol
      )
    )
    eps_star <- step_info$eps_star
    prop_star <- step_info$prop_star

    kl_after_fw <- state_objective(state)$loss
    record_trace(fw_idx, "fw", state_n_atoms(state), kl_after_fw)

    # EM refinement
    if (em_iters > 0L) {
      em_result <- em_step(
        state,
        problem,
        em_iters,
        kl_atol,
        kl_rtol,
        kl_after_fw
      )
      for (kl_em in em_result$kl_trace) {
        record_trace(fw_idx, "em", state_n_atoms(state), kl_em)
      }
      kl <- em_result$kl
    } else {
      kl <- kl_after_fw
    }
    # Oracle -- single call per iteration; result used for next iteration's step.
    fw <- oracle_step(state, problem, n_seeds = n_seeds)
    gap <- fw$E_star - 1
    oracle_theta <- fw$best_theta
    kl_ulb <- max(kl_ulb, kl - gap)
    gr <- kl - log1p(gap)
    best_gr <- max(best_gr, gr)
    outer_rows[[length(outer_rows) + 1L]] <- data.frame(
      iter = fw_idx,
      face_idx = fw$best_fi,
      gap = gap,
      eps_star = eps_star,
      prop_star = if (is.na(prop_star)) NA_real_ else prop_star,
      kl_after_fw = kl_after_fw,
      kl_after_em = kl,
      kl_ulb = kl_ulb,
      gr = gr,
      best_gr = best_gr,
      elapsed_s = proc.time()[["elapsed"]] - t_start
    )
    support_sizes <- c(support_sizes, state_n_atoms(state))
    gap_ses <- c(gap_ses, fw$se %||% 0)
    if (verbose) {
      message(sprintf(
        "Iter %d [%d atoms]: Gap %e, KL %e (delta %.1e), ULB %e, KL-ULB %e, GR %e, alpha* %.2e, prop* %.2f",
        fw_idx,
        state_n_atoms(state),
        gap,
        kl,
        kl_prev - kl,
        kl_ulb,
        kl - kl_ulb,
        gr,
        eps_star,
        if (is.na(prop_star)) NaN else prop_star
      ))
    }
    record_checkpoint(fw_idx, oracle_theta_cp = fw$best_theta)
    if (gap < gap_tol) {
      converged <- TRUE
      if (verbose) {
        message(sprintf(
          "Converged after %d FW iterations (gap = %e, kl = %e).",
          fw_idx,
          gap,
          kl
        ))
      }
    }
    kl_prev <- kl
  }

  kl_trace <- do.call(rbind, trace_rows)
  history <- do.call(rbind, outer_rows)
  metrics <- data.frame(
    iteration = history$iter,
    wall_clock = history$elapsed_s,
    gap = history$gap,
    gap_se = gap_ses,
    support_size = support_sizes,
    objective = history$kl_after_em
  )

  # --- Prune, then certify the pruned projection ---
  # Pruning drops the near-zero-weight atoms Frank-Wolfe leaves behind; the
  # certificate is computed on the pruned mixture so it always describes the
  # object we return. With the default `prune_threshold = 0` nothing is dropped
  # and the final oracle sweep is reused.
  n_live <- state_n_atoms(state)
  atoms <- state_atoms(state)[seq_len(n_live)]
  weights <- state_weights(state)[seq_len(n_live)]
  faces <- state_face_idx(state)[seq_len(n_live)]

  keep <- weights > prune_threshold
  if (!any(keep)) {
    stop("`prune_threshold` removed every atom; lower it")
  }
  atoms <- atoms[keep]
  faces <- faces[keep]
  weights <- weights[keep] / sum(weights[keep])

  if (all(keep)) {
    cert <- certify(state, problem, oracle_result = fw)
  } else {
    pruned <- mixture_state(problem$engine, length(atoms))
    for (k in seq_along(atoms)) {
      state_add_atom(pruned, atoms[[k]], faces[k])
    }
    state_set_weights(pruned, weights)
    cert <- certify(pruned, problem)
  }

  projection <- mixture_dist(
    components = do.call(cbind, atoms),
    weights = weights
  )
  ev <- e_variable(
    numerator = problem$alternative,
    projection = projection,
    family = problem$family,
    gap = cert$gap_used
  )

  list(
    e_variable = ev,
    atoms = atoms,
    weights = weights,
    atom_face_idx = faces,
    kl_trace = kl_trace,
    history = history,
    gap = cert$gap,
    oracle_theta = oracle_theta,
    kl = cert$kl,
    kl_ulb = kl_ulb,
    converged = converged,
    checkpoints = checkpoints,
    metrics = metrics,
    certificate = cert,
    state = state
  )
}
