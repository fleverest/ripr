#' Compute the reverse information projection (RIPr)
#'
#' The sole public entry point. Runs the reference hybrid schedule on a
#' [ripr_problem()]: optional initial EM, then `fw_iters` outer iterations of
#' (weight-update step, EM refinement, oracle sweep), with per-iteration
#' recording. All variant behaviour comes from which weight-update step
#' `fw_variant` selects; the loop itself is variant-agnostic.
#'
#' The result is deliberately partitioned by provenance: `projection` and
#' `e_variable` are the deliverables, `certificate` holds everything measured on
#' the fresh certification sample, and `history` / `checkpoints` hold everything
#' measured on the fit sample. No quantity appears in more than one of those.
#'
#' @param problem Problem list from [ripr_problem()].
#' @param init_atoms `d x l` matrix of initial atom locations on the null
#'   faces. Each column is one atom.
#' @param init_atom_faces Integer vector of face indices, one per column.
#' @param init_weights Optional starting mixture weights, one per column of
#'   `init_atoms`. Default `NULL` starts uniform. Pass
#'   `checkpoints$final$weights` from a previous result to resume that fit
#'   where it left off rather than re-deriving the weights by EM. Must be
#'   finite and non-negative with a positive sum; renormalised to sum to 1.
#' @param fw_iters Number of Frank-Wolfe iterations. Set to `0L` for pure EM.
#' @param em_iters Max EM iterations per outer step. Set to `0L` to skip EM.
#' @param n_seeds Random Dirichlet/Gaussian seeds per face for the oracle. The
#'   final certification sweep uses `10 * n_seeds`, matching the factor applied
#'   to `certify_draws`: a wider seed pool costs almost nothing (BFGS polishing
#'   is capped inside [oracle()]), and an under-seeded oracle undershoots the
#'   duality gap, which is the unsafe direction for a certificate.
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
#' @param certify_draws For Monte Carlo engines only, the number of fresh draws
#'   from Q used to certify the final mixture. Certification always runs on a new
#'   sample independent of the one `P_W` was fitted against, so the certified gap
#'   and its standard error are not optimistically biased by the fit. Default
#'   `NULL` uses `10 *` the engine's fit sample size. Certification is a single
#'   oracle sweep, not a per-iteration cost, so its cost is linear in this
#'   count. More draws tighten the *typical* certificate, but they do not bound
#'   its worst case: the certified gap is a maximum over the null faces of a
#'   Monte Carlo estimate, and the region the oracle can profitably exploit
#'   moves outward as the sample grows. See the accuracy section of [certify()],
#'   and check `certificate$ess`. Ignored (with a warning if set) for exact
#'   engines, which integrate the true Q and have nothing to resample.
#' @param certify_ess_min Warn when the final certificate rests on fewer than
#'   this many effective draws (Monte Carlo engines only; see [certify()]).
#'   Default 100. Set `0` to silence -- but a warned certificate is
#'   uninformative rather than merely loose, so prefer investigating it.
#' @param certify_bnb Pass a [bnb_control()] to certify the final gap with the
#'   deterministic Bernstein branch-and-bound bound instead of the heuristic
#'   oracle, when the problem qualifies (see [certifiable()]). Default `NULL`
#'   keeps the heuristic behaviour. This is what makes the returned
#'   `e_variable`'s `gap_certified` `TRUE`; it only ever *raises* the gap, so
#'   the rescaled e-variable gets more conservative, not less.
#' @param fw_variant Weight-update variant: `"line-search"` (default),
#'   `"pairwise"`, `"vanilla"`, `"fully-corrective"`, or `"li-barron"` (the
#'   exact-objective greedy step).
#' @param checkpoint_iters Integer vector of outer-iteration indices (0 =
#'   post-init) at which to snapshot the mixture. Default `NULL` records nothing.
#' @param verbose Print per-iteration progress. Default TRUE.
#' @return A list of six elements, grouped so that each one has a single
#'   provenance -- results, the certificate, and the fit record are not mixed.
#'
#' - `projection`: the fitted RIPr `P*` as a [marginal] (its mixing measure is a
#'   [finite_mixing] over the atoms), after pruning.
#' - `e_variable`: the [e_variable()] wrapping `projection` with the numerator
#'   `Q` and `certificate$gap_used`. Call [e_value()] on it to score data. Its
#'   `gap_certified` carries `certificate$gap_certified`, so the guarantee level
#'   survives being passed around on its own.
#' - `certificate`: the [certify()] result for `projection`, computed on a fresh
#'   sample independent of the fit (see `certify_draws`). Everything derived
#'   from the certification sample lives here and nowhere else, including the
#'   certified worst-case `oracle_theta` / `oracle_face` and the `n_draws` /
#'   `n_seeds` that produced it.
#' - `history`: a data.frame with one row per outer iteration, holding every
#'   fit-sample quantity: `gap`, `gap_se`, `kl_after_fw` / `kl_after_em`,
#'   `kl_ulb`, `gr` / `best_gr`, `support_size`, `face_idx`, `eps_star` /
#'   `prop_star`, and `elapsed_s`. Two list columns nest the finer records:
#'   `kl_trace` (a data.frame of the init/FW/EM steps within that iteration) and
#'   `oracle_theta` (the fit-sample oracle argmax at that iteration).
#' - `checkpoints`: a named list of mixture snapshots. `$final` is always
#'   present and describes the returned `projection` *after* pruning -- its
#'   `atoms`, `atom_face_idx`, and `weights` are exactly the `init_atoms`,
#'   `init_atom_faces`, and `init_weights` needed to resume this fit. Iteration
#'   snapshots
#'   (`$iter_0`, `$iter_3`, ...) appear only if `checkpoint_iters` asked for
#'   them, and are recorded *before* pruning.
#' - `converged`: whether the fit-sample gap fell below `gap_tol`. A statement
#'   about the fit only; it can be `TRUE` alongside a large certified gap.
#' @export
run_ripr <- function(
  problem,
  init_atoms,
  init_atom_faces,
  fw_iters,
  em_iters,
  init_weights = NULL,
  n_seeds = 200L,
  kl_atol = 1e-12,
  kl_rtol = 1e-8,
  gap_tol = 1e-8,
  ls_tol = 1e-12,
  removal_thresh = 1e-8,
  mirror_iters = 500L,
  prune_threshold = 0,
  certify_draws = NULL,
  certify_ess_min = 100,
  certify_bnb = NULL,
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
  # `state_set_weights()` does no validation, so a bad vector here would surface
  # much later as a silently wrong objective. Check it at the boundary instead.
  if (is.null(init_weights)) {
    init_weights <- rep(1 / n_init, n_init)
  } else {
    init_weights <- as.numeric(init_weights)
    if (length(init_weights) != n_init) {
      stop(sprintf(
        "length(init_weights) (%d) must equal ncol(init_atoms) (%d).",
        length(init_weights),
        n_init
      ))
    }
    if (anyNA(init_weights) || any(!is.finite(init_weights))) {
      stop("`init_weights` must be finite and non-missing.")
    }
    if (any(init_weights < 0)) {
      stop("`init_weights` must be non-negative.")
    }
    total <- sum(init_weights)
    if (total <= 0) {
      stop("`init_weights` must have a positive sum.")
    }
    init_weights <- init_weights / total
  }

  state <- mixture_state(problem$engine, n_init + fw_iters)
  t_start <- proc.time()[["elapsed"]]

  # --- History accumulators ---
  trace_rows <- list()
  outer_rows <- list()
  support_sizes <- integer(0L)
  gap_ses <- numeric(0L)
  oracle_thetas <- list()
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

  # Checkpoints are keyed by name so the always-present post-prune `final` entry
  # cannot collide with an iteration checkpoint recorded at the same `iter`.
  record_checkpoint <- function(iter, oracle_theta_cp = NULL) {
    if (is.null(checkpoint_iters) || !(iter %in% checkpoint_iters)) {
      return(invisible(NULL))
    }
    snap <- state_snapshot(state)
    checkpoints[[sprintf("iter_%d", iter)]] <<- list(
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
  state_set_weights(state, init_weights)
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
  oracle_thetas[[length(oracle_thetas) + 1L]] <- fw$best_theta
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
    oracle_thetas[[length(oracle_thetas) + 1L]] <- fw$best_theta
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

  # One row per outer iteration. The inner (init/FW/EM) KL trace and the oracle
  # argmax are nested as list columns rather than kept as parallel objects, so
  # every fit-sample record is keyed by the iteration it belongs to.
  kl_trace <- do.call(rbind, trace_rows)
  history <- do.call(rbind, outer_rows)
  history$gap_se <- gap_ses
  history$support_size <- support_sizes
  history$oracle_theta <- oracle_thetas
  history$kl_trace <- lapply(history$iter, function(i) {
    rows <- kl_trace[kl_trace$iter == i, setdiff(names(kl_trace), "iter"), drop = FALSE]
    rownames(rows) <- NULL
    rows
  })

  # --- Prune, then certify the pruned projection ---
  # Pruning drops the near-zero-weight atoms Frank-Wolfe leaves behind; the
  # certificate is computed on the pruned mixture so it always describes the
  # object we return.
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

  proj_mixing <- finite_mixing(
    components = do.call(cbind, atoms),
    weights = weights
  )

  # Certification. Exact engines integrate the true Q, so we certify against the
  # same engine (reusing the final oracle sweep when nothing was pruned). Monte
  # Carlo engines resample fresh, independent draws inside `certify` -- the
  # mixture was optimised against the fit draws, so an in-sample gap would be
  # optimistically biased. Certification is one oracle sweep rather than a
  # per-iteration cost, so `certify_draws` defaults well above the fit size.
  if (deterministic(problem$engine) && !is.null(certify_draws)) {
    warning(
      "`certify_draws` is ignored for the exact (deterministic) engine, ",
      "which integrates the true Q and has no sample to redraw."
    )
  }
  n_cert <- if (deterministic(problem$engine)) {
    NULL
  } else {
    certify_draws %||% (10L * as.integer(problem$engine@n_draws))
  }
  # Certification inherits the run's `n_seeds`, scaled up by the same factor as
  # the draw count. Widening the seed pool is nearly free: `oracle()` caps BFGS
  # polishing at its own `n_restarts` (25), so extra seeds only enlarge the
  # cheap batched screen. Under-seeding *undershoots* sup G, biasing `gap_used`
  # down -- the unsafe direction -- so the certification sweep should never
  # search less thoroughly than the fit sweeps that placed the atoms.
  n_seeds_cert <- 10L * as.integer(n_seeds)
  cert <- certify(
    proj_mixing,
    problem,
    n_draws = n_cert,
    n_seeds = n_seeds_cert,
    ess_min = certify_ess_min,
    bnb = certify_bnb,
    oracle_result = if (deterministic(problem$engine) && all(keep)) fw else NULL
  )
  if (verbose && !deterministic(problem$engine)) {
    message(sprintf(
      "Certified on %d fresh draws: gap %e (se %.2e), gap_used %e",
      n_cert,
      cert$gap,
      cert$gap_se,
      cert$gap_used
    ))
  }

  projection <- as_marginal(proj_mixing, problem$family)
  ev <- e_variable(
    numerator = problem$alternative,
    projection = projection,
    gap = cert$gap_used,
    # The guarantee level travels with the object, so a caller holding only the
    # e-variable can still tell a proven rescaling from a heuristic one.
    gap_certified = isTRUE(cert$gap_certified)
  )

  # The `final` checkpoint is always recorded, after pruning, so its atoms,
  # weights, and face indices describe exactly the returned `projection`. That
  # makes it the round-trip record: its three fields are the `init_atoms`,
  # `init_atom_faces`, and `init_weights` needed to resume this fit.
  # `oracle_theta` is
  # NA because no fit-sample sweep was run on the pruned mixture -- the
  # certified worst case lives in `certificate$oracle_theta`, on its own sample.
  checkpoints$final <- list(
    iter = NA_integer_,
    atoms = atoms,
    weights = weights,
    atom_face_idx = faces,
    oracle_theta = NA
  )

  list(
    projection = projection,
    e_variable = ev,
    certificate = cert,
    history = history,
    checkpoints = checkpoints,
    converged = converged
  )
}
