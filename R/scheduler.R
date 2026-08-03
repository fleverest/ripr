#' @include subnull.R quadrature.R mixture.R
NULL

#' Tuning for the RIPr optimiser
#'
#' `ripr_control` changes how well the answer is computed, never what is being
#' computed. `alternative`, `null`, and the engine determine the answer;
#' anything in this object may be varied without changing what the fit converges
#' to. Other operations such as pruning alter the returned mixture, which belong
#' to [ripr_finish()].
#'
#' @param fw_step How the new atom is mixed in. `"line-search"` minimises KL
#'   exactly over the mixing weight; `"fixed"` uses the classical `2/(k+2)`
#'   schedule; `"fully-corrective"` re-optimises every weight over the active
#'   set to convergence after adding the atom; `"none"` adds no atom at all,
#'   leaving EM to refine the mixture it was given. `"pairwise"` and `"away"`
#'   are recognised but not yet implemented.
#' @param em_max_iter Maximum EM sweeps after each Frank-Wolfe step. Each sweep
#'   updates the weights and then relocates the atoms within their own subnulls.
#'   The loop stops early once KL improves by less than `kl_tol`, so this is a
#'   cap rather than a count. Zero runs no EM at all, giving plain Frank-Wolfe
#'   with no corrective step; for a weights-only correction use
#'   `fw_step = "fully-corrective"` with `em_max_iter = 0`.
#' @param n_seeds Random chart seeds per oracle call.
#' @param n_restarts Number of top seeds refined by BFGS per oracle call.
#' @param gap_tol Stop when the estimated gap falls below this.
#' @param kl_tol Stop when KL improves by less than this over an iteration.
#' @param fc_tol Convergence tolerance for the corrective weight solve under
#'   `fw_step = "fully-corrective"`, measured as the KKT residual
#'   `max_c G(theta_c) - 1`. Same units as `gap_tol`.
#' @param fc_max_iter Cap on corrective weight sweeps, so a stalled solve
#'   terminates. Reaching it means the solve did not converge, which shows in
#'   the trace as a full run of `"fc"` rows.
#' @param snapshot Record full mixtures in the trace: `"outer"` once per
#'   iteration, `"all"` after every EM sweep too, `"none"` for scalars only.
#' @param verbose Print progress.
#' @return A list of control settings.
#' @export
ripr_control <- function(
  fw_step = c(
    "line-search",
    "fixed",
    "pairwise",
    "away",
    "fully-corrective",
    "none"
  ),
  em_max_iter = 10L,
  n_seeds = 200L,
  n_restarts = 25L,
  gap_tol = 1e-8,
  kl_tol = 1e-12,
  fc_tol = 1e-10,
  fc_max_iter = 500L,
  snapshot = c("none", "outer", "all"),
  verbose = FALSE
) {
  fw_step <- match.arg(fw_step)
  snapshot <- match.arg(snapshot)
  if (fw_step == "none" && em_max_iter == 0L) {
    warning(
      "with `fw_step = 'none'` and `em_max_iter = 0` nothing can change: no atom ",
      "is added, no weight is updated, and no atom moves.",
      call. = FALSE
    )
  }
  list(
    fw_step = fw_step,
    em_max_iter = as.integer(em_max_iter),
    n_seeds = as.integer(n_seeds),
    n_restarts = as.integer(n_restarts),
    gap_tol = gap_tol,
    kl_tol = kl_tol,
    fc_tol = fc_tol,
    fc_max_iter = as.integer(fc_max_iter),
    snapshot = snapshot,
    verbose = isTRUE(verbose)
  )
}


#' State of a RIPr fit
#'
#' Atoms and weights are lists indexed by subnull: element `i` holds the atoms
#' belonging to \eqn{\Theta_{0i}}{Theta_0i}. The subnull tag matters because the
#' pieces overlap, so an atom in an intersection would otherwise be ambiguous
#' about which chart to search in, and a parallel index vector would make that
#' an invariant to check rather than one to hold structurally. Weights are
#' normalised across the whole list, not within a subnull.
#'
#' @param atoms List of `(d, n_i)` matrices, one per subnull.
#' @param weights List of numeric vectors matching `atoms`, summing to 1 overall.
#' @param alternative The alternative \eqn{Q}{Q}.
#' @param null A [null_model].
#' @param engine A resolved [quadrature].
#' @param control From [ripr_control()].
#' @param trace Data frame, one row per recorded event.
#' @param snapshots List of recorded mixtures.
#' @param iter Completed outer iterations.
#' @param converged Did the fit stop because the estimated gap met `gap_tol`?
#' @return A `ripr_state`.
#' @export
ripr_state <- new_class(
  "ripr_state",
  properties = list(
    atoms = class_list,
    weights = class_list,
    alternative = outcome_distribution,
    null = null_model,
    engine = quadrature,
    control = class_list,
    trace = class_any,
    snapshots = class_list,
    iter = class_numeric,
    converged = class_logical
  ),
  validator = function(self) {
    if (length(self@atoms) != length(self@weights)) {
      return("`weights` must have one element per element of `atoms`")
    }
    if (length(self@atoms) != length(self@null@subnulls)) {
      return("`atoms` must have one element per subnull")
    }
    sizes <- vapply(self@atoms, ncol, integer(1))
    if (!identical(sizes, vapply(self@weights, length, integer(1)))) {
      return("each weight vector must match the column count of its atoms")
    }
    total <- sum(unlist(self@weights))
    if (sum(sizes) > 0L && abs(total - 1) > 1e-8) {
      return("weights must sum to 1 across all subnulls")
    }
    NULL
  }
)


# --- Flat views ---------------------------------------------------------------
# The list-of-matrices shape is the honest one and the one a user reads, but the
# arithmetic wants a single (d, C) matrix. These convert at the boundary.

#' All atoms as one `(d, C)` matrix
#'
#' Columns are grouped by subnull, in subnull order, and a new atom is appended
#' within its own subnull's block. So the flat ordering is not chronological:
#' an atom added at iteration 5 may sit before one added at iteration 2. Use
#' [flat_subnull()] to recover which block a column came from.
#' @param state A [ripr_state].
#' @return `(d, C)` numeric matrix.
#' @export
flat_atoms <- function(state) {
  keep <- vapply(state@atoms, ncol, integer(1)) > 0L
  if (!any(keep)) {
    return(matrix(numeric(0), nrow = 0L, ncol = 0L))
  }
  do.call(cbind, state@atoms[keep])
}


#' All weights as one vector, aligned with [flat_atoms()]
#' @param state A [ripr_state].
#' @return Numeric vector summing to 1.
#' @export
flat_weights <- function(state) unlist(state@weights, use.names = FALSE)


#' Which subnull each column of [flat_atoms()] belongs to
#' @param state A [ripr_state].
#' @return Integer vector.
#' @export
flat_subnull <- function(state) {
  rep(seq_along(state@atoms), vapply(state@atoms, ncol, integer(1)))
}


#' Redistribute a flat weight vector back into the per-subnull list
#' @keywords internal
#' @noRd
unflatten_weights <- function(state, w) {
  sizes <- vapply(state@atoms, ncol, integer(1))
  split(unname(w), rep(seq_along(sizes), sizes))[as.character(seq_along(
    sizes
  ))] |>
    lapply(\(x) if (is.null(x)) numeric(0) else x) |>
    stats::setNames(NULL)
}


#' Redistribute a flat `(d, C)` atom matrix back into the per-subnull list
#' @keywords internal
#' @noRd
unflatten_atoms <- function(state, mat) {
  sizes <- vapply(state@atoms, ncol, integer(1))
  ends <- cumsum(sizes)
  starts <- ends - sizes + 1L
  lapply(seq_along(sizes), function(i) {
    if (sizes[i] == 0L) {
      matrix(numeric(0), nrow = nrow(mat), ncol = 0L)
    } else {
      mat[, starts[i]:ends[i], drop = FALSE]
    }
  })
}


# --- Core quantities ----------------------------------------------------------

#' Log mixture density at the engine's nodes
#' @keywords internal
#' @noRd
log_p_at_nodes <- function(state, ld = NULL) {
  if (is.null(ld)) {
    ld <- compile_engine(state@engine)
  }
  row_logsumexp(add_by_col(ld(flat_atoms(state)), log(flat_weights(state))))
}


#' `KL(Q || P_W)` under the engine's quadrature rule
#' @keywords internal
#' @noRd
kl_divergence <- function(state, log_p = NULL, ld = NULL) {
  if (is.null(log_p)) {
    log_p <- log_p_at_nodes(state, ld)
  }
  expect_q(state@engine, state@engine@log_q - log_p)
}


#' The linear oracle objective `G(theta) = E_theta[p_Q / P_W]`
#'
#' Also its gradient, which is
#' \eqn{E_Q[(p_\theta / P_W) s_\theta]}{E_Q[(p_theta / P_W) score]} by
#' differentiating under the integral.
#' @keywords internal
#' @noRd
fw_objective <- function(state, log_p, ld) {
  engine <- state@engine
  family <- engine@family
  w <- exp(engine@log_w)

  objective(
    value = function(theta) {
      exp(log_expect_q(engine, as.vector(ld(matrix(theta, ncol = 1L))) - log_p))
    },
    grad = function(theta) {
      ratio <- exp(as.vector(ld(matrix(theta, ncol = 1L))) - log_p)
      as.vector(crossprod(score(family, theta, engine@nodes), w * ratio))
    },
    value_batch = function(theta_mat) {
      exp(col_logsumexp(ld(theta_mat) - log_p + engine@log_w))
    }
  )
}


# --- Initialisation -----------------------------------------------------------

#' @keywords internal
#' @noRd
empty_trace <- function() {
  data.frame(
    iter = integer(0),
    phase = character(0),
    inner = integer(0),
    kl = numeric(0),
    gap = numeric(0),
    oracle_value = numeric(0),
    subnull = integer(0),
    step_size = numeric(0),
    support_size = integer(0),
    max_weight = numeric(0),
    stringsAsFactors = FALSE
  )
}


#' Begin a RIPr fit
#'
#' @param alternative The alternative \eqn{Q}{Q}, an [outcome_distribution].
#' @param null A [null_model].
#' @param engine An engine spec, e.g. [exact_engine()].
#' @param atoms Optional list of `(d, n_i)` matrices, one per subnull. `NULL`
#'   places one atom per subnull by projecting the alternative's reference point,
#'   which is the sensible default and what the examples use. Empty subnulls are
#'   `ncol = 0` matrices, which the loop handles without a special case.
#' @param weights Optional list matching `atoms`; defaults to uniform.
#' @param control From [ripr_control()].
#' @return A [ripr_state] with no iterations run.
#' @export
ripr_init <- function(
  alternative,
  null,
  engine = exact_engine(),
  atoms = NULL,
  weights = NULL,
  control = ripr_control()
) {
  resolved <- resolve_engine(engine, alternative, null@family)

  if (is.null(atoms)) {
    # The alternative's modal parameter when there is one: the projection lands
    # near where W_1 puts its mass, which is where the RIPr will be. Falls back
    # to the family's canonical point, since Q need not be a mixture at all.
    ref <- if (S7_inherits(alternative, mixture)) {
      mode_parameter(alternative@mixing)
    } else {
      reference_parameter(null@family)
    }
    atoms <- lapply(null@subnulls, \(s) matrix(init_point(s, ref), ncol = 1L))
  }
  if (length(atoms) != length(null@subnulls)) {
    stop(
      "`atoms` must be a list with one element per subnull (",
      length(null@subnulls),
      "), not ",
      length(atoms),
      ".",
      call. = FALSE
    )
  }
  atoms <- lapply(atoms, as.matrix)

  sizes <- vapply(atoms, ncol, integer(1))
  if (sum(sizes) == 0L) {
    stop("at least one subnull must carry an atom.", call. = FALSE)
  }
  if (is.null(weights)) {
    weights <- lapply(sizes, \(n) rep(1 / sum(sizes), n))
  }

  state <- ripr_state(
    atoms = atoms,
    weights = weights,
    alternative = alternative,
    null = null,
    engine = resolved,
    control = control,
    trace = empty_trace(),
    snapshots = list(),
    iter = 0,
    converged = FALSE
  )
  record(state, phase = "init", inner = 0L)
}


#' Append one row to the trace, and a snapshot if the control asks for it
#' @keywords internal
#' @noRd
record <- function(
  state,
  phase,
  inner,
  kl = NA_real_,
  gap = NA_real_,
  oracle_value = NA_real_,
  subnull = NA_integer_,
  step_size = NA_real_,
  ld = NULL
) {
  if (is.na(kl)) {
    kl <- kl_divergence(state, ld = ld)
  }
  w <- flat_weights(state)
  row <- data.frame(
    iter = as.integer(state@iter),
    phase = phase,
    inner = as.integer(inner),
    kl = kl,
    gap = gap,
    oracle_value = oracle_value,
    subnull = as.integer(subnull),
    step_size = step_size,
    support_size = length(w),
    max_weight = if (length(w)) max(w) else NA_real_,
    stringsAsFactors = FALSE
  )
  state@trace <- rbind(state@trace, row)

  wanted <- switch(
    state@control$snapshot,
    none = FALSE,
    outer = phase %in% c("init", "outer"),
    all = TRUE
  )
  if (wanted) {
    state@snapshots <- c(
      state@snapshots,
      list(list(
        iter = state@iter,
        phase = phase,
        inner = inner,
        atoms = state@atoms,
        weights = state@weights
      ))
    )
  }
  state
}


# --- One iteration ------------------------------------------------------------

#' Run Frank-Wolfe and EM iterations
#'
#' Each iteration sweeps the oracle over every subnull, adds the best new atom
#' according to `fw_step`, then runs EM sweeps.
#'
#' The oracle is seeded with the current atoms: we optimise the objective from
#' around each existing atom, refined via BFGS. Atoms in a near-optimal mixture
#' should have a near-zero gap, hopefully saving us from hitting local optima
#' far from the global maximum.
#'
#' The log-likelihood is compiled once per iteration and passed down, so no
#' inner call recompiles it.
#'
#' @param state A [ripr_state].
#' @param n Iterations to run.
#' @return The updated [ripr_state].
#' @export
ripr_step <- function(state, n = 1L) {
  for (i in seq_len(n)) {
    state <- ripr_iterate(state)
    if (state@converged) break
  }
  state
}


#' @keywords internal
#' @noRd
ripr_iterate <- function(state) {
  ctl <- state@control
  ld <- compile_engine(state@engine)
  log_p <- log_p_at_nodes(state, ld)
  obj <- fw_objective(state, log_p, ld)

  # 1. Linear oracle over each subnull. Seeded with the current atoms.
  seeds <- flat_atoms(state)
  results <- lapply(state@null@subnulls, function(s) {
    maximise_over(
      s,
      obj,
      seeds = seeds,
      n_seeds = ctl$n_seeds,
      n_restarts = ctl$n_restarts
    )
  })
  values <- vapply(results, \(r) r$value, numeric(1))
  best <- which.max(values)
  gap <- values[best] - 1

  state@iter <- state@iter + 1
  if (ctl$verbose) {
    message(sprintf("iter %d: gap %.3e (subnull %d)", state@iter, gap, best))
  }

  # Stop before adding anything. A converged mixture that still gains an atom
  # inflates the reported support size for no benefit.
  if (gap < ctl$gap_tol) {
    state@converged <- TRUE
    return(record(
      state,
      phase = "outer",
      inner = 0L,
      gap = gap,
      oracle_value = values[best],
      subnull = best,
      ld = ld
    ))
  }

  # 2. Frank-Wolfe step.
  state <- fw_update(state, results[[best]]$theta, best, gap, log_p, ld)

  # 3. EM sweeps.
  before <- kl_divergence(state, ld = ld)
  for (k in seq_len(ctl$em_max_iter)) {
    state <- em_step(state, ld)
    after <- kl_divergence(state, ld = ld)
    state <- record(state, phase = "em", inner = k, kl = after, ld = ld)
    if (before - after < ctl$kl_tol) {
      break
    }
    before <- after
  }

  record(state, phase = "outer", inner = 0L, gap = gap, subnull = best)
}


# --- Frank-Wolfe variants -----------------------------------------------------

#' Mix the new atom into the current iterate
#'
#' Dispatches on `control$fw_step`. Every variant here changes only *how fast*
#' the iterate moves, never what it converges to, so the choice belongs in
#' [ripr_control()].
#' @keywords internal
#' @noRd
fw_update <- function(state, theta_new, subnull_index, gap, log_p, ld) {
  ctl <- state@control

  gamma <- switch(
    ctl$fw_step,
    "line-search" = fw_line_search(state, theta_new, log_p, ld),
    # Jaggi's schedule with k counting from 0. `state@iter` has already been
    # incremented, so subtract 1 (k+1 instead of k+2) to get the step index.
    "fixed" = 2 / (state@iter + 1),
    # A line search is only a warm start here: the corrective solve below
    # re-optimises every weight jointly, so the entry weight need only be
    # positive and sensible.
    "fully-corrective" = fw_line_search(state, theta_new, log_p, ld),
    "none" = 0,
    "pairwise" = stop(
      "`fw_step = 'pairwise'` is not implemented. Pairwise steps move mass ",
      "between an existing atom and the new one (Lacoste-Julien and Jaggi, ",
      "2015); use 'line-search' for now.",
      call. = FALSE
    ),
    "away" = stop(
      "`fw_step = 'away'` is not implemented. Away steps remove mass from the ",
      "worst active atom; use 'line-search' for now.",
      call. = FALSE
    )
  )

  if (gamma <= 0) {
    return(record(
      state,
      phase = "fw",
      inner = 0L,
      gap = gap,
      subnull = subnull_index,
      step_size = 0,
      ld = ld
    ))
  }

  state <- add_atom(state, theta_new, subnull_index, gamma)
  state <- record(
    state,
    phase = "fw",
    inner = 0L,
    gap = gap,
    oracle_value = gap + 1,
    subnull = subnull_index,
    step_size = gamma,
    ld = ld
  )

  if (ctl$fw_step == "fully-corrective") {
    state <- correct_weights(state, ld)
  }
  state
}


#' Append an atom to a subnull and rescale all weights by `1 - gamma`
#' @keywords internal
#' @noRd
add_atom <- function(state, theta, subnull_index, gamma) {
  weights <- lapply(state@weights, \(w) w * (1 - gamma))
  weights[[subnull_index]] <- c(weights[[subnull_index]], gamma)

  atoms <- state@atoms
  atoms[[subnull_index]] <- cbind(
    atoms[[subnull_index]],
    theta,
    deparse.level = 0
  )

  S7::set_props(state, atoms = atoms, weights = weights)
}


#' Re-optimise every weight over the active set, to convergence
#'
#' The corrective half of fully-corrective Frank-Wolfe. With the atoms fixed,
#' minimising KL over the weights is convex, and `w <- w * G` is the MM
#' algorithm for it.
#'
#' Convergence is judged on the KKT residual `max_c G(theta_c) - 1`, not on the
#' change in KL. KL is flat at the optimum, so a KL-based update asks for a
#' residual of order `sqrt(tol)` and needs far more sweeps to get there -- EM
#' converges linearly here, and slowly near the simplex boundary. The residual
#' is also in the same units as `gap_tol`.
#'
#' The atoms do not move, so their log-likelihoods are computed once rather than
#' once per sweep.
#' @keywords internal
#' @noRd
correct_weights <- function(state, ld) {
  ctl <- state@control
  engine <- state@engine
  ld_atoms <- ld(flat_atoms(state))

  w <- flat_weights(state)
  log_p <- row_logsumexp(add_by_col(ld_atoms, log(w)))

  for (k in seq_len(ctl$fc_max_iter)) {
    g <- exp(col_logsumexp(ld_atoms - log_p + engine@log_w))
    if (max(g) - 1 < ctl$fc_tol) {
      break
    }

    # sum_c w_c G_c = 1 identically, so this needs no renormalisation.
    w <- w * g
    log_p <- row_logsumexp(add_by_col(ld_atoms, log(w)))
    state@weights <- unflatten_weights(state, w)
    state <- record(
      state,
      phase = "fc",
      inner = k,
      kl = expect_q(engine, engine@log_q - log_p)
    )
  }
  state
}


#' Exact line search for the Frank-Wolfe mixing weight
#'
#' Minimises `KL(Q || (1 - gamma) P_W + gamma p_new)` over `gamma` in `[0, 1]`.
#' One-dimensional and convex in `gamma`, so `optimize` is enough.
#' @keywords internal
#' @noRd
fw_line_search <- function(state, theta_new, log_p, ld) {
  log_new <- as.vector(ld(matrix(theta_new, ncol = 1L)))
  log_q <- state@engine@log_q
  engine <- state@engine

  kl_at <- function(gamma) {
    if (gamma <= 0) {
      return(expect_q(engine, log_q - log_p))
    }
    if (gamma >= 1) {
      return(expect_q(engine, log_q - log_new))
    }
    mixed <- row_logsumexp(cbind(log_p + log1p(-gamma), log_new + log(gamma)))
    expect_q(engine, log_q - mixed)
  }
  stats::optimize(kl_at, interval = c(0, 1), tol = 1e-12)$minimum
}


# --- EM -----------------------------------------------------------------------

#' One EM sweep: weights, then atoms
#'
#' Weights first, because the atom M-step conditions on the responsibilities and
#' those are sharper once the weights have been updated.
#'
#' There is no option to move only the weights. Iterated, that is exactly the
#' fully-corrective step, which [ripr_control()] already exposes as
#' `fw_step = "fully-corrective"` and solves to a convergence test rather than a
#' fixed sweep count.
#' @keywords internal
#' @noRd
em_step <- function(state, ld) {
  em_atom_step(em_weight_step(state, ld), ld)
}


#' Responsibilities of each component at each node
#'
#' `(M, C)` matrix whose rows sum to 1: the E-step.
#' @keywords internal
#' @noRd
em_responsibilities <- function(state, ld) {
  log_comp <- add_by_col(ld(flat_atoms(state)), log(flat_weights(state)))
  exp(log_comp - row_logsumexp(log_comp))
}


#' One EM sweep on the mixture weights
#'
#' The M-step for mixture weights with components fixed:
#' \eqn{w_c \leftarrow E_Q[w_c p_c / P_W]}{w_c <- E_Q[w_c p_c / P_W]}. The
#' updated weights sum to 1 automatically, since the responsibilities sum to 1
#' at every node, and the step cannot increase KL.
#'
#' Iterated to convergence this is the exact fully-corrective step, being the MM
#' algorithm for the same convex subproblem; truncated at `em_max_iter` it is an
#' approximation to it.
#' @keywords internal
#' @noRd
em_weight_step <- function(state, ld) {
  engine <- state@engine
  w <- flat_weights(state)
  log_comp <- add_by_col(ld(flat_atoms(state)), log(w))
  log_p <- row_logsumexp(log_comp)
  # Responsibility-weighted mass per component, integrated under Q.
  new_w <- new_w <- exp(col_logsumexp(log_comp - log_p + engine@log_w))
  new_w <- new_w / sum(new_w)
  state@weights <- unflatten_weights(state, new_w)
  state
}


#' One EM sweep on the atom positions
#'
#' The M-step for the components with weights fixed: each atom maximises its own
#' responsibility-weighted log-likelihood
#' \eqn{\sum_i \omega_i r_{ic} \log p_\theta(x_i)}{sum_i omega_i r_ic log p_theta(x_i)}
#' over its own subnull. Constrained, so it goes through [maximise_over()] in
#' that subnull's chart rather than a closed form.
#'
#' Deliberately a *local* refinement: the only starting point is the atom's
#' current position, with no random seeding. Exploration is the Frank-Wolfe
#' oracle's job, and it already searches each subnull globally with `n_seeds`
#' starts and `n_restarts` refinements. An atom stuck in a poor region is fixed
#' by another Frank-Wolfe iteration placing a new atom, not by this step
#' teleporting the old one -- which would also make atoms jump between sweeps
#' for no gain.
#'
#' Because the inner search is heuristic, the sweep is only approximately
#' monotone: it cannot do worse than the chart image of the current atom, which
#' differs from the atom itself by the round-trip error.
#' @keywords internal
#' @noRd
em_atom_step <- function(state, ld) {
  engine <- state@engine
  family <- engine@family
  atoms_flat <- flat_atoms(state)
  if (ncol(atoms_flat) == 0L) {
    return(state)
  }

  resp <- em_responsibilities(state, ld)
  omega <- exp(engine@log_w)
  idx <- flat_subnull(state)

  moved <- vapply(
    seq_len(ncol(atoms_flat)),
    function(c_i) {
      wt <- omega * resp[, c_i]
      obj <- objective(
        value = function(theta) {
          sum(wt * as.vector(ld(matrix(theta, ncol = 1L))))
        },
        grad = function(theta) {
          as.vector(crossprod(score(family, theta, engine@nodes), wt))
        },
        value_batch = function(theta_mat) {
          as.vector(crossprod(ld(theta_mat), wt))
        }
      )
      maximise_over(
        state@null@subnulls[[idx[c_i]]],
        obj,
        seeds = atoms_flat[, c_i, drop = FALSE],
        n_seeds = 0L,
        n_restarts = 1L
      )$theta
    },
    numeric(nrow(atoms_flat))
  )

  state@atoms <- unflatten_atoms(state, matrix(moved, nrow = nrow(atoms_flat)))
  state
}


# --- Finishing ----------------------------------------------------------------

#' Complete a fit and return its mixing measure and projection
#'
#' Pruning lives here rather than in [ripr_control()] because it changes the
#' answer: it drops atoms, so it is something the caller asks for explicitly
#' rather than a tuning knob.
#'
#' The result makes no claim about validity. It is
#' \eqn{\widehat{W}_0}{W0_hat} and \eqn{\widehat{P}^*}{P_star_hat} and nothing
#' more; an e-variable requires a certified bound obtained separately.
#'
#' @param state A [ripr_state].
#' @param prune Drop atoms with weight at or below this, then renormalise.
#' @return A list with `W0` (a [finite_mixing]), `P_star` (a [mixture]),
#'   `kl`, `gap_est`, `converged`, `atoms`, `weights`, `trace`, and `snapshots`.
#' @export
ripr_finish <- function(state, prune = 1e-6) {
  w <- flat_weights(state)
  keep <- w > prune
  if (!any(keep)) {
    stop(
      "no atom has weight above `prune` (",
      prune,
      "); the largest is ",
      signif(max(w), 3),
      ".",
      call. = FALSE
    )
  }

  mixing <- finite_mixing(
    components = flat_atoms(state)[, keep, drop = FALSE],
    weights = w[keep] / sum(w[keep])
  )
  gaps <- state@trace$gap[!is.na(state@trace$gap)]

  list(
    W0 = mixing,
    P_star = mixture(mixing, state@engine@family),
    kl = kl_divergence(state),
    gap_est = if (length(gaps)) utils::tail(gaps, 1L) else NA_real_,
    converged = state@converged,
    atoms = state@atoms,
    weights = state@weights,
    subnull = flat_subnull(state)[keep],
    trace = state@trace,
    snapshots = state@snapshots
  )
}
