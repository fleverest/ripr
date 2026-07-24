#' Candidate mixture state W
#'
#' Owns the atoms, weights, face assignments, and the per-atom log-density
#' buffer columns (with a matching q-positive restriction). The interior is a
#' mutable environment: atoms are appended and replaced in place during
#' optimisation.
#'
#' @param engine An `expectation_engine`; expectations and Q-mass quantities are
#'   drawn from it.
#' @param max_atoms Integer. Buffer capacity: the maximum number of atoms this
#'   state can ever hold.
#' @return A `mixture_state`.
#' @export
mixture_state <- new_class(
  "mixture_state",
  properties = list(
    env = class_any
  ),
  constructor = function(engine, max_atoms) {
    M <- n_outcomes(engine)
    M_f <- engine@M_f

    e <- new.env(parent = emptyenv())
    e$engine <- engine
    e$max_atoms <- as.integer(max_atoms)
    e$lc <- matrix(-Inf, nrow = M, ncol = max_atoms)
    e$lc_f <- matrix(-Inf, nrow = M_f, ncol = max_atoms)
    e$n_live <- 0L
    e$atoms <- list()
    e$face_idx <- integer(0L)
    e$weights <- numeric(0L)

    new_object(S7_object(), env = e)
  }
)

#' Build a mixture_state from explicit atoms, faces, and weights on an engine
#'
#' Rehydrates a fitted mixture against a (possibly different) engine, so its
#' per-outcome log-density columns are recomputed over that engine's outcome set.
#' Used to move a fitted `P_W` onto a fresh certification sample.
#' @param engine An `expectation_engine`.
#' @param atoms List of length-d parameter vectors.
#' @param faces Integer face index per atom.
#' @param weights Numeric mixture weights (same length as `atoms`).
#' @return A `mixture_state`.
#' @keywords internal
build_mixture_state <- function(engine, atoms, faces, weights) {
  state <- mixture_state(engine, length(atoms))
  for (k in seq_along(atoms)) {
    state_add_atom(state, atoms[[k]], faces[k])
  }
  state_set_weights(state, weights)
  state
}

# Evaluate and cache atom k's log-density column (full and q-positive).
write_atom_col <- function(e, k, theta) {
  lm <- eval_log_density(e$engine, theta)
  e$lc[, k] <- lm
  e$lc_f[, k] <- lm[q_positive_mask(e$engine)]
}

#' Number of live atoms
#' @param state A `mixture_state`.
#' @return Integer.
#' @keywords internal
state_n_atoms <- new_generic("state_n_atoms", "state", function(state) S7::S7_dispatch())

method(state_n_atoms, mixture_state) <- function(state) {
  state@env$n_live
}

#' Atom list (length-d parameter vectors)
#' @param state A `mixture_state`.
#' @return List of parameter vectors.
#' @keywords internal
state_atoms <- new_generic("state_atoms", "state", function(state) S7::S7_dispatch())

method(state_atoms, mixture_state) <- function(state) {
  state@env$atoms
}

#' Mixture weights
#' @param state A `mixture_state`.
#' @return Numeric weight vector.
#' @keywords internal
state_weights <- new_generic("state_weights", "state", function(state) S7::S7_dispatch())

method(state_weights, mixture_state) <- function(state) {
  state@env$weights
}

#' Face index per atom
#' @param state A `mixture_state`.
#' @return Integer vector.
#' @keywords internal
state_face_idx <- new_generic("state_face_idx", "state", function(state) S7::S7_dispatch())

method(state_face_idx, mixture_state) <- function(state) {
  state@env$face_idx
}

#' Replace the weight vector (atoms unchanged)
#' @param state A `mixture_state`.
#' @param w New weight vector.
#' @return `state`, invisibly.
#' @keywords internal
state_set_weights <- new_generic("state_set_weights", "state", function(state, w) S7::S7_dispatch())

method(state_set_weights, mixture_state) <- function(state, w) {
  state@env$weights <- w
  invisible(state)
}

#' Append an atom (weights are managed separately by the caller's step rule)
#' @param state A `mixture_state`.
#' @param theta Parameter vector of the new atom.
#' @param face_idx Face index of the new atom.
#' @return `state`, invisibly.
#' @keywords internal
state_add_atom <- new_generic("state_add_atom", "state", function(state, theta, face_idx) S7::S7_dispatch())

method(state_add_atom, mixture_state) <- function(state, theta, face_idx) {
  e <- state@env
  if (e$n_live >= e$max_atoms) {
    stop("mixture_state is full: increase max_atoms")
  }
  e$n_live <- e$n_live + 1L
  write_atom_col(e, e$n_live, theta)
  e$atoms[[e$n_live]] <- theta
  e$face_idx[e$n_live] <- as.integer(face_idx)
  invisible(state)
}

#' Overwrite atom k in place
#' @param state A `mixture_state`.
#' @param k Atom index.
#' @param theta New parameter vector.
#' @param face_idx Optional new face index.
#' @return `state`, invisibly.
#' @keywords internal
state_replace_atom <- new_generic("state_replace_atom", "state", function(state, k, theta, face_idx = NULL) S7::S7_dispatch())

method(state_replace_atom, mixture_state) <- function(
  state,
  k,
  theta,
  face_idx = NULL
) {
  e <- state@env
  write_atom_col(e, k, theta)
  e$atoms[[k]] <- theta
  if (!is.null(face_idx)) {
    e$face_idx[k] <- as.integer(face_idx)
  }
  invisible(state)
}

# The live prefix of the full log-density buffer, (M, n_live).
live_lc <- function(e) e$lc[, seq_len(e$n_live), drop = FALSE]

# The live prefix of the q-positive log-density buffer, (M_f, n_live).
live_lc_f <- function(e) e$lc_f[, seq_len(e$n_live), drop = FALSE]

#' Log mixture density log P_W over the engine's outcome set
#'
#' Pass `w` to evaluate at trial weights without mutating the state.
#' @param state A `mixture_state`.
#' @param w Optional trial weights (defaults to the state's weights).
#' @return Length-M numeric vector.
#' @keywords internal
state_log_p_mixture <- new_generic("state_log_p_mixture", "state", function(state, w = NULL) S7::S7_dispatch())

method(state_log_p_mixture, mixture_state) <- function(state, w = NULL) {
  e <- state@env
  if (is.null(w)) {
    w <- e$weights
  }
  w_log <- log(pmax(w, 1e-300))
  row_logsumexp(sweep(live_lc(e), 2L, w_log, "+"))
}

#' KL objective and weight gradient
#'
#' Returns `list(loss = KL(Q || P_W), grad = dKL/dw)`.
#' @param state A `mixture_state`.
#' @param w Optional trial weights (defaults to the state's weights).
#' @return `list(loss = , grad = )`.
#' @export
state_objective <- new_generic("state_objective", "state", function(state, w = NULL) S7::S7_dispatch())

method(state_objective, mixture_state) <- function(state, w = NULL) {
  e <- state@env
  if (is.null(w)) {
    w <- e$weights
  }
  engine <- e$engine
  qpm <- q_positive_mask(engine)
  log_Pw_f <- state_log_p_mixture(state, w)[qpm]

  diff <- sweep(live_lc_f(e), 1L, log_Pw_f, "-")
  w_mat <- nan_to_zero(exp(diff)) * engine@q_mass_f

  list(
    loss = -entropy_q(engine) -
      expect(engine, log_Pw_f, support = "q_positive"),
    grad = -colSums(w_mat)
  )
}

#' Log responsibilities log r_k(x) of each live atom for each outcome
#' @param state A `mixture_state`.
#' @param log_Pw Length-M mixture log-density.
#' @param w Optional trial weights (defaults to the state's weights).
#' @return `(M, n_live)` matrix of log responsibilities.
#' @keywords internal
state_responsibilities <- new_generic("state_responsibilities", "state", function(state, log_Pw, w = NULL) S7::S7_dispatch())

method(state_responsibilities, mixture_state) <- function(
  state,
  log_Pw,
  w = NULL
) {
  e <- state@env
  if (is.null(w)) {
    w <- e$weights
  }
  w_log <- log(pmax(w, 1e-300))
  log_r <- sweep(sweep(live_lc(e), 2L, w_log, "+"), 1L, log_Pw, "-")
  log_r <- nan_to_neginf(log_r)
  log_r - row_logsumexp(log_r)
}

#' `G(theta_k) = E_Q[p_k / P_W]` for every live atom at once
#' @param state A `mixture_state`.
#' @param log_Pw Length-M mixture log-density.
#' @return Numeric vector of length `n_live`.
#' @keywords internal
state_atom_ratios <- new_generic("state_atom_ratios", "state", function(state, log_Pw) S7::S7_dispatch())

method(state_atom_ratios, mixture_state) <- function(state, log_Pw) {
  e <- state@env
  expect_ratio_batch(e$engine, live_lc(e), log_Pw)
}

#' Cached log density column of atom k
#' @param state A `mixture_state`.
#' @param k Atom index.
#' @return Length-M numeric vector.
#' @keywords internal
state_atom_log_density <- new_generic("state_atom_log_density", "state", function(state, k) S7::S7_dispatch())

method(state_atom_log_density, mixture_state) <- function(state, k) {
  state@env$lc[, k]
}

#' Immutable snapshot of the candidate mixture
#'
#' The checkpoint payload: copies of atoms, weights, and face indices for the
#' live prefix.
#' @param state A `mixture_state`.
#' @return `list(atoms = , weights = , atom_face_idx = )`.
#' @export
state_snapshot <- new_generic("state_snapshot", "state", function(state) S7::S7_dispatch())

method(state_snapshot, mixture_state) <- function(state) {
  e <- state@env
  live <- seq_len(e$n_live)
  list(
    atoms = e$atoms[live],
    weights = e$weights[live],
    atom_face_idx = e$face_idx[live]
  )
}
