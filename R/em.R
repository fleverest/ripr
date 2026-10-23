# EM sweeps: reweight, then move each atom to its new position.

#' One EM sweep: weights, then atoms
#'
#' Weights first, since the atom M-step conditions on responsibilities.
#' @keywords internal
#' @noRd
em_sweep <- function(state, ld) {
  wt <- exp(state@engine@log_w) * em_responsibilities(state, ld)
  em_atom_step(em_weight_step(state, wt), ld, wt)
}


#' Responsibilities `r_ic = w_c p_c(x_i) / P(x_i)`; rows sum to 1
#' @keywords internal
#' @noRd
em_responsibilities <- function(state, ld) {
  ld_all <- ld(state@mixing@atoms)
  w <- state@mixing@weights
  exp(add_by_col(ld_all - mixture_log_p(ld_all, w), log(w)))
}


#' The M-step for the weights, \eqn{w_c \leftarrow E_Q[r_c]}{w_c <- E_Q[r_c]}
#'
#' Equals `weight_sweep()`, since \eqn{E_Q[r_c] = w_c G(\theta_c)}{E_Q[r_c] = w_c G(theta_c)}.
#' `wt` is the `(M, C)` responsibilities scaled by the quadrature weights.
#' @keywords internal
#' @noRd
em_weight_step <- function(state, wt) {
  new_w <- colSums(wt)
  # Sums to 1 already, since the rows of `wt` sum to the quadrature weights.
  set_mixture(state, weights = new_w / sum(new_w))
}


#' The M-step for the atoms
#'
#' Each atom maximises its responsibility-weighted log-likelihood within its
#' own part, locally from where it is.
#' @keywords internal
#' @noRd
em_atom_step <- function(state, ld, wt) {
  engine <- state@engine
  family <- engine@family
  atoms <- state@mixing@atoms
  if (nrow(atoms) == 0L) {
    return(state)
  }

  moved <- lapply(
    seq_len(nrow(atoms)),
    function(c_i) {
      # Drop nodes with zero responsibility: their log density may be `-Inf`,
      # and `0 * -Inf` is NaN rather than 0. Responsibilities are fixed, so no
      # value changes, and the objective stays finite for atoms on the
      # boundary of their support (where M-step optima genuinely land).
      keep <- wt[, c_i] > 0
      w_c <- wt[keep, c_i]
      nodes_c <- engine@nodes[keep, , drop = FALSE]
      obj <- objective(
        value = function(theta) {
          sum(w_c * as.vector(ld(matrix(theta, nrow = 1L)))[keep])
        },
        grad = function(theta) {
          as.vector(crossprod(score(family, theta, nodes_c), w_c))
        },
        value_batch = function(theta_mat) {
          as.vector(crossprod(ld(theta_mat)[keep, , drop = FALSE], w_c))
        }
      )
      res <- maximise_over(
        parts(state@null@region)[[state@part[c_i]]],
        obj,
        seeds = atoms[c_i, , drop = FALSE],
        n_seeds = 0L,
        n_restarts = 1L
      )
      # A step that found nothing finite learned nothing: keep the atom.
      if (is.finite(res$value)) res$theta else atoms[c_i, ]
    }
  )
  set_mixture(
    state,
    atoms = matrix(
      unlist(moved, use.names = FALSE),
      nrow = nrow(atoms),
      byrow = TRUE
    )
  )
}
