#' @include null.R quadrature.R numerics.R state.R weights.R em.R oracles.R step_rules.R
NULL

# Oracles and step rules
#
# Oracles propose where a new atom should go. Step rules decide whether or not
# to add it and how to update the weights accordingly. EM sweeps also move
# atoms but do not change the support size by adding or removing atoms.
#
# Oracles and steps work on `(ld_all, w)`: the log density of every atom at the
# quadrature nodes, and the weights. Atoms do not move, so the caller builds
# `ld_all`; `log_p` is rebuilt after each weight update. The candidate atom is
# always the last column of `ld_all`, entering with weight zero, which is where
# `add_atom()` puts it.

# --- Step verbs ---------------------------------------------------------------

#' Search every cell; return the best `theta`, `value` and `part`
#'
#' `seeds` defaults to the current atoms: `sum_c w_c G(theta_c) = 1` forces
#' `max_c G(theta_c) >= 1`, so including them keeps the maximum at least one.
#' @keywords internal
#' @noRd
search_null <- function(state, obj, seeds = state@mixing@atoms) {
  ctl <- state@control
  found <- lapply(
    state@null@cells,
    \(s) {
      maximise_over(
        s,
        obj,
        seeds = seeds,
        n_seeds = ctl$n_seeds,
        n_restarts = ctl$n_restarts
      )
    }
  )
  best <- which.max(vapply(found, \(f) f$value, numeric(1)))
  c(found[[best]], list(part = state@null@cell_part[[best]]))
}


#' Write a planned step back to the state
#'
#' The candidate joins only if the step weighted it. Incumbents driven towards
#' zero stay, since only an oracle can grow the support back.
#' @keywords internal
#' @noRd
commit_step <- function(state, theta, part, planned) {
  stepped <- if (planned$uses_candidate) {
    add_atom(state, theta, part, planned$weights)
  } else {
    set_mixture(state, weights = planned$weights[-length(planned$weights)])
  }
  drop_empty(stepped)
}


# --- Support identification ---------------------------------------------------

#' Log density with atom `c_i` removed and the rest renormalised
#' @keywords internal
#' @noRd
log_p_without <- function(log_p, ld_c, w_c) {
  # `w_c p_c <= P` holds exactly but can fail by an ulp when one atom carries
  # nearly all the mass, sending `log1p` to NaN rather than -Inf.
  share <- pmin(w_c * exp(ld_c - log_p), 1)
  log_p + log1p(-share) - log1p(-min(w_c, 1 - .Machine$double.eps))
}


#' Zero every atom whose removal does not increase KL
#'
#' One pass, lightest first, each applied before the next test.
#'
#' **Only safe once the atoms have stopped moving.** Nothing restores a zeroed
#' atom (an oracle will not re-propose a point with `G < 1`), so mid-fit it
#' ratchets KL up.
#' @keywords internal
#' @noRd
identify_support <- function(ld_all, w, engine) {
  log_p <- mixture_log_p(ld_all, w)
  kl <- kl_at(engine, log_p)

  for (c_i in order(w)) {
    active <- which(w > 0)
    if (w[c_i] <= 0 || length(active) < 2L) {
      next
    }
    trial <- log_p_without(log_p, ld_all[, c_i], w[c_i])
    kl_trial <- kl_at(engine, trial)
    if (kl_trial <= kl) {
      w[c_i] <- 0
      w <- w / sum(w)
      log_p <- trial
      kl <- kl_trial
    }
  }
  w
}
