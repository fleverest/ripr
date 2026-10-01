#' @include null.R quadrature.R numerics.R state.R
NULL

# Oracles and step rules
#
# Oracles propose where a new atom should go. Step rules decide whether or not
# to add it and how to update the weights accordingly. EM sweeps (at the foot)
# also move atoms but do not change the support size by adding or removing
# atoms.
#
# Oracles and steps work on `(ld_all, w)`: the log density of every atom at the
# quadrature nodes, and the weights. Atoms do not move, so the caller builds
# builds `ld_all`; `log_p` is rebuilt after each weight update. The candidate
# atom is always the last column of `ld_all`, entering with weight zero, which
# is where `add_atom()` puts it.

# --- Oracles ------------------------------------------------------------------

#' KL Minimisation Oracles
#'
#' The RIPr fit minimises \eqn{KL(Q \| P)}{KL(Q || P)} over mixtures `P` whose
#' components lie in the null. An *oracle* proposes the next component to
#' bring in; [fw_step()] and [lb_step()] differ only in how far ahead they look.
#'
#' Write \eqn{P_i}{P_i} for the current mixture and
#' \eqn{G(\theta) = E_Q[P_\theta / P_i] = E_\theta[Q / P_i]}{G(theta) = E_Q[P_theta / P_i] = E_theta[Q / P_i]}.
#'
#' The Frank--Wolfe \insertCite{Jaggi2013}{ripr} linear oracle [fw_step()] asks
#' how fast KL falls as an infinitesimal amount of mass moves towards
#' \eqn{P_\theta}{P_theta}:
#'
#' \deqn{\left.\frac{d}{d\epsilon} KL\!\left(Q \,\|\, (1-\epsilon) P_i + \epsilon P_\theta\right)\right|_{\epsilon=0} = 1 - G(\theta).}{d/d(eps) KL(Q || (1 - eps) P_i + eps P_theta) =  1 - G(theta) at eps=0.}
#'
#' Equivalently, it finds where in the null the likelihood ratio
#' \eqn{R_i=Q/P_i}{R_i=Q/P_i} most fails to be an e-variable. The attained
#' \eqn{\sup_{\theta \in \Theta_0} G(\theta) = 1 + \mathrm{gap}}{sup G = 1 + gap}
#' bounds the KL suboptimality by Frank--Wolfe duality. Since the search is
#' non-convex, the reported gap is only a lower bound on the true gap. It is
#' the gap of the mixture the step started *from*, so it fills the previous
#' row's `gap_after`.
#'
#' The Li--Barron \insertCite{LiBarron1999}{ripr} greedy oracle [lb_step()]
#' does not linearise; it picks the atom minimising KL after choosing its
#' weight:
#'
#' \deqn{\theta^* = \arg\min_{\theta \in \Theta_0} KL\!\left(Q \,\|\, (1 - w(\theta)) P_i + w(\theta) P_\theta\right),}{theta* = argmin_theta KL(Q || (1 - w(theta)) P_i + w(theta) P_theta),}
#'
#' where \eqn{w(\theta)}{w(theta)} comes from the weight update. Each
#' candidate may therefore cost a line search, or a full weight solve under
#' correction. It has no duality bound; `record_gap = TRUE` estimates the
#' Frank--Wolfe gap with an extra linear-oracle sweep, without adding its atom.
#'
#' # What a trace row records
#'
#' Every verb writes one row per step to `state@trace`. A row spans a step, so
#' its columns split by which mixture they measure.
#'
#' \describe{
#'   \item{`step`, `phase`}{The row's index ([ripr_init()]'s row is 0) and the
#'   verb that wrote it (`"init"`, `"fw"`, `"lb"`, `"em"` or `"weight"`).
#'   Count steps of a kind with e.g. `cumsum(trace$phase == "fw")`.}
#'   \item{`oracle_value`, `oracle_theta`}{What the step's own search found at
#'   the mixture it stepped *from*. For [fw_step()], `oracle_theta` is the
#'   linear oracle's maximiser and `oracle_value` is `NA` (its value is the
#'   previous row's `gap_after + 1`). For [lb_step()] `oracle_value` is the
#'   Li--Barron objective, not a gap; for [weight_step()] it is the pre-sweep
#'   support gap plus one, with no `oracle_theta`. When `part` is not `NA`, the
#'   atom went to `oracle_theta`.}
#'   \item{`kl`, `gap_after`, `gap_after_theta`, `gap_after_part`}{The mixture
#'   the row *produced*, and the linear oracle's gap, maximiser and part for
#'   it. Filled by `record_gap = TRUE` or the next [fw_step()]; otherwise `NA`.}
#'   \item{`gap_after_elapsed`}{Seconds the search that filled `gap_after`
#'   took.}
#'   \item{`part`, `step_size`, `direction`, `support_size`, `max_weight`}{
#'   Where an oracle step put its atom (`NA` if none), the step length and
#'   direction, and the size and heaviest weight of the produced mixture.}
#'   \item{`elapsed`}{Seconds of the step rule's work, excluding diagnostics
#'   (`record_gap`, `snapshot`, `until`). An [fw_step()] row includes the
#'   linear-oracle search it consumed, even if cached. Per step; use
#'   `cumsum(trace$elapsed)` for the total.}
#' }
#'
#' Both `theta` columns are list columns holding the family's parameter, with
#' `NA` for rows that recorded none.
#'
#' `state@oracle` caches the linear oracle's result for the current mixture
#' only: a step that changes the mixture clears it and `record_gap = TRUE`
#' sets it. [fw_step()] uses the cache when present, so an fw run pays one
#' search per step. The final row's `gap_after` stays `NA` unless
#' `record_gap = TRUE` or `until` stopped the run; `ripr_finish(record_gap =
#' TRUE)` measures the mixture it returns.
#'
#' # Which to use
#'
#' Prefer [fw_step()]: it is cheaper per iteration with similar guarantees
#' (none of which strictly hold, since the oracle is heuristic).
#'
#' Expect the two to disagree sharply on the gap even if the KL is similar.
#' Frank--Wolfe places its atom exactly at the worst-case \eqn{\theta}{theta},
#' absorbing the point that defines the gap; Li--Barron places it wherever
#' post-step KL is smallest, leaving that point untouched.
#'
#' @name oracles
#' @references
#'   \insertAllCited{}
#' @seealso [fw_step()], [lb_step()]
NULL


#' The linear oracle: maximise `G(theta)` over the null (see [oracles])
#'
#' Gradient is \eqn{E_\theta[(Q/P_i)\, s_\theta]}{E_theta[(Q/P_i) s_theta]},
#' differentiating under the integral.
#' @keywords internal
#' @noRd
linear_oracle <- function(state, log_p, ld) {
  engine <- state@engine
  family <- engine@family
  w <- exp(engine@log_w)

  objective(
    value = function(theta) {
      exp(log_expect_q(engine, as.vector(ld(matrix(theta, nrow = 1L))) - log_p))
    },
    grad = function(theta) {
      ratio <- exp(as.vector(ld(matrix(theta, nrow = 1L))) - log_p)
      as.vector(crossprod(score(family, theta, engine@nodes), w * ratio))
    },
    value_batch = function(theta_mat) atom_g(ld(theta_mat), log_p, engine)
  )
}


#' What a step towards `theta` would do, without doing it
#'
#' Returns a memoised function of `theta` giving `weights`, `log_p`, `kl`,
#' `gamma`, `direction`, `uses_candidate` and `ld_new`, so an oracle can score
#' candidates with the machinery that later takes the step.
#'
#' `correct` re-solves every weight after the step, which may zero the
#' candidate or an incumbent.
#' @keywords internal
#' @noRd
plan_step <- function(
  state,
  log_p,
  ld,
  directions = "forward",
  size = "line-search",
  gamma_fixed = NULL,
  correct = FALSE
) {
  engine <- state@engine
  ctl <- state@control
  ld_atoms <- ld(state@mixing@atoms)
  w <- c(state@mixing@weights, 0)

  memoise_last(function(theta) {
    ld_new <- as.vector(ld(matrix(theta, nrow = 1L)))
    ld_all <- cbind(ld_atoms, ld_new, deparse.level = 0)

    res <- apply_step(
      ld_all,
      w,
      log_p,
      engine,
      directions = directions,
      size = size,
      gamma_fixed = gamma_fixed
    )
    if (correct) {
      res$weights <- solve_weights(
        ld_all,
        res$weights,
        engine,
        tol = ctl$fc_tol,
        max_iter = ctl$fc_max_iter
      )
      res$log_p <- mixture_log_p(ld_all, res$weights)
      res$kl <- kl_at(engine, res$log_p)
      res$uses_candidate <- res$weights[length(res$weights)] > 0
    }
    res$ld_new <- ld_new
    res
  })
}


#' Gradient of the Li--Barron objective
#'
#' By the envelope theorem the weight can be held fixed:
#' \deqn{\partial_\theta E_Q[\log P] = E_Q\!\left[\frac{w(\theta) P_\theta}{P}\, s_\theta\right],}{d/dtheta E_Q[log P] = E_Q[(w(theta) P_theta / P) s_theta],}
#' so no differentiation through the line search or weight solve is needed.
#' Zero when the step left the candidate unweighted.
#' @keywords internal
#' @noRd
lb_gradient <- function(state, theta, planned) {
  w_new <- utils::tail(planned$weights, 1L)
  if (w_new <= 0) {
    return(numeric(length(theta)))
  }
  engine <- state@engine
  share <- w_new * exp(planned$ld_new - planned$log_p)
  as.vector(crossprod(
    score(engine@family, theta, engine@nodes),
    exp(engine@log_w) * share
  ))
}

#' The Li--Barron nonlinear oracle (see [oracles]); `value` is `-kl`
#'
#' `value_batch` scores every seed with a full step (`n_seeds` line searches or
#' weight solves).
#' @keywords internal
#' @noRd
nonlinear_oracle <- function(
  state,
  log_p,
  ld,
  size = "line-search",
  gamma_fixed = NULL,
  correct = FALSE
) {
  after <- plan_step(
    state,
    log_p,
    ld,
    size = size,
    gamma_fixed = gamma_fixed,
    correct = correct
  )

  objective(
    value = function(theta) -after(theta)$kl,
    grad = function(theta) lb_gradient(state, theta, after(theta)),
    value_batch = function(theta_mat) {
      -vapply(
        seq_len(nrow(theta_mat)),
        \(i) after(theta_mat[i, ])$kl,
        numeric(1)
      )
    }
  )
}


# --- Step rules ---------------------------------------------------------------

#' `G(theta_c) = E_Q[p_c / P]` for every column of `ld_all`
#'
#' The KL gradient in the weights is `-G`, so every first-order quantity the
#' step layer needs comes from this vector.
#' @keywords internal
#' @noRd
atom_g <- function(ld_all, log_p, engine) {
  exp(col_logsumexp(ld_all - log_p + engine@log_w))
}


#' The active atom with the smallest `G` (the away vertex), or `NULL` if fewer
#' than two are active, since emptying the only one would empty the mixture
#' @keywords internal
#' @noRd
worst_of <- function(g, w) {
  active <- which(w > 0)
  if (length(active) < 2L) {
    return(NULL)
  }
  active[which.min(g[active])]
}


# Paths through weight space. Each `path_*` returns
# `list(direction, value, gamma_max, w_of, log_p_at)` (first-order value, the
# map `gamma -> w`, and log density along it) so one line search serves all.
# The candidate is the last entry of `w`.

#' Forward: \eqn{w \leftarrow (1-\gamma) w + \gamma e_{new}}{w <- (1 - gamma) w + gamma e_new}
#' @keywords internal
#' @noRd
path_forward <- function(ld_all, w, log_p, value) {
  new_idx <- length(w)
  ld_new <- ld_all[, new_idx]
  list(
    direction = "forward",
    value = value,
    gamma_max = 1,
    w_of = function(gamma) {
      out <- (1 - gamma) * w
      out[new_idx] <- gamma
      out
    },
    log_p_at = function(gamma) {
      if (gamma <= 0) {
        return(log_p)
      }
      if (gamma >= 1) {
        return(ld_new)
      }
      row_logsumexp(cbind(log_p + log1p(-gamma), ld_new + log(gamma)))
    }
  )
}


#' Clamp negatives and renormalise onto the simplex
#' @keywords internal
#' @noRd
normalise_weights <- function(w) {
  w <- pmax(w, 0)
  w / sum(w)
}


#' Pairwise: move mass to the candidate from the worst active atom
#'
#' Capped at that atom's weight. With one active atom this is `path_forward`.
#' @references
#'   \insertRef{LacosteJulienJaggi2015}{ripr}
#' @keywords internal
#' @noRd
path_pairwise <- function(ld_all, w, log_p, worst, value) {
  move <- function(gamma) {
    out <- w
    out[worst] <- out[worst] - gamma
    out[length(out)] <- gamma
    normalise_weights(out)
  }
  list(
    direction = "pairwise",
    value = value,
    gamma_max = w[worst],
    w_of = move,
    log_p_at = function(gamma) mixture_log_p(ld_all, move(gamma))
  )
}


#' Away: \eqn{w \leftarrow (1+\gamma) w - \gamma e_v}{w <- (1 + gamma) w - gamma e_v}
#'
#' Leaves the candidate unused, so only offered alongside another direction.
#' @references
#'   \insertRef{LacosteJulienJaggi2015}{ripr}
#' @keywords internal
#' @noRd
path_away <- function(ld_all, w, log_p, worst, value) {
  # Cap is w_v/(1 - w_v): beyond it the worst atom's weight would go negative.
  gamma_max <- w[worst] / (1 - w[worst])
  w_of <- function(gamma) {
    out <- (1 + gamma) * w
    # At the cap the algebra can leave a positive residual; a drop step must
    # leave exactly zero, since `drop_empty()` keys on it.
    out[worst] <- if (gamma >= gamma_max) 0 else out[worst] - gamma
    normalise_weights(out)
  }
  list(
    direction = "away",
    value = value,
    gamma_max = gamma_max,
    w_of = w_of,
    log_p_at = function(gamma) mixture_log_p(ld_all, w_of(gamma))
  )
}


#' The directions a named Frank--Wolfe variant may move in
#' @keywords internal
#' @noRd
variant_directions <- function(variant) {
  switch(
    variant,
    standard = "forward",
    "away-step" = c("forward", "away"),
    pairwise = "pairwise"
  )
}


#' The paths on offer this step, each with its first-order value
#'
#' Since \eqn{\nabla f = -G}{grad f = -G} and
#' \eqn{\sum_c w_c G_c = 1}{sum_c w_c G_c = 1}: forward is \eqn{G_s - 1}{G_s - 1}
#' (the Frank--Wolfe gap), away \eqn{1 - G_v}{1 - G_v}, pairwise
#' \eqn{G_s - G_v}{G_s - G_v}. Only away can be unavailable (below two active
#' atoms), and is then dropped.
#' @keywords internal
#' @noRd
step_paths <- function(directions, ld_all, w, log_p, engine) {
  g <- atom_g(ld_all, log_p, engine)
  g_new <- g[length(g)]
  worst <- worst_of(g, w)
  paths <- lapply(directions, function(d) {
    switch(
      d,
      forward = path_forward(ld_all, w, log_p, g_new - 1),
      # `worst` is NULL below two active atoms; pairwise then uses the single
      # active atom and coincides with forward. See `path_pairwise`.
      pairwise = {
        v <- if (is.null(worst)) which(w > 0)[1L] else worst
        path_pairwise(ld_all, w, log_p, v, g_new - g[v])
      },
      away = if (!is.null(worst)) {
        path_away(ld_all, w, log_p, worst, 1 - g[worst])
      }
    )
  })
  # `away` is only ever offered alongside `forward`, so dropping it when
  # there is no second active atom still leaves a path.
  Filter(Negate(is.null), paths)
}


#' The open-loop step size `2 / (k + 2)`, for `size = "fixed"`
#'
#' `k` starts at the initial support size (warm starts) and counts every
#' oracle step (`fw` and `lb` rows). Jaggi's `2/(k+2)` from 0 and Li--Barron's
#' `2/(k+1)` from 1 are the same sequence.
#' @references
#'   \insertRef{Jaggi2013}{ripr}
#'
#'   \insertRef{LiBarron1999}{ripr}
#' @keywords internal
#' @noRd
schedule_gamma <- function(state) {
  rows <- state@trace_rows
  phase <- vapply(rows, .subset2, character(1), "phase")
  k <- sum(phase %in% c("fw", "lb")) + rows[[1L]]$support_size
  2 / (k + 2)
}


#' Minimise KL along a step path
#'
#' KL is convex along every path, so a path whose first-order `value` is not
#' positive cannot descend and the step is 0. Otherwise `gamma = 0` stays in
#' range and wins ties, so no step can increase KL or move for nothing.
#' @keywords internal
#' @noRd
line_search <- function(path, engine) {
  gamma_max <- path$gamma_max
  if (path$value <= 0 || !is.finite(gamma_max) || gamma_max <= 0) {
    return(0)
  }
  kl <- \(gamma) kl_at(engine, path$log_p_at(gamma))
  found <- stats::optimize(kl, interval = c(0, gamma_max), tol = 1e-12)
  best <- if (kl(gamma_max) <= found$objective) gamma_max else found$minimum
  if (kl(best) < kl(0)) best else 0
}


#' Take one step, without a state or trace
#'
#' `size = "fixed"` caps `gamma_fixed` at the path's own maximum, since
#' pairwise and away cap below 1. Returns `list(weights, log_p, kl, gamma,
#' direction, uses_candidate)`, `weights` of length `C + 1`.
#' @keywords internal
#' @noRd
apply_step <- function(
  ld_all,
  w,
  log_p,
  engine,
  directions = "forward",
  size = "line-search",
  gamma_fixed = NULL
) {
  paths <- step_paths(directions, ld_all, w, log_p, engine)
  path <- paths[[which.max(vapply(paths, \(p) p$value, numeric(1)))]]

  gamma <- if (size == "fixed") {
    min(gamma_fixed, path$gamma_max)
  } else {
    line_search(path, engine)
  }
  stepped <- path$log_p_at(gamma)
  weights <- pmax(path$w_of(gamma), 0)
  list(
    weights = weights,
    log_p = stepped,
    kl = kl_at(engine, stepped),
    gamma = gamma,
    direction = path$direction,
    # False for away, and for any search that put nothing on the candidate.
    uses_candidate = weights[length(weights)] > 0
  )
}


# --- Fully corrective weights -------------------------------------------------

#' One multiplicative sweep \eqn{w_c \leftarrow w_c G(\theta_c)}{w_c <- w_c G(theta_c)}
#'
#' The exact M-step for the weights (an MM algorithm), monotone by construction.
#' `residual`, \eqn{\max_c G_c - 1}{max_c G_c - 1} before the sweep, is the
#' Frank--Wolfe gap over the current support: an upper bound on the KL still
#' available from reweighting.
#' @keywords internal
#' @noRd
weight_sweep <- function(ld_all, w, log_p, engine) {
  g <- atom_g(ld_all, log_p, engine)
  # `sum_c w_c G_c = 1` identically, so the result needs no renormalisation.
  list(weights = w * g, residual = max(g) - 1)
}


#' Re-optimise every weight over the current atoms
#'
#' Minimises `KL(Q || P_w)` over the simplex with the atoms fixed, a smooth
#' convex problem, by SLSQP. Each row of the likelihoods is scaled by its
#' maximum, which shifts the objective by a constant, so nothing overflows.
#' `tol` is SLSQP's relative tolerance on the weights and `max_iter` its cap on
#' evaluations.
#' @keywords internal
#' @noRd
solve_weights <- function(ld_all, w, engine, tol, max_iter) {
  row_max <- matrixStats::rowMaxs(ld_all)
  a <- exp(ld_all - ifelse(is.finite(row_max), row_max, 0))
  q <- exp(engine@log_w)
  x <- nloptr::slsqp(
    w,
    fn = function(x) -sum(q * log(as.vector(a %*% x))),
    gr = function(x) -as.vector(crossprod(a, q / as.vector(a %*% x))),
    lower = numeric(length(w)),
    heq = function(x) sum(x) - 1,
    heqjac = function(x) matrix(1, nrow = 1L, ncol = length(x)),
    control = list(xtol_rel = tol, maxeval = max_iter)
  )$par
  # SLSQP meets its bounds to rounding; snap those to exact zeros so they drop.
  x[x < 1e-12 * max(x)] <- 0
  normalise_weights(x)
}


# --- EM -----------------------------------------------------------------------

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
