# Step rules: how far to move towards a candidate atom, and which weights
# give up mass (forward, pairwise, away), on `(ld_all, w)`.

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
