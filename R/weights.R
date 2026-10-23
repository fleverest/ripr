# Weight updates on a fixed support.

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
#' convex problem, from `w`, with the solver settings in `control`
#' ([ripr_control()]).
#' @keywords internal
#' @noRd
solve_weights <- function(ld_all, w, engine, control) {
  slsqp_weights(ld_all, w, engine, control$fc_tol, control$fc_max_iter)
}


#' Solve for weights via SLSQP
#'
#' Each row of the likelihoods is scaled by its maximum, which shifts the
#' objective by a constant, so nothing overflows. `tol` is SLSQP's relative
#' tolerance on the weights and `max_iter` its cap on evaluations.
#' @keywords internal
#' @noRd
slsqp_weights <- function(ld_all, w, engine, tol, max_iter) {
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
