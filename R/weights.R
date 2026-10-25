# Weight updates on a fixed support.

#' Re-solve every weight with the atoms fixed, then drop the empty atoms
#'
#' The fully-corrective update for verbs that re-solve the weights after
#' moving or adding atoms: merge coincident atoms (two copies of one atom make
#' the hessian degenerate), then `solve_weights()` from the current weights,
#' kept only if it does not raise KL. Returns `list(state, log_p)`, `log_p` being
#' the new mixture's at the engine's nodes.
#' @keywords internal
#' @noRd
correct_state <- function(state, ld) {
  state <- merge_state(state)
  engine <- state@engine
  ld_all <- ld(state@mixing@atoms)
  w <- state@mixing@weights
  log_p <- mixture_log_p(ld_all, w)
  if (length(w) > 1L) {
    solved <- solve_weights(ld_all, w, engine, state@control)
    solved_log_p <- mixture_log_p(ld_all, solved)
    if (kl_at(engine, solved_log_p) <= kl_at(engine, log_p)) {
      w <- solved
      log_p <- solved_log_p
    }
  }
  list(state = drop_empty(set_mixture(state, weights = w)), log_p = log_p)
}


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
#' convex problem. `control$weight_solver` picks the method: `"sqp"` runs
#' `sqp_weights()` to a KKT residual of `weight_sqp_tol`, falling back to
#' SLSQP from where it stopped (keeping whichever has the lower KL) if it
#' stops short; `"slsqp"` runs `slsqp_weights()` alone. Each takes at most
#' `weight_max_iter` iterations or evaluations.
#' @keywords internal
#' @noRd
solve_weights <- function(ld_all, w, engine, control) {
  if (identical(control$weight_solver, "slsqp")) {
    return(slsqp_weights(
      ld_all,
      w,
      engine,
      control$weight_slsqp_tol,
      control$weight_max_iter
    ))
  }
  fit <- sqp_weights(
    ld_all,
    w,
    engine,
    tol = control$weight_sqp_tol,
    max_iter = control$weight_max_iter
  )
  if (fit$kkt <= control$weight_sqp_tol) {
    return(fit$weights)
  }
  fallback <- slsqp_weights(
    ld_all,
    fit$weights,
    engine,
    control$weight_slsqp_tol,
    control$weight_max_iter
  )
  kl <- \(x) kl_at(engine, mixture_log_p(ld_all, x))
  if (kl(fallback) < kl(fit$weights)) fallback else fit$weights
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


#' Weights by exact-Hessian SQP
#'
#' Solves the equivalent problem
#' \deqn{\min_{w \ge 0} -\sum_x q(x) \log (A w)_x + \sum_c w_c,}{min_{w >= 0} -sum_x q(x) log (A w)_x + sum_c w_c,}
#' whose minimiser lies on the simplex: at a KKT point `G_c = 1` on the
#' support, so `sum_c w_c = sum_c w_c G_c = sum_x q(x) = 1`. This is the
#' formulation and the method of mix-SQP \insertCite{Kim2020}{ripr} from the
#' NPMLE literature, without its low-rank approximation of `A`: `sqp_nonneg()`
#' with gradient `1 - G` and Hessian `A^T diag(q / (A w)^2) A`, so its KKT
#' residual is `max(|G_c - 1|` on the support, `G_c - 1` off it`)`.
#'
#' Rows of `A` are scaled by their maximum (a constant shift; `G` is
#' unchanged), outcomes with `q = 0` are dropped, and a ridge of `1e-12` of the
#' largest diagonal keeps the model well-posed when two atoms have nearly the
#' same density. Returns `list(weights, kkt, iterations)`; `kkt` above `tol`
#' means it stopped short (iteration cap, or no descent found).
#' @keywords internal
#' @noRd
sqp_weights <- function(ld_all, w, engine, tol, max_iter) {
  n <- ncol(ld_all)
  q <- exp(engine@log_w)
  keep <- q > 0
  row_max <- matrixStats::rowMaxs(ld_all)
  a <- exp(ld_all - ifelse(is.finite(row_max), row_max, 0))[
    keep,
    ,
    drop = FALSE
  ]
  q <- q[keep]

  w <- pmax(w, 0)
  if (sum(w) <= 0) {
    w <- rep(1 / n, n)
  }
  # Start where every outcome with q > 0 has positive mixture density.
  if (any(a %*% w <= 0)) {
    w <- 0.5 * w / sum(w) + 0.5 / n
  }

  fit <- sqp_nonneg(
    w,
    fn = function(w) {
      y <- as.vector(a %*% w)
      if (any(y <= 0)) Inf else -sum(q * log(y)) + sum(w)
    },
    gr = function(w) 1 - as.vector(crossprod(a, q / as.vector(a %*% w))),
    hess = function(w) {
      h <- crossprod(a * (sqrt(q) / as.vector(a %*% w)))
      diag(h) <- diag(h) + 1e-12 * max(diag(h))
      h
    },
    tol = tol,
    max_iter = max_iter
  )
  list(
    weights = fit$par / sum(fit$par),
    kkt = fit$kkt,
    iterations = fit$iterations
  )
}
