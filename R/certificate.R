#' Inflate a Monte Carlo duality-gap estimate to a one-sided upper bound
#'
#' The single, auditable inflation rule for stochastic engines: the estimated
#' gap `gap_hat` is replaced by `gap_hat + z_conf * se`, the upper end of a
#' one-sided normal confidence interval at level `conf`. Exact engines have
#' `se = 0`, so the rule reduces to the identity and the exact and Monte Carlo
#' certificate paths share one formula.
#'
#' @param gap_hat Numeric. Estimated duality gap `E_star - 1`.
#' @param se Numeric. Standard error of the gap estimate (0 for exact engines).
#' @param conf Confidence level of the one-sided bound. Default 0.95.
#' @return Numeric: the inflated gap.
#' @keywords internal
inflate_gap <- function(gap_hat, se, conf = 0.95) {
  gap_hat + qnorm(conf) * se
}

#' Certify a projection: duality gap and e-variable growth rate
#'
#' The common evaluation tool for every algorithm and engine. Takes a fitted
#' projection `P*` (a [finite_mixing] of atoms on the null faces, or the
#' [marginal] wrapping one -- e.g. the `projection` from a [run_ripr()] result),
#' runs the Frank-Wolfe oracle over all null faces, and converts the gap into the
#' certified quantities:
#'
#' - `kl_lower_bound = KL - gap`: the Frank-Wolfe lower bound on the optimal KL
#'   divergence.
#' - `growth_rate = KL - log1p(gap)`: the guaranteed e-power of the rescaled
#'   e-variable `(Q / P*) / (1 + gap)`, whose null expectation is at most 1 by
#'   construction.
#'
#' Monte Carlo certification does not reuse the passed engine's own draws: it
#' resamples fresh outcomes from the engine's Q, with the engine's own `n_draws`
#' (override with `n_draws`). By default (`split = TRUE`) it draws **two**
#' independent samples -- one to locate the worst-case `theta*` (the oracle
#' argmax over the faces), a second to estimate the gap and its standard error at
#' that fixed `theta*`. Because `theta*` is then estimated on data independent of
#' the data that selected it, the gap estimate is free of the selection
#' ("winner's curse") bias that a single sample incurs. With `split = FALSE` one
#' sample does both, at half the draws but with an optimistically biased gap and
#' an over-tight `gap_se`. Certification is therefore non-deterministic for Monte
#' Carlo engines (governed by the global RNG). Exact engines integrate the true Q
#' with `gap_se = 0`; `split` and `n_draws` do not apply.
#'
#' Pass a different `engine` to certify against a different Q integration -- e.g.
#' an [exact_engine()], or another `mc_engine()`. The default is the problem's
#' own engine (used only as the Q/family specification to resample from).
#'
#' For Monte Carlo engines the gap is inflated via [inflate_gap()] at level
#' `conf`; the returned `gap` is the raw estimate and `gap_used` the inflated
#' value entering both certified quantities. Certification refuses to proceed if
#' any face reports an inexact oracle: an undershooting oracle would silently
#' overstate the certificate.
#'
#' @param projection The candidate `P*`: a [finite_mixing] over the null-face
#'   atoms, or a [marginal] wrapping one.
#' @param problem Problem list from [ripr_problem()]; supplies the null faces and
#'   the default engine.
#' @param engine Optional `expectation_engine` giving the Q/family to certify
#'   against, overriding `problem$engine`. Must be over the same family and Q.
#' @param split Draw two independent samples -- one to select `theta*`, one to
#'   estimate the gap at it -- removing the selection bias. Default `TRUE`. With
#'   `FALSE`, a single sample does both. Ignored for exact engines.
#' @param n_draws Draws per resample (Monte Carlo engines only). Default `NULL`
#'   uses the engine's own `n_draws`; `split = TRUE` draws this many twice.
#' @param oracle_result Optional precomputed [oracle_step()] result; honoured
#'   only for exact engines (which do not resample).
#' @param n_seeds Oracle seeds per face.
#' @param conf One-sided confidence level for the Monte Carlo gap inflation.
#'   Default 0.95.
#' @return List with `kl`, `gap`, `gap_se`, `gap_used`, `kl_lower_bound`,
#'   `growth_rate`, `oracle_theta`, `oracle_face`, and `conf`.
#' @export
certify <- function(
  projection,
  problem,
  engine = NULL,
  split = TRUE,
  n_draws = NULL,
  oracle_result = NULL,
  n_seeds = 200L,
  conf = 0.95
) {
  exactness <- vapply(problem$null, oracle_exactness, character(1L))
  if (any(exactness != "exact")) {
    stop(
      "cannot certify: face(s) ",
      paste(which(exactness != "exact"), collapse = ", "),
      " report an inexact oracle, so the duality gap may be underestimated. ",
      "Inexact-oracle certification is not implemented."
    )
  }

  mixing <- if (S7_inherits(projection, marginal)) {
    projection@mixing
  } else {
    projection
  }
  if (!S7_inherits(mixing, finite_mixing)) {
    stop(
      "`projection` must be a `finite_mixing` (or a `marginal` wrapping one)"
    )
  }

  engine <- engine %||% problem$engine
  atoms <- lapply(
    seq_len(ncol(mixing@components)),
    function(k) mixing@components[, k]
  )

  # Rebuild the projection's per-outcome columns over `eng`'s outcome set (face
  # indices are irrelevant to certification).
  build <- function(eng) {
    build_mixture_state(eng, atoms, rep(1L, length(atoms)), mixing@weights)
  }
  swap <- function(eng) {
    p <- problem
    p$engine <- eng
    p$family <- eng@family
    p
  }
  finish <- function(kl, gap, gap_se, oracle_theta, oracle_face) {
    gap_used <- inflate_gap(gap, gap_se, conf = conf)
    list(
      kl = kl,
      gap = gap,
      gap_se = gap_se,
      gap_used = gap_used,
      kl_lower_bound = kl - gap_used,
      growth_rate = kl - log1p(gap_used),
      oracle_theta = oracle_theta,
      oracle_face = oracle_face,
      conf = conf
    )
  }

  # Exact engines integrate the true Q: no sampling, `split`/`n_draws` moot,
  # `se = 0`. Reuse a supplied oracle_result when given.
  if (deterministic(engine)) {
    state <- build(engine)
    orr <- oracle_result %||% oracle_step(state, swap(engine), n_seeds = n_seeds)
    return(finish(
      state_objective(state)$loss,
      orr$E_star - 1,
      orr$se %||% 0,
      orr$best_theta,
      orr$best_fi
    ))
  }

  n <- n_draws %||% engine@n_draws

  if (!split) {
    # One fresh sample selects theta* and estimates the gap at it.
    e1 <- resample_engine(engine, n)
    state <- build(e1)
    orr <- oracle_step(state, swap(e1), n_seeds = n_seeds)
    return(finish(
      state_objective(state)$loss,
      orr$E_star - 1,
      orr$se %||% 0,
      orr$best_theta,
      orr$best_fi
    ))
  }

  # Two independent samples: select theta* on the first, estimate the gap (and
  # its se) at that fixed theta* on the second.
  e_sel <- resample_engine(engine, n)
  e_est <- resample_engine(engine, n)
  orr <- oracle_step(build(e_sel), swap(e_sel), n_seeds = n_seeds)
  theta_star <- orr$best_theta

  state_est <- build(e_est)
  log_Pw <- state_log_p_mixture(state_est)
  ratio <- exp(nan_to_neginf(eval_log_density(e_est, theta_star) - log_Pw))
  finish(
    state_objective(state_est)$loss,
    expect(e_est, ratio) - 1,
    expect_se(e_est, ratio),
    theta_star,
    orr$best_fi
  )
}
