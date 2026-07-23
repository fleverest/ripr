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

#' Certify a candidate mixture: duality gap and e-variable growth rate
#'
#' The common evaluation tool for every algorithm and engine. Runs the
#' Frank-Wolfe oracle over all null faces (or consumes a caller-supplied oracle
#' result to avoid a duplicate sweep), and converts the gap into the certified
#' quantities:
#'
#' - `kl_lower_bound = KL - gap`: the Frank-Wolfe lower bound on the optimal KL
#'   divergence.
#' - `growth_rate = KL - log1p(gap)`: the guaranteed e-power of the rescaled
#'   e-variable `(Q / P_W) / (1 + gap)`, whose null expectation is at most 1 by
#'   construction.
#'
#' For Monte Carlo engines the gap is first inflated via [inflate_gap()] at
#' level `conf`; the returned `gap` is the raw estimate and `gap_used` the
#' inflated value entering both certified quantities. Exact engines have
#' `gap_se = 0`, making the two identical.
#'
#' Certification refuses to proceed if any face reports an inexact oracle: an
#' undershooting oracle would silently overstate the certificate.
#'
#' @param state A `mixture_state` for the candidate W.
#' @param problem Problem list from [ripr_problem()].
#' @param oracle_result Optional result of a fresh [oracle_step()] at this
#'   state; when NULL the oracle is run here with `n_seeds` seeds per face.
#' @param n_seeds Oracle seeds per face when `oracle_result` is NULL.
#' @param conf One-sided confidence level for the Monte Carlo gap inflation.
#'   Default 0.95.
#' @return List with `kl`, `gap`, `gap_se`, `gap_used`, `kl_lower_bound`,
#'   `growth_rate`, `oracle_theta`, `oracle_face`, and `conf`.
#' @export
certify <- function(
  state,
  problem,
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
  if (is.null(oracle_result)) {
    oracle_result <- oracle_step(state, problem, n_seeds = n_seeds)
  }

  kl <- state_objective(state)$loss
  gap <- oracle_result$E_star - 1
  gap_se <- oracle_result$se %||% 0
  gap_used <- inflate_gap(gap, gap_se, conf = conf)

  list(
    kl = kl,
    gap = gap,
    gap_se = gap_se,
    gap_used = gap_used,
    kl_lower_bound = kl - gap_used,
    growth_rate = kl - log1p(gap_used),
    oracle_theta = oracle_result$best_theta,
    oracle_face = oracle_result$best_fi,
    conf = conf
  )
}
