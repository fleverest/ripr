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
#' For Monte Carlo engines, certification draws one fresh sample from Q --
#' independent of the sample `P*` was fitted against, with the engine's own
#' `n_draws` (override with `n_draws`) -- and takes the oracle maximum over the
#' null faces. That maximum is deliberately upward-biased for `sup_theta G`: it
#' is the conservative direction, and `gap_used` inflates it further by
#' `qnorm(conf) * gap_se`. `gap_used` is the only quantity that may be used to
#' rescale an e-variable. Certification is therefore non-deterministic for Monte
#' Carlo engines (governed by the global RNG). Exact engines integrate the true Q
#' with `gap_se = 0`; `estimate` and `n_draws` do not apply.
#'
#' `estimate = TRUE` additionally draws a second independent sample and returns
#' `gap_est`, an unbiased point estimate of the gap at the selected `theta*`. It
#' is a diagnostic -- for comparing algorithm variants, monitoring convergence,
#' or deciding whether to keep optimising -- and is biased *downward* relative to
#' `sup_theta G`. **It must never be used to rescale an e-variable.**
#'
#' Pass a different `engine` to certify against a different Q integration -- e.g.
#' an [exact_engine()], or another `mc_engine()`. The default is the problem's
#' own engine (used only as the Q/family specification to resample from).
#'
#' For Monte Carlo engines the gap is inflated via [inflate_gap()] at level
#' `conf`; the returned `gap` is the raw estimate and `gap_used` the inflated
#' value entering both certified quantities.
#'
#' # What `gap_used` does and does not cover
#'
#' `gap_used` accounts for **Monte Carlo error only** -- the sampling error in
#' estimating `G` at a fixed `theta`, reported by `gap_se` and inflated at level
#' `conf`. It does **not** account for optimisation error. [oracle()] maximises
#' a generally non-convex objective by multi-start BFGS, so `E_star` is a lower
#' bound on the true face maximum: if the search misses the global optimum,
#' `gap` is too small, `gap_used` is too small, and the rescaled e-variable is
#' under-corrected. That error is one-sided in the unsafe direction and is
#' currently neither bounded nor estimated.
#'
#' The practical mitigations are to raise `n_seeds` (recorded in the return
#' value for exactly this reason) and to treat a certificate as evidence rather
#' than proof on nulls whose faces have a multimodal `G`. A deterministic upper
#' bound on `sup_theta G` would remove the caveat; none is implemented.
#'
#' There is a second, opposite failure mode, and unlike the first it is
#' detectable. `Var[G_hat(theta)]` grows like `exp(|theta - a|^2)` in the
#' distance from the nearest atom `a`, so a `theta` far from the support needs
#' roughly `exp(|theta - a|^2)` draws before its estimate means anything. A
#' capable oracle will find such a `theta` and report a large `gap` that is
#' pure sampling noise. This direction is *safe* -- the certificate only gets
#' looser -- but it can make it vacuous, and **more draws do not fix it**: the
#' radius the oracle can profitably exploit grows with the sample.
#'
#' The returned `ess` is the Kish effective sample size of the importance ratio
#' at the selected `theta*`: the number of draws the certified gap is
#' effectively built from. Healthy certificates report `ess` within an order of
#' magnitude of `n_draws`; noise-driven ones collapse to single digits, and
#' `certify()` warns below `ess_min`. Treat a warned certificate as
#' uninformative rather than as a large true gap.
#'
#' @param projection The candidate `P*`: a [finite_mixing] over the null-face
#'   atoms, or a [marginal] wrapping one.
#' @param problem Problem list from [ripr_problem()]; supplies the null faces and
#'   the default engine.
#' @param engine Optional `expectation_engine` giving the Q/family to certify
#'   against, overriding `problem$engine`. Must be over the same family and Q.
#' @param estimate Additionally draw a second independent sample and return the
#'   unbiased diagnostic `gap_est` (see Details). Default `FALSE`. Never affects
#'   `gap_used`. Ignored for exact engines.
#' @param n_draws Draws per resample (Monte Carlo engines only). Default `NULL`
#'   uses the engine's own `n_draws`; drawn twice only when `estimate = TRUE`.
#' @param oracle_result Optional precomputed [oracle_step()] result; honoured
#'   only for exact engines (which do not resample).
#' @param n_seeds Oracle seeds per face.
#' @param conf One-sided confidence level for the Monte Carlo gap inflation.
#'   Default 0.95.
#' @param ess_min Warn when the effective sample size behind the certified gap
#'   falls below this (Monte Carlo engines only). Default 100. Set `0` to
#'   silence. See the accuracy section.
#' @return List with `kl`, `gap`, `gap_se`, `gap_used`, `kl_lower_bound`,
#'   `growth_rate`, `oracle_theta`, `oracle_face`, `conf`, `gap_est`, `ess`, and
#'   the settings that produced it: `n_draws` (`NA` for exact engines) and
#'   `n_seeds`. `gap_est` is `NA_real_` unless `estimate = TRUE`; it is an
#'   unbiased diagnostic estimate of the gap at the selected `theta*`, biased
#'   downward relative to `sup_theta G`, and **must never be used to rescale an
#'   e-variable** -- use `gap_used` for that.
#'
#'   `n_seeds` is recorded because it governs an error the certificate does not
#'   otherwise account for: `gap_se` quantifies the Monte Carlo error from
#'   finite `n_draws`, but an under-seeded oracle *undershoots* `sup_theta G`,
#'   biasing `gap_used` downward -- the unsafe direction.
#' @export
certify <- function(
  projection,
  problem,
  engine = NULL,
  estimate = FALSE,
  n_draws = NULL,
  oracle_result = NULL,
  n_seeds = 200L,
  conf = 0.95,
  ess_min = 100
) {
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
  finish <- function(
    kl,
    gap,
    gap_se,
    oracle_theta,
    oracle_face,
    gap_est = NA_real_,
    n_draws_used = NA_integer_,
    ess = NA_real_
  ) {
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
      conf = conf,
      gap_est = gap_est,
      n_draws = n_draws_used,
      n_seeds = as.integer(n_seeds),
      ess = ess
    )
  }

  # Exact engines integrate the true Q: no sampling, `estimate`/`n_draws` moot,
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

  # One fresh sample both selects theta* and supplies the certified gap: the
  # oracle maximum is upward-biased for sup_theta G, which is the safe
  # direction for a one-sided bound.
  e_cert <- resample_engine(engine, n)
  state <- build(e_cert)
  orr <- oracle_step(state, swap(e_cert), n_seeds = n_seeds)
  gap <- orr$E_star - 1
  gap_se <- orr$se %||% 0
  ess <- orr$ess %||% NA_real_
  if (is.finite(ess) && ess < ess_min) {
    warning(sprintf(
      paste0(
        "certified gap rests on an effective sample size of %.0f (of %d ",
        "draws): the oracle selected a theta* that almost no draw from Q ",
        "supports, so `gap` reflects sampling noise more than sup_theta G. ",
        "The certificate stays conservative but is likely far too loose. ",
        "See the accuracy section of ?certify."
      ),
      ess, n
    ), call. = FALSE)
  }

  # Diagnostic only: a second independent sample gives an unbiased estimate of
  # the gap at the theta* selected above. It never feeds `gap_used`.
  gap_est <- NA_real_
  if (estimate) {
    e_est <- resample_engine(engine, n)
    state_est <- build(e_est)
    log_Pw <- state_log_p_mixture(state_est)
    ratio <- exp(nan_to_neginf(
      eval_log_density(e_est, orr$best_theta) - log_Pw
    ))
    gap_est <- expect(e_est, ratio) - 1
  }

  finish(
    state_objective(state)$loss,
    gap,
    gap_se,
    orr$best_theta,
    orr$best_fi,
    gap_est = gap_est,
    n_draws_used = as.integer(n),
    ess = ess
  )
}
