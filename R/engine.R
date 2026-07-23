#' Expectation engine: the seam between algorithms and the sample space
#'
#' An `expectation_engine` owns everything needed to take expectations under the
#' alternative distribution Q: for finite sample spaces this is the full support
#' enumeration and the vector of Q-masses; for continuous spaces a fixed set of
#' Monte Carlo draws from Q (common random numbers). All downstream objective
#' quantities -- the Frank-Wolfe objective `G(theta) = E_Q[p_theta / p_W]`, its
#' gradients, the KL divergence, and the duality gap -- are written against the
#' generics below and never against a concrete engine.
#'
#' Every engine stores a working set of `M` outcomes and a log-weight vector
#' `log_q_mass` such that `E_Q[f] = sum_i exp(log_q_mass_i) f(outcome_i)`. Exact
#' engines enumerate the support with true Q-masses; Monte Carlo engines store
#' draws from Q with uniform weights `1/M`, making the same expressions the plain
#' Monte Carlo estimators. The generics accept base R vectors/matrices defined
#' over the engine's outcome rows, all in log space.
#'
#' The properties are populated by the concrete engine constructors
#' ([exact_engine()], [mc_engine()]); this abstract base is not constructed
#' directly.
#' @param family,alternative The problem components the engine integrates.
#' @param id Character identifier for the engine.
#' @param M,M_f Outcome counts: all outcomes, and the positive-Q-mass subset.
#' @param outcomes The `(M, d)` outcome set the engine integrates over.
#' @param log_q_mass,q_mass,q_mass_f Q-mass over the outcome set: log, linear,
#'   and the positive-mass subset.
#' @param finite_q Logical mask of the positive-Q-mass outcomes.
#' @param H_Q Shannon entropy of Q.
#' @export
expectation_engine <- new_class(
  "expectation_engine",
  abstract = TRUE,
  properties = list(
    family = class_any,
    alternative = class_any,
    id = class_character,
    M = class_numeric,
    outcomes = class_any,
    log_q_mass = class_any,
    q_mass = class_any,
    finite_q = class_any,
    q_mass_f = class_any,
    M_f = class_numeric,
    H_Q = class_numeric
  )
)

#' `E_Q[f(X)]` for a pointwise integrand
#'
#' Contract: `f` is a length-M vector over the engine's outcome set; returns a
#' numeric scalar. `0 * Inf` products are treated as 0. With
#' `support = "q_positive"`, `f` is instead a length-M_f vector over only the
#' outcomes with positive Q-mass (see [q_positive_mask()]).
#' @param engine An `expectation_engine`.
#' @param f Integrand over the outcome set.
#' @param ... Ignored.
#' @param support `"all"` or `"q_positive"`.
#' @return Numeric scalar.
#' @export
expect <- new_generic("expect", "engine", function(engine, f, ..., support = "all") {
  S7::S7_dispatch()
})

#' `E_Q[exp(log_num - log_den)]` via a stable log-space sum
#'
#' Contract: `log_num` and `log_den` are length-M log-density vectors. Returns a
#' numeric scalar. NaNs from `-Inf - -Inf` are treated as log(0). This is the
#' Frank-Wolfe objective `G(theta)` when `log_num = log p_theta` and
#' `log_den = log P_W`.
#' @param engine An `expectation_engine`.
#' @param log_num Length-M log numerator.
#' @param log_den Length-M log denominator.
#' @return Numeric scalar.
#' @keywords internal
expect_ratio <- new_generic(
  "expect_ratio",
  "engine",
  function(engine, log_num, log_den) {
    S7::S7_dispatch()
  }
)

#' Batched [expect_ratio()]: one value per column of `log_num_mat`
#'
#' Contract: `log_num_mat` is `(M, C)`; `log_den` is length-M. Returns a numeric
#' vector of length C.
#' @param engine An `expectation_engine`.
#' @param log_num_mat `(M, C)` log numerators.
#' @param log_den Length-M log denominator.
#' @return Numeric vector of length C.
#' @keywords internal
expect_ratio_batch <- new_generic(
  "expect_ratio_batch",
  "engine",
  function(engine, log_num_mat, log_den) {
    S7::S7_dispatch()
  }
)

#' `E_Q[exp(log_num - log_den) * f(X)]` for a vector-valued integrand
#'
#' Contract: `f_mat` is `(M, D)`; returns a numeric vector of length D. NaN
#' products (zero weight times infinite integrand) are zeroed. This is the
#' gradient of `G(theta)` when `f_mat` is the score of `p_theta`.
#' @param engine An `expectation_engine`.
#' @param log_num Length-M log numerator.
#' @param log_den Length-M log denominator.
#' @param f_mat `(M, D)` integrand matrix.
#' @return Numeric vector of length D.
#' @keywords internal
expect_ratio_grad <- new_generic(
  "expect_ratio_grad",
  "engine",
  function(engine, log_num, log_den, f_mat) {
    S7::S7_dispatch()
  }
)

#' `E_Q[exp(log_w(X)) * f(X)]` for a log-space reweighting
#'
#' Contract: `log_w` is a length-M vector of log-weights; `f` is length-M
#' (returning a scalar) or `(M, D)` (returning a length-D vector). Weights with
#' zero Q-mass contribute nothing. This is the EM M-step objective
#' (`f = log p_theta`) and its gradient (`f = em_score`).
#' @param engine An `expectation_engine`.
#' @param log_w Length-M log-weights.
#' @param f Length-M vector or `(M, D)` matrix.
#' @return Numeric scalar or length-D vector.
#' @keywords internal
expect_weighted <- new_generic(
  "expect_weighted",
  "engine",
  function(engine, log_w, f) {
    S7::S7_dispatch()
  }
)

#' Standard error of the [expect()] estimate
#'
#' Contract: exactly 0 for exact (enumerating) engines; Monte Carlo engines
#' return the estimated standard error of the mean of `f`.
#' @param engine An `expectation_engine`.
#' @param f Integrand over the outcome set.
#' @param ... Ignored.
#' @return Numeric scalar.
#' @keywords internal
expect_se <- new_generic("expect_se", "engine", function(engine, f, ...) {
  S7::S7_dispatch()
})

#' Is every expectation this engine produces deterministic and error-free?
#'
#' Contract: `TRUE` for enumerating engines, `FALSE` for Monte Carlo engines.
#' @param engine An `expectation_engine`.
#' @return Logical scalar.
#' @keywords internal
deterministic <- new_generic("deterministic", "engine", function(engine) S7::S7_dispatch())

#' Log density of a candidate parameter over the engine's outcome set
#'
#' Contract: a length-M vector of `log p_theta` evaluated exactly where the
#' engine integrates. Algorithms must use this (not the family directly) so
#' integrands and weights always align.
#' @param engine An `expectation_engine`.
#' @param theta Parameter vector.
#' @return Length-M numeric vector.
#' @keywords internal
eval_log_density <- new_generic("eval_log_density", "engine", function(engine, theta) S7::S7_dispatch())

#' Batched [eval_log_density()]: `(M, C)` matrix for parameter columns
#' @param engine An `expectation_engine`.
#' @param theta_mat `(d, C)` parameter columns.
#' @return `(M, C)` matrix.
#' @keywords internal
eval_log_density_batch <- new_generic("eval_log_density_batch", "engine", function(engine, theta_mat) S7::S7_dispatch())

#' Score of a candidate parameter over the engine's outcome set
#' @param engine An `expectation_engine`.
#' @param theta Parameter vector.
#' @return `(M, d)` matrix.
#' @keywords internal
eval_score <- new_generic("eval_score", "engine", function(engine, theta) S7::S7_dispatch())

#' EM M-step score over the engine's outcome set (see [em_score()])
#' @param engine An `expectation_engine`.
#' @param theta Parameter vector.
#' @return `(M, d)` matrix.
#' @keywords internal
eval_em_score <- new_generic("eval_em_score", "engine", function(engine, theta) S7::S7_dispatch())

#' Stable identifier for the engine
#' @param engine An `expectation_engine`.
#' @return Character scalar.
#' @keywords internal
engine_id <- new_generic("engine_id", "engine", function(engine) S7::S7_dispatch())

#' Number of outcomes in the engine's working set
#' @param engine An `expectation_engine`.
#' @return Integer.
#' @keywords internal
n_outcomes <- new_generic("n_outcomes", "engine", function(engine) S7::S7_dispatch())

#' Mask of outcomes with positive Q-mass
#' @param engine An `expectation_engine`.
#' @return Logical vector of length M.
#' @keywords internal
q_positive_mask <- new_generic("q_positive_mask", "engine", function(engine) S7::S7_dispatch())

#' Shannon entropy of Q over the outcome set
#'
#' Contract: `H(Q) = -E_Q[log q(X)]`, computed once at construction.
#' `KL(Q || P_W) = -H(Q) - E_Q[log P_W]`.
#' @param engine An `expectation_engine`.
#' @return Numeric scalar.
#' @keywords internal
entropy_q <- new_generic("entropy_q", "engine", function(engine) S7::S7_dispatch())

# --- Shared estimator arithmetic (dispatches for every concrete engine) ---

method(expect, expectation_engine) <- function(engine, f, ..., support = "all") {
  if (identical(support, "q_positive")) {
    sum(nan_to_zero(engine@q_mass_f * f))
  } else {
    sum(nan_to_zero(engine@q_mass * f))
  }
}

method(expect_ratio, expectation_engine) <- function(engine, log_num, log_den) {
  log_terms <- nan_to_neginf(log_num + engine@log_q_mass - log_den)
  exp(logsumexp_vec(log_terms))
}

method(expect_ratio_batch, expectation_engine) <- function(
  engine,
  log_num_mat,
  log_den
) {
  log_terms <- sweep(
    as.matrix(log_num_mat),
    1L,
    engine@log_q_mass - log_den,
    "+"
  )
  exp(col_logsumexp(nan_to_neginf(log_terms)))
}

method(expect_ratio_grad, expectation_engine) <- function(
  engine,
  log_num,
  log_den,
  f_mat
) {
  weights_x <- exp(nan_to_neginf(log_num + engine@log_q_mass - log_den))
  colSums(nan_to_zero(as.matrix(f_mat) * weights_x))
}

method(expect_weighted, expectation_engine) <- function(engine, log_w, f) {
  weights_x <- nan_to_zero(exp(engine@log_q_mass + log_w))
  if (is.null(dim(f))) {
    sum(nan_to_zero(weights_x * f))
  } else {
    colSums(nan_to_zero(as.matrix(f) * weights_x))
  }
}

method(expect_se, expectation_engine) <- function(engine, f, ...) {
  0
}

method(deterministic, expectation_engine) <- function(engine) {
  TRUE
}

method(eval_log_density, expectation_engine) <- function(engine, theta) {
  log_density(engine@family, theta)
}

method(eval_log_density_batch, expectation_engine) <- function(
  engine,
  theta_mat
) {
  log_density_batch(engine@family, theta_mat)
}

method(eval_score, expectation_engine) <- function(engine, theta) {
  score(engine@family, theta)
}

method(eval_em_score, expectation_engine) <- function(engine, theta) {
  em_score(engine@family, theta)
}

method(engine_id, expectation_engine) <- function(engine) {
  engine@id
}

method(n_outcomes, expectation_engine) <- function(engine) {
  engine@M
}

method(q_positive_mask, expectation_engine) <- function(engine) {
  engine@finite_q
}

method(entropy_q, expectation_engine) <- function(engine) {
  engine@H_Q
}

#' Exact expectation engine over an enumerated finite support
#'
#' Enumerates the family's full support once, evaluates Q's mass on it once, and
#' serves every expectation as a dense weighted sum in log space.
#'
#' @param family A `family` with a finite, enumerable support.
#' @param alternative An `alternative`; its induced outcome distribution is the
#'   Q every expectation integrates against.
#' @return An `exact_engine`.
#' @export
exact_engine <- new_class(
  "exact_engine",
  parent = expectation_engine,
  constructor = function(family, alternative) {
    family <- as_family(family)
    alternative <- as_alternative(alternative)
    support_x <- support(family)
    M <- nrow(support_x)
    log_q_mass <- q_log_density(alternative, family, support_x)
    q_mass <- exp(log_q_mass)
    finite_q <- q_mass > 0
    q_mass_f <- q_mass[finite_q]
    M_f <- length(q_mass_f)
    H_Q <- -sum(nan_to_zero(q_mass_f * log_q_mass[finite_q]))
    new_object(
      S7_object(),
      family = family,
      alternative = alternative,
      id = paste0("exact-M", M),
      M = M,
      outcomes = support_x,
      log_q_mass = log_q_mass,
      q_mass = q_mass,
      finite_q = finite_q,
      q_mass_f = q_mass_f,
      M_f = M_f,
      H_Q = H_Q
    )
  }
)

#' Monte Carlo expectation engine with common random numbers
#'
#' For families with continuous (or impractically large) sample spaces. Draws
#' `n_draws` outcomes from Q once at construction and reuses the same draw set
#' for every candidate theta, line search, and iteration. Expectations are plain
#' Monte Carlo means; [expect_se()] returns the usual `sd / sqrt(n_draws)`
#' standard error and [deterministic()] is `FALSE`, which switches the
#' certificate layer to the inflated-gap rule.
#'
#' @param family A `family` whose alternative can be sampled from.
#' @param alternative The Q to sample; must implement [q_sample()].
#' @param n_draws Number of common-random-number draws.
#' @param seed Integer seed fully determining the draw set.
#' @return An `mc_engine`.
#' @export
mc_engine <- new_class(
  "mc_engine",
  parent = expectation_engine,
  properties = list(
    n_draws = class_numeric,
    seed = class_numeric
  ),
  constructor = function(family, alternative, n_draws, seed) {
    family <- as_family(family)
    alternative <- as_alternative(alternative)
    n_draws <- as.integer(n_draws)
    outcomes <- q_sample(alternative, family, n_draws, seed = seed)
    log_q_mass <- rep(-log(n_draws), n_draws)
    q_mass <- exp(log_q_mass)
    finite_q <- q_mass > 0
    log_q_at_draws <- q_log_density(alternative, family, outcomes)
    H_Q <- -mean(nan_to_zero(log_q_at_draws))
    new_object(
      S7_object(),
      family = family,
      alternative = alternative,
      id = paste0("mc-N", n_draws, "-seed", seed),
      M = n_draws,
      outcomes = outcomes,
      log_q_mass = log_q_mass,
      q_mass = q_mass,
      finite_q = finite_q,
      q_mass_f = q_mass,
      M_f = n_draws,
      H_Q = H_Q,
      n_draws = as.numeric(n_draws),
      seed = as.numeric(seed)
    )
  }
)

method(expect_se, mc_engine) <- function(engine, f, ...) {
  sd(nan_to_zero(f)) / sqrt(engine@n_draws)
}

method(deterministic, mc_engine) <- function(engine) {
  FALSE
}

method(eval_log_density, mc_engine) <- function(engine, theta) {
  log_density(engine@family, theta, x = engine@outcomes)
}

method(eval_log_density_batch, mc_engine) <- function(engine, theta_mat) {
  log_density_batch(engine@family, theta_mat, x = engine@outcomes)
}

method(eval_score, mc_engine) <- function(engine, theta) {
  score(engine@family, theta, x = engine@outcomes)
}

method(eval_em_score, mc_engine) <- function(engine, theta) {
  # No multinomial-style correction applies off the simplex: the plain score is
  # the M-step gradient for location families.
  score(engine@family, theta, x = engine@outcomes)
}
