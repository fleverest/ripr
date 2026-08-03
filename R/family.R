#' Sampling family: the observation model
#'
#' A `sampling_family` owns the model \eqn{p_\theta(x)}{p_theta(x)} and nothing
#' else; no knowledge of null hypotheses, alternatives, or any optimisation.
#' Families supply a compiled log-likelihood, the score, a sampler, and support
#' enumeration (for finite families).
#'
#' @export
sampling_family <- new_class("sampling_family", abstract = TRUE)


#' Dimension of the parameter vector
#' @param family A [sampling_family].
#' @return Integer parameter dimension.
#' @export
param_dim <- new_generic("param_dim", "family", function(family) {
  S7::S7_dispatch()
})


#' Compile the log-likelihood function for a fixed set of outcomes
#'
#' The single density method a family must implement; [log_density_batch()] and
#' [log_density()] are thin wrappers over it. Returns a function of `theta_mat`
#' alone.
#'
#' Compiling lets a family precompute whatever depends on `x` alone, which pays
#' off because the outcomes stay fixed for a whole fit while `theta` moves
#' thousands of times. Closing over `x` rather than taking a cache argument also
#' means precomputed constants cannot be paired with the wrong outcomes.
#'
#' @param family A [sampling_family].
#' @param x `(M, K)` matrix of outcomes, where `M` is the number of outcomes and
#'   `K` the dimension of the sample space.
#' @return A function of `theta_mat`, a `(d, C)` matrix of parameter columns,
#'   returning the `(M, C)` matrix of log densities at `x`.
#' @export
compile_loglik <- new_generic("compile_loglik", "family", function(family, x) {
  S7::S7_dispatch()
})

#' Batched log density over parameter columns
#'
#' Recompiles on every call, so in a loop over many `theta` at fixed `x`,
#' use [compile_loglik()] once and call its result instead.
#'
#' @param family A [sampling_family].
#' @param theta_mat `(d, C)` matrix of parameter columns, where `d` is the
#'   dimension of the parameter and `C` the number of columns. The same `x` is
#'   used for every column.
#' @param x `(M, K)` matrix of outcomes.
#' @return `(M, C)` matrix of log densities.
#' @export
log_density_batch <- function(family, theta_mat, x) {
  compile_loglik(family, x)(theta_mat)
}


#' Log density `log p_theta(x)`
#'
#' @param family A [sampling_family].
#' @param theta Parameter vector of length [param_dim()].
#' @param x `(M, K)` matrix of outcomes, or a length-`K` vector for one outcome.
#' @return Length-`M` numeric vector.
#' @export
log_density <- function(family, theta, x) {
  as.vector(log_density_batch(family, matrix(theta, ncol = 1L), x))
}


#' Score `d log p_theta(x) / d theta`
#'
#' Per-outcome contributions in the family's own parameter coordinates, with no
#' constraint projection applied. Applying the Jacobian of a parametrisation
#' belongs to whatever owns that parametrisation.
#'
#' @param family A [sampling_family].
#' @param theta Parameter vector.
#' @param x `(M, K)` matrix of outcomes.
#' @return `(M, d)` matrix.
#' @export
score <- new_generic("score", "family", function(family, theta, x) {
  S7::S7_dispatch()
})


#' Draw observations from `p_theta`
#'
#' Named `draw` rather than `simulate` to avoid masking [stats::simulate()].
#' @param family A [sampling_family].
#' @param theta Parameter vector.
#' @param n_obs Number of draws.
#' @return `(n_obs, k)` numeric matrix.
#' @export
draw <- new_generic("draw", "family", function(family, theta, n_obs) {
  S7::S7_dispatch()
})

#' Enumerated support of a finite family
#'
#' Every outcome with positive probability under some parameter. Families with
#' infinite sample spaces error: exact integration over the sample space is
#' undefined for them, not merely slow.
#' @param family A [sampling_family].
#' @return `(M, K)` matrix of outcomes.
#' @export
support <- new_generic("support", "family", function(family) S7::S7_dispatch())


method(support, sampling_family) <- function(family) {
  stop(
    "`",
    S7_class(family)@name,
    "` has no enumerable support. Use a Monte ",
    "Carlo or quadrature engine instead of an exact one.",
    call. = FALSE
  )
}


#' Does this family have a finite, enumerable sample space?
#' @param family A [sampling_family].
#' @return `TRUE` or `FALSE`.
#' @export
is_finite_support <- new_generic(
  "is_finite_support",
  "family",
  function(family) {
    S7::S7_dispatch()
  }
)

method(is_finite_support, sampling_family) <- function(family) FALSE

#' A reference point for the parameter space.
#'
#' Used as a fallback for initialising the atoms for a RIPr optimisation run,
#' where the alternative does not take the form of a mixture.
#' @keywords internal
#' @noRd
reference_parameter <- new_generic("reference_parameter", "family", \(family) {
  S7::S7_dispatch()
})


# --- Multinomial --------------------------------------------------------------

#' Log multinomial coefficient `log(n! / prod x_j!)` per count vector
#' @param x `(M, K)` matrix of count vectors.
#' @param n Total trials.
#' @keywords internal
#' @noRd
log_multinom_coef <- function(x, n) {
  lgamma(n + 1) - rowSums(lgamma(as.matrix(x) + 1))
}


#' Enumerate every count vector with `K` categories summing to `n`
#'
#' Stars and bars: each count vector is a choice of `K - 1` bar positions among
#' `n + K - 1` slots, so enumeration costs `O(M * K)`.
#'
#' @param n Total trials.
#' @param k Number of categories.
#' @return `(M, K)` integer matrix, `M = choose(n + K - 1, K - 1)`.
#' @keywords internal
#' @noRd
enumerate_counts <- function(n, k) {
  if (k == 1L) {
    return(matrix(as.integer(n), nrow = 1L))
  }
  bars <- utils::combn(n + k - 1L, k - 1L)
  t(diff(rbind(0L, bars, n + k)) - 1L)
}


#' Multinomial sampling family
#'
#' The `n`-trial, `K`-category multinomial. The support enumeration and its log
#' multinomial coefficients are built once at construction and stored as plain
#' data rather than closures, so a serialised family carries no environment.
#'
#' @param n_trials Integer. Trials per observation.
#' @param k Integer. Number of categories.
#' @return A `multinomial_family`.
#' @export
multinomial_family <- new_class(
  "multinomial_family",
  parent = sampling_family,
  properties = list(
    n_trials = class_numeric,
    k = class_numeric,
    outcomes = class_any,
    log_coef = class_numeric
  ),
  constructor = function(n_trials, k) {
    n_trials <- as.integer(n_trials)
    k <- as.integer(k)
    stopifnot(
      "`n_trials` must be a single non-negative integer" = length(n_trials) ==
        1L &&
        !is.na(n_trials) &&
        n_trials >= 0L,
      "`k` must be a single integer >= 1" = length(k) == 1L &&
        !is.na(k) &&
        k >= 1L
    )
    outcomes <- enumerate_counts(n_trials, k)
    new_object(
      S7_object(),
      n_trials = n_trials,
      k = k,
      outcomes = outcomes,
      log_coef = log_multinom_coef(outcomes, n_trials)
    )
  }
)


method(param_dim, multinomial_family) <- function(family) as.integer(family@k)


method(support, multinomial_family) <- function(family) family@outcomes


method(is_finite_support, multinomial_family) <- function(family) TRUE


method(compile_loglik, multinomial_family) <- function(family, x) {
  x <- as_outcome_matrix(x)
  log_coef <- log_multinom_coef(x, family@n_trials)

  function(theta_mat) {
    # matmul_0_ninf, not %*%: a zero category probability gives log(0) = -Inf,
    # and a zero count against it must contribute 0 rather than NaN.
    matmul_0_ninf(x, log(as.matrix(theta_mat))) + log_coef
  }
}


method(score, multinomial_family) <- function(family, theta, x) {
  nan_to_zero(div_by_col(as_outcome_matrix(x), theta))
}


method(draw, multinomial_family) <- function(family, theta, n_obs) {
  t(stats::rmultinom(n_obs, size = family@n_trials, prob = theta))
}


method(reference_parameter, multinomial_family) <- function(family) {
  rep(1 / family@k, family@k)
}
