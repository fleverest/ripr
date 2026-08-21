#' @include family.R mixing_measure.R mixture.R quadrature.R
NULL

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
#' @examples
#' multinomial_family(n_trials = 20L, k = 3L)
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


method(as_outcomes, multinomial_family) <- function(family, x) {
  x <- check_outcome_shape(x, outcome_dim(family))
  if (any(x < 0) || any(x != trunc(x))) {
    stop(
      "multinomial outcomes must be non-negative whole numbers.",
      call. = FALSE
    )
  }
  totals <- rowSums(x)
  if (any(totals != family@n_trials)) {
    stop(
      "multinomial outcomes must sum to `n_trials` (",
      family@n_trials,
      "); got ",
      paste(unique(totals[totals != family@n_trials]), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  x
}


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
