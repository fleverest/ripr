# Multinomial likelihood arithmetic (base R). The support enumeration, the
# log multinomial coefficient, and the log-PMF matrix multiply the exact
# engine and the multinomial family are built on. All log-space; the
# `matmul_0_ninf` helper in numerics.R handles the `0 * -Inf` products that
# arise when a category probability is exactly zero.

#' Log multinomial coefficient `log(n! / prod x_i!)` per count vector
#'
#' @param X `(M, m)` matrix of count vectors.
#' @param n Total number of trials.
#' @return Length-`M` numeric vector.
#' @keywords internal
#' @noRd
log_multinom_coef <- function(X, n) {
  X <- as.matrix(X)
  lgamma(n + 1) - rowSums(lgamma(X + 1))
}

#' All valid count vectors for a multinomial with `m` categories
#'
#' Enumerates every integer vector `(x_1, ..., x_m)` with `x_i >= 0` and
#' `sum(x_i) = n` via stars and bars: each count vector corresponds to a choice
#' of `m - 1` bar positions among `n + m - 1` slots, so the enumeration costs
#' `O(M * m)` -- proportional to the output.
#'
#' @param n Total number of trials.
#' @param m Number of categories. Default: 3.
#' @return `(M, m)` integer matrix where `M = choose(n + m - 1, m - 1)`.
#' @keywords internal
#' @noRd
build_counts_matrix <- function(n, m = 3) {
  if (m == 1L) {
    return(matrix(as.integer(n), nrow = 1L))
  }
  bars <- combn(n + m - 1L, m - 1L) # (m - 1, M) increasing bar positions
  boundaries <- rbind(0L, bars, n + m)
  t(diff(boundaries) - 1L) # x_i = b_i - b_{i-1} - 1, with b_0 = 0, b_m = n + m
}

#' Log PMF of a multinomial for `M` count vectors and `C` component
#' distributions
#'
#' Computes `log Multinomial(x; n, p_c)` for every `(x, p_c)` pair via a single
#' `-Inf`-safe matrix multiplication in log space.
#'
#' @param X `(M, m)` matrix of count vectors from `build_counts_matrix()`.
#' @param log_comp_probs `(m, C)` matrix of log category probabilities for each
#'   of `C` component distributions.
#' @param n Total number of trials.
#' @return `(M, C)` matrix of log PMF values.
#' @keywords internal
#' @noRd
mnom_logpmf <- function(X, log_comp_probs, n) {
  X <- as.matrix(X)
  lc <- log_multinom_coef(X, n)
  sweep(matmul_0_ninf(X, log_comp_probs), 1L, lc, "+")
}

#' Multinomial likelihood interface
#'
#' Constructs the likelihood interface for the multinomial distribution with
#' `n` draws over `K` categories: the enumerated support, and closures for the
#' log-PMF (single and batched) and the score. All base R.
#'
#' @param n Integer. Total number of draws.
#' @param K Integer. Number of categories.
#' @return A named list with `support` (`(M, K)` matrix), `M`, `log_pmf`,
#'   `log_pmf_batch`, `score`, `n`, `K`.
#' @keywords internal
#' @noRd
make_multinomial_likelihood <- function(n, K) {
  support_mat <- build_counts_matrix(n, K)
  M <- nrow(support_mat)
  log_base <- log_multinom_coef(support_mat, n)

  log_pmf <- function(theta) {
    as.vector(matmul_0_ninf(support_mat, log(theta))) + log_base
  }

  log_pmf_batch <- function(theta_mat) {
    sweep(matmul_0_ninf(support_mat, log(theta_mat)), 1L, log_base, "+")
  }

  score <- function(theta) {
    sweep(support_mat, 2L, theta, "/") - n
  }

  list(
    support = support_mat,
    M = M,
    log_pmf = log_pmf,
    log_pmf_batch = log_pmf_batch,
    score = score,
    n = n,
    K = K
  )
}
