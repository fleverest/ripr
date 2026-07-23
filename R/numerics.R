# Internal base R numerics shared across the package. Nothing here is
# RIPr-specific: these are the log-space reductions and the -Inf-safe matrix
# multiply that the torch-free rewrite leans on. Kept unexported.

#' Replace NaN entries with zero
#' @keywords internal
#' @noRd
nan_to_zero <- function(x) {
  x[is.nan(x)] <- 0
  x
}

#' Replace NaN entries with -Inf (log of zero)
#' @keywords internal
#' @noRd
nan_to_neginf <- function(x) {
  x[is.nan(x)] <- -Inf
  x
}

#' Numerically stable log-sum-exp of a vector
#'
#' Returns `-Inf` when every entry is `-Inf`.
#' @keywords internal
#' @noRd
logsumexp_vec <- function(v) {
  m <- max(v)
  if (!is.finite(m)) {
    return(m)
  }
  m + log(sum(exp(v - m)))
}

#' Column-wise log-sum-exp: reduce each column over its rows
#'
#' `mat` is `(M, C)`; the result is a length-`C` vector with one log-sum-exp
#' per column (the base R analogue of `logsumexp(mat, dim = 1)`).
#' @keywords internal
#' @noRd
col_logsumexp <- function(mat) {
  mat <- as.matrix(mat)
  m <- apply(mat, 2L, max)
  m_safe <- ifelse(is.finite(m), m, 0)
  out <- m_safe + log(colSums(exp(sweep(mat, 2L, m_safe, "-"))))
  out[!is.finite(m)] <- m[!is.finite(m)]
  out
}

#' Row-wise log-sum-exp: reduce each row over its columns
#'
#' `mat` is `(N, M)`; the result is a length-`N` vector with one log-sum-exp
#' per row (the base R analogue of `logsumexp(mat, dim = 2)`).
#' @keywords internal
#' @noRd
row_logsumexp <- function(mat) {
  mat <- as.matrix(mat)
  m <- apply(mat, 1L, max)
  m_safe <- ifelse(is.finite(m), m, 0)
  out <- m_safe + log(rowSums(exp(sweep(mat, 1L, m_safe, "-"))))
  out[!is.finite(m)] <- m[!is.finite(m)]
  out
}

#' Matrix multiplication treating `0 * -Inf` as `0`
#'
#' Standard `%*%` yields `NaN` for any `0 * -Inf` product. This variant zeroes
#' those contributions instead (the correct behaviour when a zero count means
#' "this log-probability is never used"), while still propagating `-Inf` into
#' any entry where a strictly positive weight meets a `-Inf`.
#'
#' @param A `(M, K)` numeric matrix (or vector), assumed non-negative.
#' @param B `(K, N)` numeric matrix (or vector), may contain `-Inf` entries.
#' @return `(M, N)` numeric matrix.
#' @keywords internal
#' @noRd
matmul_0_ninf <- function(A, B) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  neg_inf <- is.infinite(B) & B < 0
  b_safe <- B
  b_safe[neg_inf] <- 0
  result <- A %*% b_safe
  ninf_hits <- (A != 0) %*% neg_inf
  result[ninf_hits > 0] <- -Inf
  result
}

#' Coerce an outcome argument to an `(N, d)` matrix
#'
#' A bare length-d vector is treated as a single outcome (one row).
#' @keywords internal
#' @noRd
as_outcome_matrix <- function(x) {
  if (is.null(dim(x))) {
    matrix(x, nrow = 1L)
  } else {
    as.matrix(x)
  }
}
