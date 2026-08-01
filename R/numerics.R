# Internal numerics. Nothing here is RIPr-specific: log-space reductions and an
# -Inf-safe matrix multiplication. Everything here works in log space.
#
# Dimension convention, used consistently throughout this file and the ones
# that consume it:
#
#   M  rows    | #outcomes, or quadrature nodes
#   C  columns | mixture components (#atoms)
#   K          | parameter dimension (for multinomial, #categories)
#
# So a log-density matrix is (M, C): one row per outcome, one column per atom.

#' Replace NaN entries with zero
#'
#' `0 / 0` arises wherever a zero count meets a zero category probability. The
#' outcome has probability zero under that parameter, so its contribution is
#' zero, not undefined.
#' @keywords internal
#' @noRd
nan_to_zero <- function(x) {
  if (!anyNA(x)) {
    return(x)
  }
  x[is.nan(x)] <- 0
  x
}


#' Replace NaN entries with -Inf (log of zero)
#' @keywords internal
#' @noRd
nan_to_neginf <- function(x) {
  if (!anyNA(x)) {
    return(x)
  }
  x[is.nan(x)] <- -Inf
  x
}


#' Offset each column of an `(M, C)` matrix by the matching entry of `w`
#'
#' Adds the scalar `w[j]` to every entry of column `j`. The column-wise
#' analogue of plain recycling (which already handles the row-wise case as
#' `mat + v`), without `sweep`.
#' @param mat `(M, C)` numeric matrix.
#' @param w Length-`C` numeric vector, one offset per column.
#' @return `(M, C)` numeric matrix.
#' @keywords internal
#' @noRd
add_by_col <- function(mat, w) {
  mat + rep(w, each = nrow(mat))
}


#' Scale each column of an `(M, C)` matrix by the matching entry of `w`
#'
#' Divides every entry of column `j` by the scalar `w[j]`.
#' @param mat `(M, C)` numeric matrix.
#' @param w Length-`C` numeric vector, one divisor per column.
#' @return `(M, C)` numeric matrix.
#' @keywords internal
#' @noRd
div_by_col <- function(mat, w) {
  mat / rep(w, each = nrow(mat))
}


#' Numerically stable log-sum-exp of a vector
#'
#' Returns `-Inf` when every entry is `-Inf`.
#' @keywords internal
#' @noRd
logsumexp_vec <- function(v) {
  matrixStats::logSumExp(v)
}


#' Row-wise log-sum-exp: reduce each row over its columns
#'
#' `mat` is `(M, C)`; returns a length-`M` vector, one log-sum-exp per row.
#' @keywords internal
#' @noRd
row_logsumexp <- function(mat) {
  matrixStats::rowLogSumExps(as.matrix(mat))
}


#' Column-wise log-sum-exp: reduce each column over its rows
#'
#' `mat` is `(M, C)`; returns a length-`C` vector, one log-sum-exp per column.
#' @keywords internal
#' @noRd
col_logsumexp <- function(mat) {
  matrixStats::colLogSumExps(as.matrix(mat))
}


#' Log-sum-exp of `v + log_w`, i.e. `log(sum_i w_i exp(v_i))`
#'
#' @param v Numeric vector of log values.
#' @param log_w Numeric vector of log weights, same length.
#' @keywords internal
#' @noRd
logsumexp_weighted <- function(v, log_w) {
  logsumexp_vec(v + log_w)
}


#' Matrix multiplication treating `0 * -Inf` as `0`
#'
#' Standard `%*%` yields `NaN` for any `0 * -Inf` product. This variant zeroes
#' those contributions (correct when a zero count means "this log-probability is
#' never used") while still propagating `-Inf` wherever a strictly positive
#' weight meets a `-Inf`.
#'
#' The `-Inf` bookkeeping needs a second matrix multiply to find which output
#' entries are poisoned, so it is skipped entirely when `b` has no `-Inf` at
#' all. This is the case in most scenarios, so the guard halves the cost of
#' this function, which is used very frequently during optimisation.
#'
#' @param a `(M, K)` numeric matrix (or vector), assumed non-negative.
#' @param b `(K, C)` numeric matrix (or vector), may contain `-Inf`.
#' @return `(M, C)` numeric matrix.
#' @keywords internal
#' @noRd
matmul_0_ninf <- function(a, b) {
  a <- as.matrix(a)
  b <- as.matrix(b)
  neg_inf <- is.infinite(b) & b < 0
  if (!any(neg_inf)) {
    return(a %*% b)
  }
  b_safe <- b
  b_safe[neg_inf] <- 0
  out <- a %*% b_safe
  out[(a != 0) %*% neg_inf > 0] <- -Inf
  out
}


#' Coerce an outcome argument to an `(M, K)` matrix
#'
#' A bare length-`K` vector is a single outcome, i.e. one row.
#' @keywords internal
#' @noRd
as_outcome_matrix <- function(x) {
  if (is.null(dim(x))) matrix(x, nrow = 1L) else as.matrix(x)
}


#' Is this a probability vector on the simplex?
#' @keywords internal
#' @noRd
is_prob_vector <- function(p, tol = 1e-9) {
  is.numeric(p) && !anyNA(p) && all(p >= 0) && abs(sum(p) - 1) <= tol
}
