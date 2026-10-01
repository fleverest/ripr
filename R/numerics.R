# Internal numerics: log-space reductions and an -Inf-safe matrix product.
#
# Dimension convention, used throughout the package: every collection of points
# is a matrix with one point per row (outcomes and quadrature nodes `(M, K)`,
# atoms `(C, d)`, vertices and generators `(V, d)`, chart coordinates
# `(N, n_par)`). A bare vector is one point, i.e. one row (`as_row_matrix()`).
#
#   M | #outcomes, or quadrature nodes
#   C | mixture components (#atoms), or parameter points in a batch
#   K | outcome dimension (for multinomial, #categories)
#   d | parameter dimension
#
# The exception is the log-density matrix `compile_loglik(family, x)(theta_mat)`,
# which is `(M, C)`: its columns are indexed by the rows of `theta_mat`.

#' Replace NaN entries with zero
#'
#' `0 / 0` arises where a zero count meets a zero probability; that outcome has
#' probability zero, so its contribution is zero, not undefined.
#' @keywords internal
#' @noRd
nan_to_zero <- function(x) {
  if (!anyNA(x)) {
    return(x)
  }
  x[is.nan(x)] <- 0
  x
}


#' Offset column `j` of `mat` by `w[j]`
#'
#' On an `(M, C)` log-density matrix this adds a per-parameter constant; on a
#' point matrix it translates every point by `w`.
#' @keywords internal
#' @noRd
add_by_col <- function(mat, w) {
  stopifnot(
    "`w` must have one entry per column of `mat`" = length(w) == ncol(mat)
  )
  mat + rep(w, each = nrow(mat))
}


#' Divide column `j` of `mat` by `w[j]`
#' @keywords internal
#' @noRd
div_by_col <- function(mat, w) {
  stopifnot(
    "`w` must have one entry per column of `mat`" = length(w) == ncol(mat)
  )
  mat / rep(w, each = nrow(mat))
}


#' Log-sum-exp of each row of an `(M, C)` matrix; length `M`
#' @keywords internal
#' @noRd
row_logsumexp <- function(mat) {
  matrixStats::rowLogSumExps(as.matrix(mat))
}


#' Log-sum-exp of each column of an `(M, C)` matrix; length `C`
#' @keywords internal
#' @noRd
col_logsumexp <- function(mat) {
  matrixStats::colLogSumExps(as.matrix(mat))
}


#' `tcrossprod()` treating `0 * -Inf` as `0`
#'
#' `(M, K)` non-negative `a` by `(C, K)` `b` gives `(M, C)`. A zero weight on a
#' `-Inf` log-probability contributes zero; a positive weight still gives
#' `-Inf`. The second product that finds those entries is skipped when `b` has
#' no `-Inf`.
#' @keywords internal
#' @noRd
tcrossprod_0_ninf <- function(a, b) {
  a <- as_row_matrix(a)
  b <- as_row_matrix(b)
  neg_inf <- is.infinite(b) & b < 0
  if (!any(neg_inf)) {
    return(tcrossprod(a, b))
  }
  b_safe <- b
  b_safe[neg_inf] <- 0
  out <- tcrossprod(a, b_safe)
  out[tcrossprod(a != 0, neg_inf) > 0] <- -Inf
  out
}


#' Coerce a point, or points, to a one-point-per-row matrix
#' @keywords internal
#' @noRd
as_row_matrix <- function(x) {
  if (is.null(dim(x))) matrix(x, nrow = 1L) else as.matrix(x)
}
