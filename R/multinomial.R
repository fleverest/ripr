#' @include family.R mixing.R distribution.R quadrature.R
NULL

#' Log multinomial coefficient `log(n! / prod x_j!)` per row of `x`.
#' @keywords internal
#' @noRd
log_multinom_coef <- function(x, n) {
  lgamma(n + 1) - rowSums(lgamma(as.matrix(x) + 1))
}


#' Multinomial sampling family
#'
#' The `n_trials`-trial, `k`-category multinomial over [count_space()], with
#' parameter space equal to the full probability simplex.
#'
#' @param n_trials Integer. Trials per observation.
#' @param k Integer. Number of categories.
#' @return A `multinomial_family`.
#' @examples
#' multinomial_family(n_trials = 20L, k = 3L)
#' @export
multinomial_family <- new_class(
  "multinomial_family",
  parent = parametric_family,
  properties = list(
    n_trials = class_numeric,
    k = class_numeric
  ),
  constructor = function(n_trials, k) {
    space <- count_space(n_trials = n_trials, k = k)
    new_object(
      at_theta,
      sample_space = space,
      parameter_space = simplex_region(vertices = diag(space@k)),
      n_trials = space@n_trials,
      k = space@k
    )
  }
)


method(compile_loglik, multinomial_family) <- function(family, x) {
  x <- as_row_matrix(x)
  log_coef <- log_multinom_coef(x, family@n_trials)

  function(theta_mat) {
    # Not tcrossprod because a zero count of a category with probability zero
    # contributes nothing to the likelihood. Otherwise we would get NaN.
    tcrossprod_0_ninf(x, log(as_row_matrix(theta_mat))) + log_coef
  }
}


method(score, multinomial_family) <- function(family, theta, x) {
  nan_to_zero(div_by_col(as_row_matrix(x), theta))
}


method(kernel_draw, multinomial_family) <- function(family, theta_mat) {
  # Conditional binomials, `x_j | x_<j ~ Binom(n - sum x_<j, p_j / (1 - sum
  # p_<j))`, vectorised over rows so only `k - 1` `rbinom()` calls.
  theta_mat <- as_row_matrix(theta_mat)
  m <- nrow(theta_mat)
  k <- ncol(theta_mat)
  out <- matrix(0L, nrow = m, ncol = k)
  remaining <- rep(as.integer(family@n_trials), m)
  unspent <- rep(1, m)
  for (j in seq_len(k - 1L)) {
    share <- ifelse(unspent > 0, pmin(1, pmax(0, theta_mat[, j] / unspent)), 0)
    out[, j] <- stats::rbinom(m, remaining, share)
    remaining <- remaining - out[, j]
    unspent <- unspent - theta_mat[, j]
  }
  out[, k] <- remaining
  out
}
