#' @include family.R mixing.R distribution.R quadrature.R
NULL

#' Log multinomial coefficient `log(n! / prod x_j!)` per count vector
#' @param x `(M, K)` matrix of count vectors.
#' @param n Total trials.
#' @keywords internal
#' @noRd
log_multinom_coef <- function(x, n) {
  lgamma(n + 1) - rowSums(lgamma(as.matrix(x) + 1))
}


#' Multinomial sampling family
#'
#' The `n`-trial, `K`-category multinomial over [count_space()]. Its log
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
  parent = parametric_family,
  properties = list(
    n_trials = class_numeric,
    k = class_numeric
  ),
  constructor = function(n_trials, k) {
    space <- count_space(n = n_trials, k = k)
    new_object(
      at_theta,
      sample_space = space,
      # The full k-simplex, spanned by its vertices {e_1, ..., e_k}.
      parameter_space = simplex_region(vertices = diag(space@k)),
      n_trials = space@n,
      k = space@k
    )
  }
)


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


method(kernel_draw, multinomial_family) <- function(family, theta_mat) {
  # The conditional-binomial form:
  # `x_j | x_1..x_(j-1) ~ Binom(n - sum, p_j / (1 - sum p))`.
  # `rbinom()` is vectorised over both size and prob, so just `k - 1` calls.
  theta_mat <- as.matrix(theta_mat)
  m <- ncol(theta_mat)
  k <- nrow(theta_mat)
  out <- matrix(0L, nrow = m, ncol = k)
  remaining <- rep(as.integer(family@n_trials), m)
  unspent <- rep(1, m)
  for (j in seq_len(k - 1L)) {
    share <- ifelse(unspent > 0, pmin(1, pmax(0, theta_mat[j, ] / unspent)), 0)
    out[, j] <- stats::rbinom(m, remaining, share)
    remaining <- remaining - out[, j]
    unspent <- unspent - theta_mat[j, ]
  }
  out[, k] <- remaining
  out
}
