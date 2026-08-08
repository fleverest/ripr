#' @include family.R mixing_measure.R mixture.R quadrature.R
NULL


# Multivariate normal helpers. `chol_l` is the lower-triangular factor L with
# cov = L L^T; whitening solves L z = x - mean.

#' Log density of `N(mean, L L^T)` at each row of `x`
#' @keywords internal
#' @noRd
mvn_log_density <- function(x, mean, chol_l) {
  z <- forwardsolve(chol_l, t(x) - mean)
  -0.5 *
    length(mean) *
    log(2 * pi) -
    sum(log(diag(chol_l))) -
    0.5 * colSums(z^2)
}


#' Validate and factor a covariance matrix
#' @keywords internal
#' @noRd
as_covariance <- function(sigma, d, what = "sigma") {
  if (is.null(sigma)) {
    sigma <- diag(d)
  }
  sigma <- as.matrix(sigma)
  if (nrow(sigma) != d || ncol(sigma) != d) {
    stop("`", what, "` must be ", d, " by ", d, ".", call. = FALSE)
  }
  if (!isTRUE(all.equal(sigma, t(sigma)))) {
    stop("`", what, "` must be symmetric.", call. = FALSE)
  }
  if (any(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values <= 0)) {
    stop("`", what, "` must be positive definite.", call. = FALSE)
  }
  sigma
}


#' Gaussian sampling family with known covariance
#'
#' Observations are single draws \eqn{X \sim N(\theta, \Sigma)}{X ~ N(theta, sigma)}
#' with \eqn{\Sigma}{sigma} known, so the parameter is the mean.
#'
#' The sample space is continuous, so [support()] errors: pair this family with
#' [mc_engine()] or [gh_engine()] rather than [exact_engine()]. It is also a
#' family for which a certified gap bound is unavailable.
#'
#' @param dim Integer dimension of the observation.
#' @param sigma Known covariance matrix, or `NULL` for the identity.
#' @return A `gaussian_family`.
#' @export
gaussian_family <- new_class(
  "gaussian_family",
  parent = sampling_family,
  properties = list(
    n_dim = class_numeric,
    sigma = class_any,
    chol_l = class_any,
    sigma_inv = class_any
  ),
  constructor = function(dim, sigma = NULL) {
    dim <- as.integer(dim)
    stopifnot(
      "`dim` must be a single positive integer" = length(dim) == 1L &&
        !is.na(dim) &&
        dim >= 1L
    )
    sigma <- as_covariance(sigma, dim)
    new_object(
      S7_object(),
      n_dim = dim,
      sigma = sigma,
      chol_l = t(chol(sigma)),
      sigma_inv = chol2inv(chol(sigma))
    )
  }
)


method(param_dim, gaussian_family) <- function(family) as.integer(family@n_dim)


#' The sample space is all of `R^d`, so the inherited checks -- numeric, right
#' shape, no missing values -- are nearly the whole contract. All this adds is
#' finiteness: an infinite outcome has zero density under every parameter, so a
#' likelihood ratio there is `0 / 0`, and a `NaN` is a worse answer than a
#' complaint.
#' @keywords internal
#' @noRd
method(as_outcomes, gaussian_family) <- function(family, x) {
  x <- check_outcome_shape(x, outcome_dim(family))
  if (any(!is.finite(x))) {
    stop("Gaussian outcomes must be finite.", call. = FALSE)
  }
  x
}


method(compile_loglik, gaussian_family) <- function(family, x) {
  x <- as_outcome_matrix(x)
  # log p(x | theta) expands as
  #   -0.5 d log(2 pi) - log|L| - 0.5 x' S^-1 x  +  x' S^-1 theta  -  0.5 theta' S^-1 theta
  # whose first group depends only on x and whose second is one matrix multiply
  # over all parameter columns at once, replacing C triangular solves.
  x_sinv <- x %*% family@sigma_inv
  const <- -0.5 *
    family@n_dim *
    log(2 * pi) -
    sum(log(diag(family@chol_l))) -
    0.5 * rowSums(x_sinv * x)

  function(theta_mat) {
    theta_mat <- as.matrix(theta_mat)
    quad <- 0.5 * colSums(theta_mat * (family@sigma_inv %*% theta_mat))
    add_by_col(x_sinv %*% theta_mat, -quad) + const
  }
}


method(score, gaussian_family) <- function(family, theta, x) {
  t(family@sigma_inv %*% (t(as_outcome_matrix(x)) - theta))
}


method(draw, gaussian_family) <- function(family, theta, n_obs) {
  d <- as.integer(family@n_dim)
  z <- matrix(stats::rnorm(n_obs * d), nrow = d, ncol = n_obs)
  t(family@chol_l %*% z + theta)
}


#' Gaussian mixing measure over the mean of a Gaussian family
#'
#' A prior \eqn{\mu \sim N(m, V)}{mu ~ N(m, V)} over the mean of a
#' [gaussian_family()]. Paired with that family it induces the mixture
#' \eqn{N(m, \Sigma + V)}{N(m, sigma + V)} in closed form, so no numerical
#' integration is needed to evaluate it.
#'
#' @param prior_mean Numeric prior mean vector.
#' @param prior_cov Prior covariance, symmetric positive definite.
#' @return A `gaussian_mixing`.
#' @export
gaussian_mixing <- new_class(
  "gaussian_mixing",
  parent = mixing_measure,
  properties = list(prior_mean = class_numeric, prior_cov = class_any),
  constructor = function(prior_mean, prior_cov) {
    prior_mean <- as.numeric(prior_mean)
    new_object(
      S7_object(),
      prior_mean = prior_mean,
      prior_cov = as_covariance(prior_cov, length(prior_mean), "prior_cov")
    )
  }
)


method(n_atoms, gaussian_mixing) <- function(x) NA_integer_


method(induced_log_density, list(gaussian_mixing, gaussian_family)) <- function(
  mixing,
  family,
  x
) {
  mvn_log_density(
    as_outcome_matrix(x),
    mixing@prior_mean,
    t(chol(family@sigma + mixing@prior_cov))
  )
}


method(induced_draw, list(gaussian_mixing, gaussian_family)) <- function(
  mixing,
  family,
  n_obs
) {
  d <- length(mixing@prior_mean)
  chol_l <- t(chol(family@sigma + mixing@prior_cov))
  z <- matrix(stats::rnorm(n_obs * d), nrow = d, ncol = n_obs)
  t(chol_l %*% z + mixing@prior_mean)
}


# --- Mode and reference parameter for a Gaussian mixing measures --------------

method(mode_parameter, gaussian_mixing) <- function(x) x@prior_mean


method(reference_parameter, gaussian_family) <- function(family) {
  rep(0, family@n_dim)
}


# --- Moments, for Gauss-Hermite quadrature ------------------------------------

method(
  induced_gaussian_moments,
  list(point_mixing, gaussian_family)
) <- function(
  mixing,
  family
) {
  list(mean = mixing@theta_star, cov = family@sigma)
}


method(
  induced_gaussian_moments,
  list(gaussian_mixing, gaussian_family)
) <- function(
  mixing,
  family
) {
  list(mean = mixing@prior_mean, cov = family@sigma + mixing@prior_cov)
}
