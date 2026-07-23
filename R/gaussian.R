# Multivariate normal helpers shared by the Gaussian family and its closed-form
# alternative. `chol_l` is the lower-triangular Cholesky factor L with
# cov = L L^T; whitening solves L z = x - mean.

mvn_log_density_matrix <- function(x_mat, mean, chol_l) {
  d <- length(mean)
  z <- forwardsolve(chol_l, t(x_mat) - mean)
  -0.5 * d * log(2 * pi) -
    sum(log(diag(chol_l))) -
    0.5 * colSums(z^2)
}

as_outcome_matrix <- function(x) {
  if (is.null(dim(x))) {
    matrix(x, nrow = 1L)
  } else {
    as.matrix(x)
  }
}

#' Gaussian sampling family with known covariance
#'
#' Observations are single draws `X ~ N(theta, sigma)` with `sigma` known. The
#' covariance is Cholesky-whitened once at construction, so densities and scores
#' cost one triangular solve; with the default `sigma = NULL` (identity) the
#' family is already white.
#'
#' The sample space is continuous, so [support()] signals an error: pair this
#' family with a Monte Carlo engine.
#'
#' @param dim Integer dimension of the observation.
#' @param sigma Known covariance matrix, or `NULL` for the identity.
#' @return A `gaussian_family`.
#' @export
gaussian_family <- new_class(
  "gaussian_family",
  parent = family,
  properties = list(
    n_dim = class_numeric,
    sigma = class_any,
    chol_l = class_any,
    sigma_inv = class_any
  ),
  constructor = function(dim, sigma = NULL) {
    dim <- as.integer(dim)
    if (is.null(sigma)) {
      sigma <- diag(dim)
    }
    if (
      !isTRUE(all.equal(sigma, t(sigma))) ||
        any(eigen(sigma, only.values = TRUE)$values <= 0)
    ) {
      stop("`sigma` must be a symmetric positive-definite matrix")
    }
    chol_l <- t(chol(sigma))
    new_object(
      S7_object(),
      n_dim = as.numeric(dim),
      sigma = sigma,
      chol_l = chol_l,
      sigma_inv = chol2inv(chol(sigma))
    )
  }
)

method(log_density, gaussian_family) <- function(family, theta, x = NULL) {
  if (is.null(x)) {
    stop(
      "gaussian_family has no finite support; supply outcomes `x` ",
      "(engines pass their draw set automatically)"
    )
  }
  mvn_log_density_matrix(as_outcome_matrix(x), theta, family@chol_l)
}

method(log_density_batch, gaussian_family) <- function(family, theta_mat, x = NULL) {
  if (is.null(x)) {
    stop(
      "gaussian_family has no finite support; supply outcomes `x` ",
      "(engines pass their draw set automatically)"
    )
  }
  x_mat <- as_outcome_matrix(x)
  cols <- lapply(seq_len(ncol(theta_mat)), function(j) {
    mvn_log_density_matrix(x_mat, theta_mat[, j], family@chol_l)
  })
  do.call(cbind, cols)
}

method(score, gaussian_family) <- function(family, theta, x = NULL) {
  if (is.null(x)) {
    stop(
      "gaussian_family has no finite support; supply outcomes `x` ",
      "(engines pass their draw set automatically)"
    )
  }
  x_mat <- as_outcome_matrix(x)
  t(family@sigma_inv %*% (t(x_mat) - theta))
}

method(param_dim, gaussian_family) <- function(family) {
  as.integer(family@n_dim)
}

method(simulate, gaussian_family) <- function(family, theta, n_obs, seed = NULL) {
  with_rng_seed(seed, {
    d <- as.integer(family@n_dim)
    z <- matrix(rnorm(n_obs * d), nrow = d, ncol = n_obs)
    t(family@chol_l %*% z + theta)
  })
}

method(support, gaussian_family) <- function(family) {
  stop(
    "gaussian_family has a continuous sample space with no finite support; ",
    "use a Monte Carlo engine (mc_engine) instead of an exact engine"
  )
}

#' Gaussian-prior alternative: Q-bar with inflated covariance (closed form)
#'
#' A Gaussian prior `mu ~ N(prior_mean, prior_cov)` over the mean of a
#' [gaussian_family()] induces the marginal outcome distribution
#' `Q-bar = N(prior_mean, sigma + prior_cov)` in closed form. The density uses
#' the inflated covariance directly, and the sampler draws `mu` then `X | mu`,
#' both seeded.
#'
#' @param prior_mean Numeric prior mean vector.
#' @param prior_cov Prior covariance (symmetric positive definite).
#' @return A `gaussian_prior_alt`.
#' @export
gaussian_prior_alt <- new_class(
  "gaussian_prior_alt",
  parent = alternative,
  properties = list(
    prior_mean = class_numeric,
    prior_cov = class_any
  ),
  validator = function(self) {
    if (!is.matrix(self@prior_cov)) {
      return("`prior_cov` must be a matrix")
    }
    if (nrow(self@prior_cov) != length(self@prior_mean)) {
      return("`prior_cov` dimensions must match `prior_mean`")
    }
    NULL
  }
)

method(q_log_density, list(gaussian_prior_alt, gaussian_family)) <- function(
  alt,
  family,
  x = NULL
) {
  if (is.null(x)) {
    stop("gaussian alternatives need explicit outcomes `x`")
  }
  marginal_chol <- t(chol(family@sigma + alt@prior_cov))
  mvn_log_density_matrix(as_outcome_matrix(x), alt@prior_mean, marginal_chol)
}

method(q_sample, list(gaussian_prior_alt, gaussian_family)) <- function(
  alt,
  family,
  n_obs,
  seed = NULL
) {
  with_rng_seed(seed, {
    d <- length(alt@prior_mean)
    marginal_chol <- t(chol(family@sigma + alt@prior_cov))
    z <- matrix(rnorm(n_obs * d), nrow = d, ncol = n_obs)
    t(marginal_chol %*% z + alt@prior_mean)
  })
}
