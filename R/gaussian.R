#' @include family.R mixing.R distribution.R quadrature.R
NULL


#' Validate a covariance matrix
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
#' The sample space [real_region()] cannot be enumerated: use [mc_engine()] or
#' [gh_engine()], not [exact_engine()]. No certified gap bound is available at
#' the time of writing.
#'
#' @param d Integer dimension of the observation.
#' @param sigma Known covariance matrix, or `NULL` for the identity.
#' @return A `gaussian_family`.
#' @examples
#' gaussian_family(d = 2L)
#' gaussian_family(d = 2L, sigma = diag(c(1, 4)))
#' @export
gaussian_family <- new_class(
  "gaussian_family",
  parent = parametric_family,
  properties = list(
    d = class_numeric,
    sigma = class_any,
    sigma_inv = class_any
  ),
  constructor = function(d, sigma = NULL) {
    d <- as.integer(d)
    stopifnot(
      "`d` must be a single positive integer" = length(d) == 1L &&
        !is.na(d) &&
        d >= 1L
    )
    sigma <- as_covariance(sigma, d)
    new_object(
      at_theta,
      sample_space = real_region(d),
      parameter_space = real_region(d),
      d = d,
      sigma = sigma,
      sigma_inv = chol2inv(chol(sigma))
    )
  }
)


method(compile_loglik, gaussian_family) <- function(family, x) {
  x <- as_row_matrix(x)
  # log p(x | theta) = [terms in x only] + x' S^-1 theta - 0.5 theta' S^-1 theta.
  x_sinv <- x %*% family@sigma_inv
  const <- -0.5 *
    family@d *
    log(2 * pi) -
    0.5 * as.numeric(determinant(family@sigma)$modulus) -
    0.5 * rowSums(x_sinv * x)

  function(theta_mat) {
    theta_mat <- as_row_matrix(theta_mat)
    quad <- 0.5 * rowSums((theta_mat %*% family@sigma_inv) * theta_mat)
    add_by_col(tcrossprod(x_sinv, theta_mat), -quad) + const
  }
}


method(score, gaussian_family) <- function(family, theta, x) {
  # Row form of `sigma_inv (x - theta)`; valid as `sigma_inv` is symmetric.
  add_by_col(as_row_matrix(x), -theta) %*% family@sigma_inv
}


method(kernel_draw, gaussian_family) <- function(family, theta_mat) {
  theta_mat <- as_row_matrix(theta_mat)
  theta_mat +
    mvtnorm::rmvnorm(nrow(theta_mat), sigma = family@sigma, method = "chol")
}


#' Multivariate Gaussian distribution
#'
#' \eqn{N(m, V)}{N(m, V)} on \eqn{\mathbb{R}^d}{R^d}. As a mixing measure over
#' the mean of a [gaussian_family()] with covariance \eqn{\Sigma}{sigma}, it
#' induces the mixture \eqn{N(m, \Sigma + V)}{N(m, sigma + V)} in closed form.
#'
#' @param mean Numeric mean vector.
#' @param cov Covariance matrix, symmetric positive definite.
#' @return A `gaussian_dist`.
#' @examples
#' fam <- gaussian_family(d = 2L)
#' prior <- gaussian_dist(mean = c(0, 0), cov = diag(2))
#' log_density(fam(prior), c(0.5, 0.5))
#' @export
gaussian_dist <- new_class(
  "gaussian_dist",
  parent = continuous_dist,
  properties = list(
    mean = class_numeric,
    cov = class_any,
    sample_space = new_property(
      space,
      getter = function(self) real_region(length(self@mean))
    )
  ),
  constructor = function(mean, cov) {
    mean <- as.numeric(mean)
    new_object(
      S7_object(),
      mean = mean,
      cov = as_covariance(cov, length(mean), "cov")
    )
  }
)


#' @rdname gaussian_dist
#' @usage NULL
method(format, gaussian_dist) <- function(x, ...) {
  sprintf("gaussian_dist: mean %s", theta_label(x@mean))
}


#' @description `print()` shows the mean and covariance, summarising the
#'   covariance by its size above eight dimensions.
#' @rdname gaussian_dist
#' @usage NULL
method(print, gaussian_dist) <- function(x, ...) {
  d <- length(x@mean)
  cat("<gaussian_dist>\n")
  cat("  mean ", theta_label(x@mean), "\n", sep = "")
  if (d <= 8L) {
    cat("  covariance:\n")
    print(signif(x@cov, 4L))
  } else {
    cat("  covariance ", d, " x ", d, " matrix\n", sep = "")
  }
  invisible(x)
}


method(mixture_log_density, list(gaussian_dist, gaussian_family)) <- function(
  mixing,
  family,
  x
) {
  mvtnorm::dmvnorm(
    as_row_matrix(x),
    mixing@mean,
    family@sigma + mixing@cov,
    log = TRUE
  )
}


method(mixture_draw, list(gaussian_dist, gaussian_family)) <- function(
  mixing,
  family,
  n
) {
  mvtnorm::rmvnorm(n, mixing@mean, family@sigma + mixing@cov, method = "chol")
}


#' @rdname draw
#' @usage NULL
method(draw, gaussian_dist) <- function(dist, n) {
  mvtnorm::rmvnorm(n, dist@mean, dist@cov, method = "chol")
}


method(log_density, gaussian_dist) <- function(dist, x) {
  mvtnorm::dmvnorm(as_row_matrix(x), dist@mean, dist@cov, log = TRUE)
}


# --- Mode and reference parameter for a Gaussian mixing measures --------------

method(reference_point, gaussian_dist) <- function(x) x@mean


# --- Moments, for Gauss-Hermite quadrature ------------------------------------

method(
  mixture_gaussian_moments,
  list(finite_dist, gaussian_family)
) <- function(
  mixing,
  family
) {
  # Only a point mass leaves the mixture Gaussian.
  if (nrow(mixing@atoms) != 1L) {
    return(NULL)
  }
  list(mean = mixing@atoms[1L, ], cov = family@sigma)
}


method(
  mixture_gaussian_moments,
  list(gaussian_dist, gaussian_family)
) <- function(
  mixing,
  family
) {
  list(mean = mixing@mean, cov = family@sigma + mixing@cov)
}
