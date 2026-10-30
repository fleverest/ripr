#' @include family.R mixing.R distribution.R quadrature.R space.R dc_enclosure.R
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


# --- Components, for Gauss-Hermite quadrature ---------------------------------

method(
  mixture_gaussian_components,
  list(finite_dist, gaussian_family)
) <- function(mixing, family) {
  lapply(seq_len(nrow(mixing@atoms)), function(i) {
    list(
      mean = mixing@atoms[i, ],
      cov = family@sigma,
      weight = mixing@weights[[i]]
    )
  })
}


method(
  mixture_gaussian_components,
  list(gaussian_dist, gaussian_family)
) <- function(mixing, family) {
  list(list(mean = mixing@mean, cov = family@sigma + mixing@cov, weight = 1))
}


# The Gaussian difference-of-convex (DC) enclosure
#
# The DC branch and bound (see `dc_enclosure.R`) tiles the null with simplices
# and bounds g(theta) = E_theta[X] over each one. It does so using two things
# the family supplies: (1) a bound on the curvature of log g, and (2) bounds on
# g at the locations theta in the null.
#
# The points theta at which we need g are these vertices: they are the
# corners of the initial triangulation, plus the midpoints added during a split.
# A vertex's value is bounded once and shared by every cell that meets there.
# This section supplies those values for the Gaussian family, through
# `log_expectation_bounds()`.
#
# Write P_theta = N(theta, Sigma) for the null distributions,
# Q = sum_i pi_i N(mu_i, S_i) and P = N(nu, C), and X = q / p. Bounding g at a
# vertex is easy when P is a single Gaussian: for each i, the logs of q_i(z),
# p_theta(z) and p(z) are quadratics in z, so q_i(z) p_theta(z) / p(z) is the
# exponential of a quadratic. Completing the square makes it a constant times a
# Gaussian density, and the integral is that constant. No quadrature
# approximation is required for such an integral, and the lower and upper bounds
# therefore coincide (up to rounding):
#
#   g(theta) = sum_i pi_i |Sigma|^-1/2 |S_i|^-1/2 |C|^1/2 |Gamma_i|^-1/2
#              exp(gamma_i' Gamma_i^-1 gamma_i / 2 - k_i / 2),
#   Gamma_i = Sigma^-1 + S_i^-1 - C^-1,
#   gamma_i = Sigma^-1 theta + S_i^-1 mu_i - C^-1 nu,
#   k_i = theta' Sigma^-1 theta + mu_i' S_i^-1 mu_i - nu' C^-1 nu.
#
# Completing the square needs Gamma_i positive definite. When it is not, p
# has lighter tails than q_i p_theta in some direction, and the integral
# does not exist, so g(theta) is infinite. With no P (i.e. X is a likelihood
# X = q), we drop the C term and replace |C|^1/2 by (2 pi)^(-d/2) for the
# reference measure instead.
#
# A P with several atoms has no closed form at the vertices, because 1 / p
# is then not the exponential of a quadratic.

method(log_expectation_bounds, list(gaussian_family, likelihood_ratio)) <-
  function(family, x) {
    ratio <- gaussian_ratio(x, family)
    if (is.null(ratio)) {
      return(paste0(
        "`X = Q / P` (or `Q`) does not have `Q` a Gaussian or finite Gaussian ",
        "mixture and `P` a single Gaussian, so its expectation has no closed ",
        "form at the vertices. A `P` with several atoms needs bounds on each ",
        "vertex's expectation, which are not implemented"
      ))
    }
    sigma_inv <- family@sigma_inv
    terms <- gaussian_terms(sigma_inv, log_det(family@sigma), ratio)
    function(theta) {
      v <- gaussian_log_g(sigma_inv, terms, theta)
      cbind(v, v, deparse.level = 0)
    }
  }


# --- Recognising X ------------------------------------------------------------

#' The Gaussian components of `Q` and `P` in the likelihood ratio
#' `X = Q / P` (or `X = Q`), or `NULL` if they are not of that form for this
#' family's dimension
#' @keywords internal
#' @noRd
gaussian_ratio <- function(x, family) {
  q <- gaussian_components_of(x@numerator, family)
  if (is.null(q)) {
    return(NULL)
  }
  if (is.null(x@denominator)) {
    return(list(q = q, p = NULL))
  }
  p <- gaussian_components_of(x@denominator, family)
  if (length(p) != 1L) {
    return(NULL)
  }
  list(q = q, p = p[[1L]])
}


#' @keywords internal
#' @noRd
gaussian_components_of <- function(dist, family) {
  if (
    !S7_inherits(dist, mixture) ||
      !S7_inherits(dist@family, gaussian_family) ||
      dist@family@d != family@d
  ) {
    return(NULL)
  }
  mixture_gaussian_components(dist@mixing, dist@family)
}


# --- The closed form ----------------------------------------------------------

#' @keywords internal
#' @noRd
log_det <- function(m) as.numeric(determinant(m, logarithm = TRUE)$modulus)


#' Per-component constants of `log g`, from the formula above; `NULL` marks a
#' component whose integral diverges
#' @keywords internal
#' @noRd
gaussian_terms <- function(sigma_inv, log_det_sigma, ratio) {
  d <- nrow(sigma_inv)
  p <- ratio$p
  if (is.null(p)) {
    C_inv <- matrix(0, d, d)
    C_inv_nu <- numeric(d)
    nu_C_inv_nu <- 0
    half_log_det_C <- -0.5 * d * log(2 * pi)
  } else {
    C_inv <- chol2inv(chol(p$cov))
    C_inv_nu <- drop(C_inv %*% p$mean)
    nu_C_inv_nu <- sum(p$mean * C_inv_nu)
    half_log_det_C <- 0.5 * log_det(p$cov)
  }
  lapply(ratio$q, function(comp) {
    S_inv <- chol2inv(chol(comp$cov))
    Gamma_root <- tryCatch(
      chol(sigma_inv + S_inv - C_inv),
      error = function(e) NULL
    )
    if (is.null(Gamma_root)) {
      return(NULL)
    }
    S_inv_mu <- drop(S_inv %*% comp$mean)
    list(
      Gamma_inv = chol2inv(Gamma_root),
      gamma0 = S_inv_mu - C_inv_nu,
      const = log(comp$weight) -
        0.5 * log_det_sigma -
        0.5 * log_det(comp$cov) +
        half_log_det_C -
        sum(log(diag(Gamma_root))) -
        0.5 * (sum(comp$mean * S_inv_mu) - nu_C_inv_nu)
    )
  })
}


#' `log E_theta[X]` at each row of `theta`, exactly
#' @keywords internal
#' @noRd
gaussian_log_g <- function(sigma_inv, terms, theta) {
  theta <- as_row_matrix(theta)
  sigma_inv_theta <- theta %*% sigma_inv
  half_quad_theta <- -0.5 * rowSums(sigma_inv_theta * theta)
  per_term <- vapply(
    terms,
    function(t) {
      if (is.null(t)) {
        return(rep(Inf, nrow(theta)))
      }
      gamma_i <- add_by_col(sigma_inv_theta, t$gamma0)
      t$const +
        0.5 * rowSums((gamma_i %*% t$Gamma_inv) * gamma_i) +
        half_quad_theta
    },
    numeric(nrow(theta))
  )
  row_logsumexp(matrix(per_term, nrow = nrow(theta)))
}
