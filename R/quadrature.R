#' @include family.R distribution.R
NULL

#' Quadrature rule for expectations under the alternative
#'
#' Every quantity the optimiser touches is an expectation under \eqn{Q}{Q}: the
#' objective \eqn{KL(Q \| P_W)}{KL(Q || P_W)}, the linear oracle
#' \eqn{G(\theta) = E_Q[p_\theta / P_W]}{G(theta) = E_Q[p_theta / P_W]}, the EM
#' responsibilities, and the gap. So an engine is one thing: a set of nodes with
#' weights approximating \eqn{E_Q[\cdot]}{E_Q[.]}. Exact enumeration, Monte
#' Carlo, and deterministic quadrature differ only in how the nodes and weights
#' are chosen.
#'
#' `log_q` is held separately from `log_w` because the two coincide only for the
#' exact engine. Under Monte Carlo the weights are `1/M` while `log_q` is the
#' density at the drawn nodes, and the objective needs both.
#'
#' Build one with [resolve_engine()] rather than the raw constructor.
#'
#' @param nodes `(M, K)` matrix of evaluation points.
#' @param log_w Length-`M` log quadrature weights; `exp(log_w)` sums to 1.
#' @param log_q Length-`M` log density of the alternative at `nodes`.
#' @param family The [parametric_family] the nodes live over.
#' @param deterministic Is the rule free of sampling error?
#' @return A `quadrature`.
#' @examples
#' fam <- multinomial_family(n_trials = 3L, k = 2L)
#' Q <- induced_distribution(fam, point_mixing(c(0.5, 0.5)))
#' resolve_engine(exact_engine(), Q, fam)
#' @export
quadrature <- new_class(
  "quadrature",
  properties = list(
    nodes = class_any,
    log_w = class_numeric,
    log_q = class_numeric,
    family = parametric_family,
    deterministic = class_logical
  ),
  validator = function(self) {
    if (!is.matrix(self@nodes)) {
      return("`nodes` must be a matrix with one outcome per row")
    }
    if (nrow(self@nodes) != length(self@log_w)) {
      return("`log_w` needs one weight per node")
    }
    if (nrow(self@nodes) != length(self@log_q)) {
      return("`log_q` needs one log_q per node")
    }
    NULL
  }
)


#' Number of quadrature nodes
#' @param engine A [quadrature].
#' @return Integer.
#' @examples
#' fam <- multinomial_family(n_trials = 3L, k = 2L)
#' Q <- induced_distribution(fam, point_mixing(c(0.5, 0.5)))
#' n_nodes(resolve_engine(exact_engine(), Q, fam))
#' @export
n_nodes <- function(engine) nrow(engine@nodes)


#' Is this rule free of sampling error?
#'
#' `TRUE` for exact enumeration and deterministic quadrature, `FALSE` for Monte
#' Carlo. A deterministic rule may still carry approximation error, but it is
#' bias rather than variance, so [expect_se()] does not describe it.
#' @param engine A [quadrature].
#' @return `TRUE` or `FALSE`.
#' @examples
#' fam <- multinomial_family(n_trials = 3L, k = 2L)
#' Q <- induced_distribution(fam, point_mixing(c(0.5, 0.5)))
#' deterministic(resolve_engine(exact_engine(), Q, fam))
#' deterministic(resolve_engine(mc_engine(100L), Q, fam))
#' @export
deterministic <- function(engine) isTRUE(engine@deterministic)


#' Compile the family's log-likelihood at the engine's nodes
#'
#' Returns a function of `theta_mat` giving the `(M, C)` matrix of log densities
#' at the nodes. Compile once at the top of a step and reuse: the nodes do not
#' move, so recompiling per evaluation repeats work that depends on them alone.
#'
#' Engines hold their nodes as plain data and compile on demand rather than
#' storing the evaluator, so a serialised engine carries no captured
#' environment.
#' @param engine A [quadrature].
#' @return A function of `theta_mat`.
#' @examples
#' fam <- multinomial_family(n_trials = 3L, k = 2L)
#' Q <- induced_distribution(fam, point_mixing(c(0.5, 0.5)))
#' engine <- resolve_engine(exact_engine(), Q, fam)
#' ll <- compile_engine(engine)
#' ll(cbind(c(0.5, 0.5), c(0.2, 0.8)))
#' @export
compile_engine <- function(engine) {
  compile_loglik(engine@family, engine@nodes)
}


#' Expectation under the alternative
#'
#' \eqn{E_Q[v] \approx \sum_i w_i v_i}{E_Q[v] = sum_i w_i v_i}, for values given
#' directly. Use this when `v` may be negative -- the KL objective, for
#' instance, integrates a log ratio.
#' @param engine A [quadrature].
#' @param v Length-`M` numeric vector of integrand values at the nodes.
#' @return Numeric scalar.
#' @examples
#' fam <- multinomial_family(n_trials = 3L, k = 2L)
#' Q <- induced_distribution(fam, point_mixing(c(0.5, 0.5)))
#' engine <- resolve_engine(exact_engine(), Q, fam)
#' expect_q(engine, engine@nodes[, 1])
#' @export
expect_q <- function(engine, v) {
  sum(exp(engine@log_w) * v)
}


#' Expectation under the alternative, in log space
#'
#' \eqn{\log E_Q[v]}{log E_Q[v]} where `v >= 0` is supplied as `log_v`. Stable
#' where [expect_q()] would overflow, which it does routinely: likelihood ratios
#' at audit scale exceed the range of a double long before they stop being
#' meaningful.
#' @param engine A [quadrature].
#' @param log_v Length-`M` numeric vector of log integrand values at the nodes.
#' @return Numeric scalar.
#' @examples
#' fam <- multinomial_family(n_trials = 3L, k = 2L)
#' Q <- induced_distribution(fam, point_mixing(c(0.5, 0.5)))
#' engine <- resolve_engine(exact_engine(), Q, fam)
#' log_expect_q(engine, log(engine@nodes[, 1] + 1))
#' @export
log_expect_q <- function(engine, log_v) {
  logsumexp_weighted(log_v, engine@log_w)
}


#' Standard error of an expectation
#'
#' The Monte Carlo standard error of [expect_q()], and exactly `0` for a
#' deterministic rule. A deterministic rule with approximation error has bias
#' rather than variance, so `0` here is not a claim of exactness.
#' @param engine A [quadrature].
#' @param v Length-`M` numeric vector of integrand values at the nodes.
#' @return Numeric scalar.
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 3L, k = 2L)
#' Q <- induced_distribution(fam, point_mixing(c(0.5, 0.5)))
#' engine <- resolve_engine(mc_engine(200L), Q, fam)
#' expect_se(engine, engine@nodes[, 1])
#' @export
expect_se <- function(engine, v) {
  if (deterministic(engine)) {
    return(0)
  }
  stats::sd(v) / sqrt(n_nodes(engine))
}


# --- Specs --------------------------------------------------------------------
# A spec is a function of (alternative, family) returning a `quadrature`, tagged
# so `resolve_engine()` can tell one from an arbitrary function. Specs are
# resolved late because the nodes depend on Q and the family, which the caller
# supplies at fit time; this also means re-resolving a spec draws a fresh
# sample, which is what certification needs.

new_engine_spec <- function(fn) {
  structure(fn, class = "ripr_engine_spec")
}


#' Exact enumeration over a finite sample space
#'
#' Nodes are the family's full sample space and the weights are the alternative's
#' own probabilities, so expectations are exact. Available only for families
#' with an enumerable sample space.
#'
#' Nodes carrying no mass under the alternative are screened out at resolution.
#' They contribute nothing, and keeping them would put `-Inf` in `log_q`, where
#' the `0 * -Inf` in the objective becomes `NaN`.
#' @return An engine spec for [resolve_engine()].
#' @examples
#' fam <- multinomial_family(n_trials = 3L, k = 2L)
#' Q <- induced_distribution(fam, point_mixing(c(0.5, 0.5)))
#' resolve_engine(exact_engine(), Q, fam)
#' @export
exact_engine <- function() {
  new_engine_spec(function(alternative, family) {
    nodes <- enumerate_space(family@sample_space)
    log_q <- log_density(alternative, nodes)
    live <- is.finite(log_q)
    quadrature(
      nodes = nodes[live, , drop = FALSE],
      log_w = log_q[live],
      log_q = log_q[live],
      family = family,
      deterministic = TRUE
    )
  })
}


#' Monte Carlo integration against draws from the alternative
#'
#' Nodes are `n_draws` draws from the alternative with uniform weights. The
#' draws are frozen when the spec is resolved, so a fit is reproducible under
#' `set.seed()`; resolving the same spec again draws afresh, which is how
#' certification obtains a sample independent of the fit.
#' @param n_draws Number of draws.
#' @return An engine spec for [resolve_engine()].
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 3L, k = 2L)
#' Q <- induced_distribution(fam, point_mixing(c(0.5, 0.5)))
#' resolve_engine(mc_engine(200L), Q, fam)
#' @export
mc_engine <- function(n_draws) {
  n_draws <- as.integer(n_draws)
  if (length(n_draws) != 1L || is.na(n_draws) || n_draws <= 0L) {
    stop("`n_draws` must be a single positive integer.", call. = FALSE)
  }
  new_engine_spec(function(alternative, family) {
    nodes <- draw(alternative, n_draws)
    quadrature(
      nodes = nodes,
      log_w = rep(-log(n_draws), n_draws),
      log_q = log_density(alternative, nodes),
      family = family,
      deterministic = FALSE
    )
  })
}


#' Resolve an engine spec against an alternative and a family
#'
#' Checks that the resulting rule integrates to one. That identity is what makes
#' \eqn{\sum_j w_j G(\theta_j) = 1}{sum_j w_j G(theta_j) = 1} hold for any
#' mixture, and hence what forces the duality gap to be non-negative, so a
#' violation is an error rather than a warning: every downstream quantity would
#' be wrong.
#'
#' @param spec An engine spec, e.g. from [exact_engine()] or [mc_engine()].
#' @param alternative The alternative \eqn{Q}{Q}, an [distribution].
#' @param family A [parametric_family].
#' @param tol Tolerance on the weight sum.
#' @return A [quadrature].
#' @examples
#' fam <- multinomial_family(n_trials = 3L, k = 2L)
#' Q <- induced_distribution(fam, point_mixing(c(0.5, 0.5)))
#' resolve_engine(exact_engine(), Q, fam)
#' @export
resolve_engine <- function(spec, alternative, family, tol = 1e-8) {
  if (!inherits(spec, "ripr_engine_spec")) {
    stop(
      "`spec` must be an engine spec, e.g. `exact_engine()` or `mc_engine(n)`.",
      call. = FALSE
    )
  }
  engine <- spec(alternative, family)
  total <- sum(exp(engine@log_w))
  if (!is.finite(total) || abs(total - 1) > tol) {
    stop(
      "quadrature weights must sum to 1, but sum to ",
      format(total, digits = 8),
      ". Every expectation, and the duality gap with them, depends on this.",
      call. = FALSE
    )
  }
  engine
}


# --- Gauss-Hermite ------------------------------------------------------------

#' Mean and covariance of a distribution, when it is Gaussian
#'
#' The interface [gh_engine()] needs: Gauss-Hermite integrates against a
#' Gaussian weight, so it applies only where the alternative is one. Returns
#' `NULL` for everything else, which is what makes the engine refuse rather than
#' silently integrate against the wrong measure.
#' @param dist An [distribution].
#' @return `list(mean = , cov = )`, or `NULL`.
#' @keywords internal
gaussian_moments <- new_generic("gaussian_moments", "dist", function(dist) {
  S7::S7_dispatch()
})

method(gaussian_moments, distribution) <- function(dist) NULL

method(gaussian_moments, induced_distribution) <- function(dist) {
  induced_gaussian_moments(dist@mixing, dist@family)
}

#' Mean and covariance of the induced mixture, when that mixture is Gaussian
#'
#' Note that these are the moments of \eqn{P_W}{P_W}, not of \eqn{W}{W}: for a
#' Gaussian prior over a Gaussian mean the covariance is
#' \eqn{\Sigma + V}{sigma + V}, where the mixing measure's own is \eqn{V}{V}.
#'
#' Deliberately not named for moments alone. Every mixture has moments, but only
#' a Gaussian one can be integrated by [gh_engine()], and a method supplied for
#' a non-Gaussian mixture would have the engine integrate against a weight that
#' does not match the measure. `NULL` is the correct answer for anything that is
#' not Gaussian, however well defined its moments are.
#'
#' @param mixing A [mixing_measure].
#' @param family A [parametric_family].
#' @return `list(mean = , cov = )`, or `NULL`.
#' @keywords internal
induced_gaussian_moments <- new_generic(
  "induced_gaussian_moments",
  c("mixing", "family"),
  function(mixing, family) S7::S7_dispatch()
)

method(
  induced_gaussian_moments,
  list(mixing_measure, parametric_family)
) <- function(mixing, family) {
  NULL
}


#' Gauss-Hermite nodes and weights for the weight function `exp(-t^2)`
#'
#' Golub-Welsch: the nodes are the eigenvalues of the symmetric tridiagonal
#' Jacobi matrix of the Hermite three-term recurrence, and the weights come from
#' the first component of each orthonormalised eigenvector. Weights sum to
#' `sqrt(pi)`.
#'
#' Golub, G. H. and Welsch, J. H. (1969). Calculation of Gauss quadrature rules.
#' Mathematics of Computation 23(106), 221-230.
#' \doi{10.1090/S0025-5718-69-99647-1}
#' @param n Number of nodes.
#' @return `list(nodes = , weights = )`.
#' @keywords internal
#' @noRd
gauss_hermite <- function(n) {
  if (n == 1L) {
    return(list(nodes = 0, weights = sqrt(pi)))
  }
  i <- seq_len(n - 1L)
  off <- sqrt(i / 2)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(i, i + 1L)] <- off
  jacobi[cbind(i + 1L, i)] <- off
  ev <- eigen(jacobi, symmetric = TRUE)
  ord <- order(ev$values)
  list(
    nodes = ev$values[ord],
    weights = sqrt(pi) * ev$vectors[1L, ord]^2
  )
}


#' Gauss-Hermite quadrature against a Gaussian alternative
#'
#' A deterministic rule for continuous families. Substituting
#' \eqn{x = m + \sqrt{2} L t}{x = m + sqrt(2) L t} turns
#' \eqn{E_{N(m, V)}[f]}{E_N(m,V)[f]} into a Gauss-Hermite integral, so the nodes
#' are an affine image of the tensor-product grid and the weights are the
#' products of the univariate weights, normalised to sum to one.
#'
#' Applies only when the alternative is Gaussian, which [gaussian_moments()]
#' decides. The grid is a tensor product, so it holds `n_nodes^d` points and is
#' unusable beyond a handful of dimensions; `max_nodes` makes that a clear error
#' rather than an exhausted session. An `n`-point rule integrates polynomials of
#' degree at most `2n - 1` exactly, so a smooth integrand needs far fewer nodes
#' than Monte Carlo needs draws.
#'
#' Unlike [mc_engine()] the error here is bias, not variance, so [expect_se()]
#' reports `0` and says nothing about accuracy. Raise `n_nodes` and compare.
#'
#' @param n_nodes Univariate nodes per dimension.
#' @param max_nodes Refuse grids larger than this.
#' @return An engine spec for [resolve_engine()].
#' @references
#'   \insertRef{GolubWelsch1969}{ripr}
#' @examples
#' fam <- gaussian_family(dim = 2L)
#' Q <- induced_distribution(fam, point_mixing(c(0, 0)))
#' resolve_engine(gh_engine(n_nodes = 10L), Q, fam)
#' @export
gh_engine <- function(n_nodes, max_nodes = 1e6) {
  n_nodes <- as.integer(n_nodes)
  if (length(n_nodes) != 1L || is.na(n_nodes) || n_nodes <= 0L) {
    stop("`n_nodes` must be a single positive integer.", call. = FALSE)
  }
  new_engine_spec(function(alternative, family) {
    mom <- gaussian_moments(alternative)
    if (is.null(mom)) {
      stop(
        "Gauss-Hermite quadrature needs a Gaussian alternative; this one is ",
        "not. Use `mc_engine()` instead.",
        call. = FALSE
      )
    }
    d <- length(mom$mean)
    total <- n_nodes^d
    if (total > max_nodes) {
      stop(
        "a ",
        n_nodes,
        "-point grid in ",
        d,
        " dimensions needs ",
        total,
        " nodes, above `max_nodes` (",
        max_nodes,
        "). Lower `n_nodes` or use ",
        "`mc_engine()`.",
        call. = FALSE
      )
    }

    gh <- gauss_hermite(n_nodes)
    grid <- as.matrix(expand.grid(rep(list(seq_len(n_nodes)), d)))
    t_mat <- matrix(gh$nodes[grid], nrow = nrow(grid), ncol = d)
    # Weights are the product of the univariate weights; dividing by pi^(d/2)
    # normalises them to sum to one, which resolve_engine() then verifies.
    log_w <- rowSums(matrix(log(gh$weights[grid]), nrow = nrow(grid))) -
      0.5 * d * log(pi)

    nodes <- t(t(chol(mom$cov)) %*% (sqrt(2) * t(t_mat)) + mom$mean)
    quadrature(
      nodes = nodes,
      log_w = log_w,
      log_q = log_density(alternative, nodes),
      family = family,
      deterministic = TRUE
    )
  })
}
