#' @include family.R distribution.R gauss.R
NULL

#' Quadrature rule for expectations under the alternative
#'
#' Every quantity the optimiser needs (objective, oracle, EM responsibilities,
#' gap) is an expectation under Q; an engine is nodes `(M, K)` with log weights
#' `log_w` (summing to 1) approximating it. `log_q`, the alternative's log
#' density at the nodes, equals `log_w` only for the exact engine. Build with
#' `resolve_engine()`, which checks the weights.
#' @keywords internal
#' @noRd
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


#' @noRd
#' @export
method(print, quadrature) <- function(x, ...) {
  cat("<", class_name(x), ">\n", sep = "")
  cat("  ", format(x@family), "\n", sep = "")
  cat(
    "  ",
    n_nodes(x),
    ngettext(n_nodes(x), " node", " nodes"),
    ", ",
    if (deterministic(x)) "deterministic" else "stochastic",
    "\n",
    sep = ""
  )
  invisible(x)
}


#' @noRd
#' @export
method(format, quadrature) <- function(x, ...) {
  sprintf(
    "%s: %d nodes over %s, %s",
    class_name(x),
    n_nodes(x),
    class_name(x@family),
    if (deterministic(x)) "deterministic" else "stochastic"
  )
}


#' Number of quadrature nodes
#' @keywords internal
#' @noRd
n_nodes <- function(engine) nrow(engine@nodes)


#' Is this rule free of sampling error?
#'
#' `FALSE` only for Monte Carlo. A deterministic rule may still be biased, but
#' no standard error describes that.
#' @keywords internal
#' @noRd
deterministic <- function(engine) isTRUE(engine@deterministic)


#' Compile the family's log-likelihood at the engine's nodes
#'
#' Returns a function of `theta_mat` giving the `(M, C)` log densities. Compile
#' once per step and reuse.
#' @keywords internal
#' @noRd
compile_engine <- function(engine) {
  compile_loglik(engine@family, engine@nodes)
}


#' `E_Q[v]` for length-`M` values `v` at the nodes, which may be negative
#' @keywords internal
#' @noRd
expect_q <- function(engine, v) {
  sum(exp(engine@log_w) * v)
}


#' `log E_Q[v]` from `log_v`, for `v >= 0`
#' @keywords internal
#' @noRd
log_expect_q <- function(engine, log_v) {
  matrixStats::logSumExp(log_v + engine@log_w)
}


# --- Specs --------------------------------------------------------------------

#' A recipe for a quadrature rule, resolved once the alternative is known
#'
#' Re-resolving draws a fresh sample, which certification relies on. `resolve`
#' is a function of `(alternative, family)`.
#' @keywords internal
#' @noRd
engine_spec <- new_class(
  "engine_spec",
  properties = list(
    label = class_character,
    resolve = class_function
  ),
  validator = function(self) {
    if (length(self@label) != 1L || is.na(self@label)) {
      return("`label` must be a single string")
    }
    if (!identical(names(formals(self@resolve)), c("alternative", "family"))) {
      return("`resolve` must be a function of `(alternative, family)`")
    }
    NULL
  }
)


#' @noRd
#' @export
method(format, engine_spec) <- function(x, ...) x@label


#' @noRd
#' @export
method(print, engine_spec) <- function(x, ...) {
  cat("<", class_name(x), "> ", format(x), "\n", sep = "")
  invisible(x)
}


#' Exact enumeration over a finite sample space
#'
#' Nodes are the family's full sample space, weighted by the alternative's
#' probabilities, so expectations are exact. Needs an enumerable sample space.
#' Outcomes with zero probability are dropped.
#' @return An engine spec, to pass as the `engine` argument of [ripr_init()].
#' @examples
#' fam <- multinomial_family(n_trials = 3L, k = 2L)
#' null <- null_model(fam, simplex_region(vertices = rbind(c(0.5, 0.5), c(0, 1))))
#' ripr_init(fam(c(0.6, 0.4)), null, engine = exact_engine())
#' @export
exact_engine <- function() {
  resolve <- function(alternative, family) {
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
  }
  engine_spec(label = "exact_engine()", resolve = resolve)
}


#' Monte Carlo integration against draws from the alternative
#'
#' Nodes are `n` draws from the alternative with equal weights, fixed for the
#' fit (reproducible under `set.seed()`). Certification re-resolves the spec to
#' get an independent sample.
#' @param n Number of draws.
#' @return An engine spec, to pass as the `engine` argument of [ripr_init()].
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 3L, k = 2L)
#' null <- null_model(fam, simplex_region(vertices = rbind(c(0.5, 0.5), c(0, 1))))
#' ripr_init(fam(c(0.6, 0.4)), null, engine = mc_engine(200L))
#' @export
mc_engine <- function(n) {
  rlang::check_number_whole(n, min = 1)
  n <- as.integer(n)
  resolve <- function(alternative, family) {
    nodes <- draw(alternative, n)
    quadrature(
      nodes = nodes,
      log_w = rep(-log(n), n),
      log_q = log_density(alternative, nodes),
      family = family,
      deterministic = FALSE
    )
  }
  engine_spec(label = sprintf("mc_engine(%d)", n), resolve = resolve)
}


#' Resolve an engine spec against an alternative and a family
#'
#' Errors unless the weights sum to one (within `tol`): that identity gives
#' `sum_j w_j G(theta_j) = 1` for any mixture, which keeps the duality gap
#' non-negative. Every spec must satisfy it.
#' @keywords internal
#' @noRd
resolve_engine <- function(spec, alternative, family, tol = 1e-8) {
  if (!S7_inherits(spec, engine_spec)) {
    stop(
      "`spec` must be an engine spec, e.g. `exact_engine()` or `mc_engine(n)`.",
      call. = FALSE
    )
  }
  engine <- spec@resolve(alternative, family)
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

#' Mean and covariance of the induced mixture, when that mixture is Gaussian
#'
#' Used by [gh_engine()]. These are the moments of \eqn{P_W}{P_W}, not
#' \eqn{W}{W} (e.g. \eqn{\Sigma + V}{sigma + V} for a Gaussian prior on a
#' Gaussian mean). `NULL` for any non-Gaussian mixture.
#'
#' @param mixing A [distribution] over the family's parameter space.
#' @param family A [parametric_family].
#' @return `list(mean = , cov = )`, or `NULL`.
#' @keywords internal
mixture_gaussian_moments <- new_generic(
  "mixture_gaussian_moments",
  c("mixing", "family"),
  function(mixing, family) S7::S7_dispatch()
)

method(
  mixture_gaussian_moments,
  list(distribution, parametric_family)
) <- function(mixing, family) {
  NULL
}


#' Gauss-Hermite quadrature against a Gaussian alternative
#'
#' A deterministic rule for continuous families, applicable only when the
#' alternative is Gaussian (as decided by [mixture_gaussian_moments()]). Nodes
#' are an affine image of a tensor-product Gauss-Hermite grid, exact for
#' polynomials of degree up to `2 * n_nodes - 1`.
#'
#' The grid has `n_nodes^d` points, so it suits only a few dimensions. Its error
#' is bias, not variance: raise `n_nodes` and compare.
#'
#' @param n_nodes Univariate nodes per dimension.
#' @param max_nodes Refuse grids larger than this.
#' @return An engine spec, to pass as the `engine` argument of [ripr_init()].
#' @references
#'   \insertRef{GolubWelsch1969}{ripr}
#' @examples
#' fam <- gaussian_family(d = 2L)
#' null <- null_model(fam, halfspace_region(normal = c(1, -1)))
#' ripr_init(fam(c(1, 0)), null, engine = gh_engine(n_nodes = 10L))
#' @export
gh_engine <- function(n_nodes, max_nodes = 1e6) {
  n_nodes <- as.integer(n_nodes)
  if (length(n_nodes) != 1L || is.na(n_nodes) || n_nodes <= 0L) {
    stop("`n_nodes` must be a single positive integer.", call. = FALSE)
  }

  resolve <- function(alternative, family) {
    mom <- if (S7_inherits(alternative, mixture)) {
      mixture_gaussian_moments(alternative@mixing, alternative@family)
    }
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

    grid <- tensor_rule(rep(list(gauss_hermite(n_nodes)), d))
    nodes <- add_by_col(sqrt(2) * grid$nodes %*% chol(mom$cov), mom$mean)
    quadrature(
      nodes = nodes,
      # Normalise the Hermite mass pi^(d/2) to one.
      log_w = grid$log_w - 0.5 * d * log(pi),
      log_q = log_density(alternative, nodes),
      family = family,
      deterministic = TRUE
    )
  }
  engine_spec(label = sprintf("gh_engine(%d)", n_nodes), resolve = resolve)
}
