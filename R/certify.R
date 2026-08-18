#' @include random_variable.R null.R bernstein.R
NULL

# Bounds on the expectation of a random variable under a particular null.
#
# Two functions, and the difference between them is the whole point.
# `sup_lb()` searches and reports what it found: a *lower* bound, useful
# for diagnosis and useless as a guarantee. `certify()` encloses and reports
# what it proved: an *upper* bound, and the only thing an e-variable may rest
# on.

#' Bounding methods available to [certify()]
#'
#' Each entry names a family of expectations and a subnull geometry it can
#' enclose. A combination missing from this table simply has no implementation
#' yet.
#'
#' The Bernstein enclosure bound for multinomial random variables may be
#' extendable to other families for which the expectation takes the form of a
#' polynomial, and for other bounded geometries (e.g. boxes).
#'
#' The multinomial is the easy case, since its basis *is* the Bernstein basis;
#' the multivariate hypergeometric and multivariate Bernoulli look workable on
#' the same lines and are not done.
#' @return A list of methods, each with `name`, `family`, `subnull` and a
#'   one-line `description`.
#' @keywords internal
#' @noRd
certify_methods <- function() {
  list(
    list(
      name = "bernstein",
      family = multinomial_family,
      subnull = simplex_null,
      description = paste(
        "Bernstein enclosure for multinomial expectations over simplices"
      )
    )
  )
}


#' Fetch a method that can certify a (family, subnull) combination, or `NULL`
#' @keywords internal
#' @noRd
certify_method <- function(family, subnull) {
  for (method in certify_methods()) {
    fits <- S7_inherits(family, method$family) &&
      S7_inherits(subnull, method$subnull)
    if (fits) {
      return(method)
    }
  }
  NULL
}


#' Name a class as it should appear in a message
#' @keywords internal
#' @noRd
class_name <- function(x) attr(S7_class(x), "name")


#' Explain that no implemented method covers this (family, subnull) combination
#'
#' Says what is missing rather than what is impossible. Certifying some other
#' pairing of family and geometry is a matter of deriving and implementing a
#' bound, not of the thing being unbounded.
#' @keywords internal
#' @noRd
unimplemented_message <- function(family, subnull) {
  geometry <- class_name(subnull)
  paste0(
    "No bounding method is implemented for ",
    class_name(family),
    " expectations over ",
    geometry,
    ".\n",
    "Certifying this requires deriving and implementing a bound on ",
    class_name(family),
    " expectations over ",
    geometry,
    ". Nothing here says one does not exist. In the meantime, `sup_lb()` ",
    "still searches, and reports a lower bound."
  )
}


# --- Estimating ---------------------------------------------------------------

#' The expectation of a random variable as a function of the parameter
#'
#' \eqn{E_\theta[X] = \sum_x P_\theta(x) X(x)}{E_theta[X] = sum_x P_theta(x) X(x)},
#' with `X` evaluated once at the whole sample space and reused.
#' @keywords internal
#' @noRd
expectation_objective <- function(family, values) {
  loglik <- compile_loglik(family, support(family))
  objective(
    value = function(theta) {
      sum(exp(as.vector(loglik(matrix(theta, ncol = 1L)))) * values)
    },
    grad = function(theta) {
      weight <- exp(as.vector(loglik(matrix(theta, ncol = 1L)))) * values
      as.vector(crossprod(score(family, theta, support(family)), weight))
    },
    value_batch = function(theta_mat) {
      as.vector(crossprod(exp(loglik(theta_mat)), values))
    }
  )
}


#' Estimate the largest null expectation of a random variable
#'
#' Multi-start local ascent, a **lower** bound on the supremum: it reports the
#' largest value it managed to find through optimisation, though a larger one
#' value may exist somewhere it did not look. Use it as a diagnoses rather than
#' treating it like a certificate on the bound. See [certify()] for proven
#' upper bounds where supported.
#'
#' Cheap by comparison, and defined wherever the search is. An unbounded subnull
#' or a continuous parameter space can still be searched, though the lower bound
#' may not be very tight.
#' @param x A [random_variable].
#' @param null A [null_model].
#' @param n_seeds,n_restarts Resolution of the search.
#' @return A list with `sup_lb` and the `theta` attaining it.
#' @seealso [certify()]
#' @export
sup_lb <- function(x, null, n_seeds = 200L, n_restarts = 25L) {
  if (!S7_inherits(x, random_variable)) {
    stop("`x` must be a `random_variable`.", call. = FALSE)
  }

  family <- null@family
  values <- x(support(family))
  obj <- expectation_objective(family, values)

  found <- lapply(
    null@subnulls,
    function(s) {
      maximise_over(s, obj, n_seeds = n_seeds, n_restarts = n_restarts)
    }
  )
  best <- which.max(vapply(found, function(f) f$value, numeric(1)))
  list(
    sup_lb = found[[best]]$value,
    theta = found[[best]]$theta,
    subnull = best
  )
}


# --- Certifying ---------------------------------------------------------------

#' Certify an upper bound on the largest null expectation
#'
#' A **proven** upper bound on
#' \eqn{\sup_{\theta \in \Theta_0} E_\theta[X]}{sup_theta E_theta[X]}, where
#' available.
#'
#' This exists to prove that a particular random variable is an e-variable for
#' `H0`: for any non-negative `X` and any scalar `b` at least this supremum,
#' `X / b` has expectation at most 1 under `H0`. So `certify(X, H0)` returning
#' `1` or less says `X` is already an e-variable for `H0`, and otherwise
#' `X / bound` is one.
#'
#' Certification is refused for any (family, subnull) combination with no known
#' bounding method.
#'
#' Currently only multinomial families with simplex subnulls are supported. The
#' upper bound is found via a branch-and-bound algorithm that recursively
#' subdivides the simplex \insertCite{Leroy2012}{ripr}, and bounds each subset
#' via the simplicial Bernstein range enclosure property
#' \insertCite{Garloff1986}{ripr}.
#'
#' Unlike [sup_lb()] this is a proven bound.
#' @param x A [random_variable].
#' @param null A [null_model].
#' @param tol Stop once the bound is within `tol` of the best value found.
#' @param max_nodes Cap on subdivisions per subnull for branch-and-bound
#'   algorithms.
#' @param max_coefficients Refuse above this many Bernstein coefficients
#'   (for bounding multinomial expectation in simplices).
#' @return A list with `sup_ub`, `sup_lb`, the `random_variable` and `null` it
#'   holds for, the `method` names that produced it (one per distinct subnull
#'   geometry), and per-subnull `bounds`, `incumbents`, `nodes` and
#'   `exhausted`.
#' @seealso [sup_lb()]
#' @references
#' \insertAllCited{}
#' @export
certify <- function(
  x,
  null,
  tol = 1e-6,
  max_nodes = 20000L,
  max_coefficients = 1024^2
) {
  if (!S7_inherits(x, random_variable)) {
    stop("`x` must be a `random_variable`.", call. = FALSE)
  }
  rlang::check_number_decimal(tol, min = 0)
  rlang::check_number_whole(max_nodes, min = 1, max = 2147483647)
  rlang::check_number_whole(max_coefficients, min = 1, max = 2147483647)

  family <- null@family
  # Get applicable bound method for each subnull.
  methods <- lapply(
    null@subnulls,
    function(s) certify_method(family, s)
  )
  # Stop if any (family, subnull) combination is not implemented. Report the
  # offending subnulls rather than the null model, and deduplicate: a plurality
  # null has one subnull per candidate and they share a geometry, so the same
  # message would otherwise repeat K - 1 times.
  unavailable <- vapply(methods, is.null, logical(1L))

  if (any(unavailable)) {
    unimpl_msgs <- null@subnulls[which(unavailable)] |>
      lapply(function(s) unimplemented_message(family, s))
    stop(
      "Cannot certify:\n",
      paste(unique(unimpl_msgs), collapse = "\n"),
      call. = FALSE
    )
  }
  method_names <- unique(vapply(methods, function(m) m$name, character(1L)))

  check_bernstein_size(family@n_trials, family@k, max_coefficients)

  # TODO: generalise for continuous support...
  # This doesn't matter yet since we have no bound for the Gaussian case. If one
  # is implemented then we would need to restructure this stopping rule.
  outcomes <- support(family)
  values <- x(outcomes)
  if (any(!is.finite(values))) {
    # If the random variable is not bounded, the expectation is not well-defined
    # and therefore no bound is coherent. This could arise, for instance, in
    # mixture likelihood ratios where the numerator and denominator are not
    # absolutely continuous.
    stop(
      "Cannot certify: the variable is not finite everywhere on the sample ",
      "space, so its null expectation is unbounded.",
      call. = FALSE
    )
  }

  # ---- Multinomial case is assumed here as we have only implemented that ----

  # The Bernstein lattice (i.e. multinomial support)
  lattice <- bernstein_lattice(family@n_trials, family@k)

  per_subnull <- lapply(null@subnulls, function(s) {
    box <- list(
      V = s@vertices,
      coef = reparametrise_to(values, lattice, s@vertices)
    )
    certify_sup(list(box), lattice, tol = tol, max_iter = max_nodes)
  })

  bounds <- vapply(per_subnull, function(r) r$bound, numeric(1L))
  incumbents <- vapply(per_subnull, function(r) r$incumbent, numeric(1L))
  nodes <- vapply(per_subnull, function(r) r$iterations, integer(1L))
  exhausted <- vapply(per_subnull, function(r) isTRUE(r$exhausted), logical(1L))

  list(
    sup_ub = max(bounds),
    sup_lb = max(incumbents),
    random_variable = x,
    null = null,
    method = method_names,
    bounds = bounds,
    incumbents = incumbents,
    nodes = nodes,
    exhausted = exhausted
  )
}
