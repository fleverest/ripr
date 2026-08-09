#' @include random_variable.R null.R bernstein.R
NULL

# Bounds on the expectation of a random variable under a particular null.
#
# Two functions, and the difference between them is the whole point.
# `sup_lb()` searches and reports what it found: a *lower* bound, useful
# for diagnosis and useless as a guarantee. `certify()` encloses and reports
# what it proved: an *upper* bound, and the only thing an e-variable may rest
# on.

# --- Size guard for Bernstein -------------------------------------------------

#' How many Bernstein coefficients a certification would need to bound a
#' multinomial expectation.
#'
#' For a multinomial random variable `X`, `E_theta[X]` takes the form of a
#' polynomial in `theta` of degree `n_trials`. Its Bernstein form has one
#' coefficient per point of the multinomial sample space.
#' @keywords internal
#' @noRd
bernstein_size <- function(n_trials, k) choose(n_trials + k - 1, k - 1)


#' The largest batch (`n_trials`) that would fit within a coefficient budget
#'
#' Reported in the refusal, since the batch size is one thing a caller can
#' actually change.
#' @keywords internal
#' @noRd
largest_batch <- function(k, max_coefficients) {
  # Binary search over n_trials over the support size
  if (bernstein_size(1L, k) > max_coefficients) {
    return(0L)
  }
  hi <- 1L
  while (bernstein_size(2L * hi, k) <= max_coefficients) {
    hi <- 2L * hi
  }
  lo <- hi
  hi <- 2L * hi
  while (lo < hi) {
    mid <- (lo + hi + 1L) %/% 2L
    if (bernstein_size(mid, k) <= max_coefficients) {
      lo <- mid
    } else {
      hi <- mid - 1L
    }
  }
  lo
}


#' Refuse a certification too large to attempt
#'
#' Checked before anything is built, since the lattice is where the memory and
#' most of the setup time go. Exact arithmetic on `choose()`, so the guard itself
#' costs nothing.
#'
#' This is a resource limit rather than a correctness one: raising it costs time
#' and memory and nothing else. The other guards, e.g. checking that the family
#' permits a bound on the expectation at all, are a correctness concern and can
#' not be overridden.
#' @keywords internal
#' @noRd
check_bernstein_size <- function(n_trials, k, max_coefficients) {
  size <- bernstein_size(n_trials, k)
  if (size <= max_coefficients) {
    return(invisible(size))
  }
  stop(
    "certifying would need ",
    format(size, big.mark = ",", scientific = FALSE),
    " Bernstein coefficients for n_trials = ",
    n_trials,
    " with k = ",
    k,
    ", above `max_coefficients` (",
    format(max_coefficients, big.mark = ",", scientific = FALSE),
    ").\n",
    "The count is choose(n_trials + k - 1, k - 1), so it is the batch size ",
    "that drives it: ",
    largest_batch(k, max_coefficients),
    " would fit.\n",
    "Raise `max_coefficients` to attempt it anyway.",
    call. = FALSE
  )
}


# --- What is implemented ------------------------------------------------------

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
certify_method <- function(family, null) {
  for (method in certify_methods()) {
    fits <- S7_inherits(family, method$family) &&
      all(vapply(
        null@subnulls,
        function(s) S7_inherits(s, method$subnull),
        logical(1)
      ))
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
unimplemented_message <- function(family, null) {
  geometries <- unique(vapply(null@subnulls, class_name, character(1)))
  available <- vapply(
    certify_methods(),
    function(m) paste0("  - ", m$description, "."),
    character(1)
  )
  paste0(
    "no bounding method is implemented for ",
    class_name(family),
    " expectations over ",
    paste(geometries, collapse = " / "),
    ".\n",
    "Implemented:\n",
    paste(available, collapse = "\n"),
    "\n",
    "Certifying this would mean deriving and implementing a bound on ",
    class_name(family),
    " expectations over ",
    paste(geometries, collapse = " / "),
    ". Nothing here says one does not ",
    "exist. In the meantime `sup_lb()` still searches, and reports a ",
    "lower bound."
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
#' `H0`: For any non-negative `X` and any scalar `b` at least this supremum,
#' `X / B` has expectation at most 1 under `H0`. So `certify(X, H0)` returning
#' `1` or less says `X` is already an e-variable for `H0`, and otherwise
#' `X / bound` is one.
#'
#' The method refuses to certify where a family with no known bounding method.
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
#' @param max_nodes Cap on subdivisions per subnull for branch-and-bound.
#' @param max_coefficients Refuse above this many Bernstein coefficients.
#' @return A list with `sup_bound`, the `random_variable` and `null` it holds
#'   for, the `method` that produced it, and per-subnull `bounds`, `incumbents`,
#'   `nodes` and whether each search `exhausted` its budget.
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
  method <- certify_method(family, null)
  if (is.null(method)) {
    stop("cannot certify: ", unimplemented_message(family, null), call. = FALSE)
  }
  check_bernstein_size(family@n_trials, family@k, max_coefficients)

  outcomes <- support(family)
  values <- x(outcomes)
  if (any(!is.finite(values))) {
    # Not a gap in what is implemented: the supremum really is infinite, and no
    # method could return anything else. A likelihood ratio does this wherever
    # its denominator vanishes and its numerator does not.
    stop(
      "cannot certify: the variable is not finite everywhere on the sample ",
      "space, so its null expectation is unbounded.",
      call. = FALSE
    )
  }

  # The enumeration must be the family's own, or every coefficient is
  # mis-indexed; passing `outcomes` makes the lattice verify that rather than
  # assume it.
  lattice <- bernstein_lattice(family@n_trials, family@k, tally = outcomes)

  per_subnull <- lapply(null@subnulls, function(s) {
    box <- list(
      V = s@vertices,
      coef = reparametrise_to(values, lattice, s@vertices)
    )
    certify_sup(list(box), lattice, tol = tol, max_iter = max_nodes)
  })

  bounds <- vapply(per_subnull, function(r) r$bound, numeric(1))
  list(
    sup_bound = max(bounds),
    random_variable = x,
    null = null,
    method = method$name,
    bounds = bounds,
    incumbents = vapply(per_subnull, function(r) r$incumbent, numeric(1)),
    nodes = vapply(per_subnull, function(r) r$iterations, numeric(1)),
    exhausted = vapply(per_subnull, function(r) isTRUE(r$exhausted), logical(1))
  )
}
