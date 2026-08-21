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
#' enclose, and the `bound_fn` that does it. A combination missing from this table
#' simply has no implementation yet.
#'
#' The Bernstein enclosure bound for multinomial random variables may be
#' extendable to other families for which the expectation takes the form of a
#' polynomial, and for other bounded geometries (e.g. boxes).
#'
#' The multinomial is the easy case, since its basis *is* the Bernstein basis;
#' the multivariate hypergeometric and multivariate Bernoulli look workable on
#' the same lines and are not done.
#'
#' # The `bound_fn` contract
#'
#' `bound_fn(x, family, subnulls, control)` returns one result per element of
#' `subnulls`, in the same order, each a list with the fields
#' `check_bound_result()` requires.
#'
#' Named `bound_fn` rather than `engine` because `ripr_engine` already means a
#' quadrature rule for expectations under the alternative ([gh_engine()],
#' [mc_engine()]). That is a different object doing a different job, and the two
#' meet in the same conversations often enough for one word to serve both.
#'
#' A `bound_fn` receives a *group* of subnulls that work for a given `(x, family)`,
#' so that enumerating the sample space, evaluating `x` on it, building the
#' lattice (for multinomial), happens just once for all matching subnulls.
#' `certify()` groups the subnulls by resolved method, so a `bound_fn` only
#' ever sees geometries it facilitates, and a null with subnulls of different
#' geometries is split across bounding methods rather than refused.
#'
#' Currently only multinomial is supported, but this architecture makes it
#' easier to extend to other families and geometries later.
#' @return A list of methods, each with `name`, `family`, `subnull`, `bound_fn`
#'   and a one-line `description`.
#' @keywords internal
#' @noRd
certify_methods <- function() {
  list(
    list(
      name = "bernstein",
      family = multinomial_family,
      subnull = simplex_null,
      bound_fn = bernstein_bound,
      description = paste(
        "Bernstein enclosure for multinomial expectations over simplices"
      )
    )
  )
}


#' Check that a bounding method returned what [certify()] needs
#'
#' This enforces a contract between `certify()` and the bounding methods.
#' @keywords internal
#' @noRd
check_bound_result <- function(results, method_name, n_subnulls) {
  if (!is.list(results) || length(results) != n_subnulls) {
    stop(
      "The `",
      method_name,
      "` bounding method returned ",
      length(results),
      " results for ",
      n_subnulls,
      " subnulls.",
      call. = FALSE
    )
  }
  numbers <- c("bound", "incumbent")
  flags <- c("converged", "budget_hit")
  for (r in results) {
    for (field in numbers) {
      if (!is.numeric(r[[field]]) || length(r[[field]]) != 1L) {
        stop(
          "The `",
          method_name,
          "` bounding method returned no scalar `",
          field,
          "`. This is a version mismatch rather than a numerical problem: ",
          "`certify()` and the bounding method disagree about what a ",
          "result looks like.",
          call. = FALSE
        )
      }
    }
    for (field in flags) {
      value <- r[[field]]
      if (!is.logical(value) || length(value) != 1L || is.na(value)) {
        stop(
          "The `",
          method_name,
          "` bounding method returned no `",
          field,
          "`. This is a version mismatch rather than a numerical problem: ",
          "`certify()` and the bounding method disagree about what a ",
          "result looks like.",
          call. = FALSE
        )
      }
    }
    if (!isTRUE(r$converged) && !isTRUE(r$budget_hit)) {
      stop(
        "The `",
        method_name,
        "` bounding method reported a search that stopped ",
        "without recording why.",
        call. = FALSE
      )
    }
    if (!is.integer(r$iterations) || length(r$iterations) != 1L) {
      stop(
        "The `",
        method_name,
        "` bounding method returned no integer `iterations`.",
        call. = FALSE
      )
    }
  }
  invisible(results)
}


#' Bernstein enclosure over simplices, for multinomial expectations
#'
#' The pmf of a multinomial *is* the degree-`n` Bernstein basis, so `x`
#' evaluated on the sample space is already a coefficient vector for
#' \eqn{E_\theta[X]}{E_theta[X]} and no conversion from a power form is needed.
#' Each subnull is reparametrised onto its own simplex and enclosed separately.
#' @keywords internal
#' @noRd
bernstein_bound <- function(x, family, subnulls, control) {
  check_bernstein_size(family@n_trials, family@k, control$max_coefficients)

  outcomes <- enumerate_space(family@sample_space)
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

  lattice <- bernstein_lattice(family@n_trials, family@k)
  lapply(subnulls, function(s) {
    box <- list(
      V = s@vertices,
      coef = reparametrise_to(values, lattice, s@vertices)
    )
    certify_sup(
      list(box),
      lattice,
      tol = control$tol,
      max_iter = control$max_nodes
    )
  })
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


#' The expectation of a random variable as a function of the parameter
#'
#' \eqn{E_\theta[X] = \sum_x P_\theta(x) X(x)}{E_theta[X] = sum_x P_theta(x) X(x)},
#' with `X` evaluated once at the whole sample space and reused.
#' @keywords internal
#' @noRd
expectation_objective <- function(family, values) {
  loglik <- compile_loglik(family, enumerate_space(family@sample_space))
  objective(
    value = function(theta) {
      sum(exp(as.vector(loglik(matrix(theta, ncol = 1L)))) * values)
    },
    grad = function(theta) {
      weight <- exp(as.vector(loglik(matrix(theta, ncol = 1L)))) * values
      as.vector(crossprod(score(family, theta, enumerate_space(family@sample_space)), weight))
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
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_null(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_null(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' X <- mixture_likelihood(mixture(point_mixing(c(0.4, 0.35, 0.25)), fam))
#' sup_lb(X, plurality)
#' @export
sup_lb <- function(x, null, n_seeds = 200L, n_restarts = 25L) {
  if (!S7_inherits(x, random_variable)) {
    stop("`x` must be a `random_variable`.", call. = FALSE)
  }

  family <- null@family
  values <- x(enumerate_space(family@sample_space))
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


#' Flatten a run's branch-and-bound nodes into a table
#'
#' One row per node ever created, with the iteration it appeared (`born`), the
#' iteration it left the active set (`retired`) and why (`fate`: `"split"`,
#' `"pruned"`, or `"active"` for one still live at the end).
#'
#' Coefficients are dropped. They are the bulk of a node (10,626 doubles each at
#' `K = 5, n = 20`, against a handful for everything else here) and nothing
#' downstream of a finished run evaluates them.
#' @keywords internal
#' @noRd
node_table <- function(result, subnull) {
  nodes <- result$history
  if (!length(nodes)) {
    return(NULL)
  }
  field <- function(name, template) {
    vapply(nodes, function(b) b[[name]] %||% template, template)
  }
  data.frame(
    subnull = as.integer(subnull),
    id = field("id", NA_integer_),
    parent = field("parent", NA_integer_),
    depth = field("depth", NA_integer_),
    born = field("born", NA_integer_),
    retired = field("retired", NA_integer_),
    fate = field("fate", NA_character_),
    upper = field("ub", NA_real_),
    volume = vapply(nodes, function(b) abs(det(b$V)), numeric(1L)),
    vertices = I(lapply(nodes, function(b) b$V)),
    row.names = NULL
  )
}


#' Record how a certification ran, for inspection and plotting
#'
#' Same computation as [certify()], reporting the branch-and-bound tree instead
#' of the certificate. One row per leaf of the search, per subnull.
#'
#' The nodes live at any iteration tile their subnull exactly, so at `K = 3` the
#' `vertices` column draws the partition of the facet directly, in barycentric
#' coordinates, at every step and not merely at the end.
#' `depth` against `born` shows the shape of the tree: a max-bound queue rule
#' can produce a chain rather than anything balanced, which is visible here and
#' nowhere else.
#'
#' Only bounding methods that use branch and bound populate this. A method that
#' bounds in closed form contributes no rows.
#' @inheritParams certify
#' @return A data frame with `subnull`, `id`, `parent`, `depth`, `born`,
#'   `retired`, `fate`, `upper`, `volume` and a `vertices` list column, plus the
#'   certificate itself in the `"certificate"` attribute and the per-iteration
#'   bound in `"trace"`.
#' @seealso [certify()]
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_null(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_null(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' X <- mixture_likelihood(mixture(point_mixing(c(0.4, 0.35, 0.25)), fam))
#' nodes <- certify_trace(X, plurality, tol = 1e-6)
#' nrow(nodes)
#' @export
certify_trace <- function(
  x,
  null,
  tol = 1e-6,
  max_nodes = 20000L,
  max_coefficients = 1024^2
) {
  result <- certify(
    x,
    null,
    tol = tol,
    max_nodes = max_nodes,
    max_coefficients = max_coefficients,
    .record = TRUE
  )
  nodes <- result$record
  traces <- result$traces
  incumbent_traces <- result$incumbent_traces
  result$record <- NULL
  result$traces <- NULL
  result$incumbent_traces <- NULL
  attr(nodes, "trace") <- traces
  attr(nodes, "incumbent_trace") <- incumbent_traces
  attr(nodes, "certificate") <- result
  nodes
}


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
#' @param .record Also return branch-and-bound nodes. Use [certify_trace()]
#'   rather than this directly.
#' @return A list with `sup_ub`, `sup_lb`, the `random_variable` and `null` it
#'   holds for, the `method` names that produced it (one per distinct subnull
#'   geometry), and per-subnull `bounds`, `incumbents`, `iterations`,
#'   `converged` and `budget_hit`. `budget_hit` flags when a search stopped at
#'   `max_nodes` with the gap still open, so its bound is valid but likely
#'   loose.
#' @seealso [sup_lb()]
#' @references
#' \insertAllCited{}
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_null(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_null(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' X <- mixture_likelihood(mixture(point_mixing(c(0.4, 0.35, 0.25)), fam))
#' cert <- certify(X, plurality, tol = 1e-6)
#' c(upper = cert$sup_ub, attained = cert$sup_lb)
#' @export
certify <- function(
  x,
  null,
  tol = 1e-6,
  max_nodes = 20000L,
  max_coefficients = 1024^2,
  .record = FALSE
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
    unimpl_msgs <- vapply(
      null@subnulls[unavailable],
      function(s) unimplemented_message(family, s),
      character(1L)
    )
    stop(
      "Cannot certify:\n",
      paste(unique(unimpl_msgs), collapse = "\n"),
      call. = FALSE
    )
  }
  method_names <- unique(vapply(methods, function(m) m$name, character(1L)))

  # Group the subnulls by resolved method, run each bound_fn once on its group,
  # then put the results back in subnull order. With one entry in the registry
  # this is a single group; the grouping is what stops that being an assumption.
  control <- list(
    tol = tol,
    max_nodes = max_nodes,
    max_coefficients = max_coefficients
  )
  per_subnull <- vector("list", length(null@subnulls))
  for (name in method_names) {
    which_subnulls <- which(
      vapply(methods, function(m) m$name, character(1L)) == name
    )
    method <- methods[[which_subnulls[[1L]]]]
    results <- method$bound_fn(
      x,
      family,
      null@subnulls[which_subnulls],
      control
    )
    check_bound_result(results, name, length(which_subnulls))
    per_subnull[which_subnulls] <- results
  }

  bounds <- vapply(per_subnull, function(r) r$bound, numeric(1L))
  incumbents <- vapply(per_subnull, function(r) r$incumbent, numeric(1L))
  iterations <- vapply(per_subnull, function(r) r$iterations, integer(1L))
  converged <- vapply(per_subnull, function(r) r$converged, logical(1L))
  budget_hit <- vapply(per_subnull, function(r) r$budget_hit, logical(1L))
  out <- list(
    sup_ub = max(bounds),
    sup_lb = max(incumbents),
    random_variable = x,
    null = null,
    method = method_names,
    bounds = bounds,
    incumbents = incumbents,
    iterations = iterations,
    converged = converged,
    budget_hit = budget_hit
  )
  if (.record) {
    # Records the tables of the node histories for each subnull. Used by
    # `certify_trace`.
    tables <- lapply(
      seq_along(per_subnull),
      function(i) node_table(per_subnull[[i]], i)
    )
    out$record <- do.call(rbind, Filter(Negate(is.null), tables))
    out$traces <- lapply(per_subnull, function(r) r$trace)
    out$incumbent_traces <- lapply(per_subnull, function(r) r$incumbent_trace)
  }
  out
}
