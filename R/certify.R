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
#' Each entry names a family of expectations and a region geometry it can
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
#' `bound_fn(x, family, parts, control)` returns one result per element of
#' `parts`, in the same order, each a list with the fields
#' `check_bound_result()` requires.
#'
#' Named `bound_fn` rather than `engine` because `ripr_engine` already means a
#' quadrature rule for expectations under the alternative ([gh_engine()],
#' [mc_engine()]). That is a different object doing a different job, and the two
#' meet in the same conversations often enough for one word to serve both.
#'
#' A `bound_fn` receives a *group* of parts that work for a given `(x, family)`,
#' so that enumerating the sample space, evaluating `x` on it, building the
#' lattice (for multinomial), happens just once for all matching parts.
#' `certify()` groups the parts by resolved method, so a `bound_fn` only
#' ever sees geometries it facilitates, and a null with parts of different
#' geometries is split across bounding methods rather than refused.
#'
#' Currently only multinomial is supported, but this architecture makes it
#' easier to extend to other families and geometries later.
#' @return A list of methods, each with `name`, `family`, `accepts`,
#'   `bound_fn` and a one-line `description`. `accepts` is a predicate on a
#'   single convex cell (a [convex_region], not a whole [region]), and a
#'   predicate rather than a class.
#' @keywords internal
#' @noRd
certify_methods <- function() {
  list(
    list(
      name = "bernstein",
      family = multinomial_family,
      accepts = bernstein_compatible,
      bound_fn = bernstein_bound,
      description = paste(
        "Bernstein enclosure for multinomial expectations over simplices"
      )
    )
  )
}


#' Can the Bernstein enclosure bound an expectation over this region?
#'
#' The three conditions `reparametrise_to()` asserts, checked before it is
#' reached. De Casteljau pushes the degree-`n` basis on the standard simplex
#' onto the region's vertices, which needs a square, non-singular,
#' simplex-valued vertex matrix.
#'
#' [simplex_region()] guarantees only the non-singular part, because the other
#' two are not properties of being a simplex: a tetrahedron in `R^3` and a
#' segment inside the 2-simplex are both simplices, and neither can carry this
#' enclosure. A lower-dimensional simplex would need a lower-degree lattice,
#' which is a separate `certify_methods()` entry rather than a loosened
#' predicate here.
#' @param space A [convex_region].
#' @param tol Tolerance for testing rank-deficiency and sum-to-one constraint.
#' @return `TRUE` or `FALSE`.
#' @keywords internal
#' @noRd
bernstein_compatible <- function(space, tol = 1e-9) {
  if (!S7_inherits(space, simplex_region)) {
    return(FALSE)
  }
  v <- space@vertices
  ncol(v) == nrow(v) &&
    all(v >= -1e-12) &&
    max(abs(colSums(v) - 1)) < tol &&
    simplex_rcond(v) > tol
}


#' Reciprocal condition number of a simplex's edge matrix
#'
#' The conditioning heuristic that used to live in `simplex_region`'s
#' validator, now a certification concern: the validator's affine-independence
#' test is exact, so a thin sliver is a genuine simplex, but
#' `reparametrise_to()` inverts the vertex matrix and an ill-conditioned one
#' cannot be enclosed reliably. Measured on the edge matrix rather than a
#' determinant, which may not exist (the vertex matrix need not be square)
#' and would not be scale invariant if it did.
#' @keywords internal
#' @noRd
simplex_rcond <- function(v) {
  edges <- v[, -1L, drop = FALSE] - v[, 1L]
  sv <- svd(edges, nu = 0L, nv = 0L)$d
  if (sv[1L] <= 0) {
    return(0)
  }
  sv[length(sv)] / sv[1L]
}


#' Why the Bernstein enclosure cannot handle this region, or `NULL` if it can
#'
#' Names the condition that fails. Currently only meaningful for a
#' [simplex_region()]; anything else is already refused by class.
#' @keywords internal
#' @noRd
bernstein_obstruction <- function(space) {
  if (bernstein_compatible(space)) {
    return(NULL)
  }
  if (!S7_inherits(space, simplex_region)) {
    return(NULL)
  }
  v <- space@vertices
  # Membership first. The count branch below reads the vertex deficit as a
  # dimension deficit, which is only true once the vertices are known to share
  # the hyperplane `sum(theta) == 1`; a tetrahedron in `R^3` has four
  # affinely independent vertices and is not lower-dimensional at all.
  if (any(v < -1e-12)) {
    return(paste0(
      "its vertices leave the standard simplex: the smallest coordinate is ",
      format(min(v))
    ))
  }
  sums <- colSums(v)
  if (max(abs(sums - 1)) >= 1e-9) {
    return(paste0(
      "its vertices leave the standard simplex: the coordinates of one sum ",
      "to ",
      format(sums[which.max(abs(sums - 1))]),
      " rather than 1"
    ))
  }
  if (ncol(v) != nrow(v)) {
    return(paste0(
      "it has ",
      ncol(v),
      " vertices in ",
      nrow(v),
      " dimensions, so it is a simplex of dimension ",
      ncol(v) - 1L,
      " inside a parameter space of dimension ",
      nrow(v) - 1L
    ))
  }
  paste0(
    "it is too ill-conditioned to enclose: the reciprocal condition number ",
    "of its edge matrix is ",
    format(simplex_rcond(v))
  )
}


#' Check that a bounding method returned what [certify()] needs
#'
#' This enforces a contract between `certify()` and the bounding methods.
#' @keywords internal
#' @noRd
check_bound_result <- function(results, method_name, n_parts) {
  if (!is.list(results) || length(results) != n_parts) {
    stop(
      "The `",
      method_name,
      "` bounding method returned ",
      length(results),
      " results for ",
      n_parts,
      " parts.",
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
#' Each part is reparametrised onto its own simplex and enclosed separately.
#' @keywords internal
#' @noRd
bernstein_bound <- function(x, family, parts, control) {
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
  lapply(parts, function(s) {
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


#' Fetch a method that can certify a (family, cell) combination, or `NULL`
#'
#' `cell` is one convex piece of a null, not the whole [region]: a null whose
#' cells differ in geometry resolves a method per cell.
#' @keywords internal
#' @noRd
certify_method <- function(family, cell) {
  for (method in certify_methods()) {
    fits <- S7_inherits(family, method$family) && method$accepts(cell)
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


#' Explain that no implemented method covers this (family, cell) combination
#'
#' Says what is missing rather than what is impossible. Certifying some other
#' pairing of family and geometry is a matter of deriving and implementing a
#' bound, not of the thing being unbounded.
#'
#' Two different refusals share this function. A [polytope_region()] or a
#' [halfspace_region()] has no method for its *geometry*, and naming the class
#' says everything. A [simplex_region()] that the Bernstein enclosure cannot
#' take is a different matter: the geometry is supported, this particular region
#' is the wrong shape for it, and naming the class would send the reader looking
#' for a missing method that is not missing. Those get the failing condition
#' instead, and a reminder that only certification is affected.
#'
#' The second wording is reserved for families some method actually claims. A
#' region obstruction is not the reason a `gaussian_family` fails to certify,
#' whatever shape its region happens to be.
#' @keywords internal
#' @noRd
unimplemented_message <- function(family, cell) {
  geometry <- class_name(cell)
  # A region obstruction is only the reason when some method claims this
  # family. Otherwise the family is what is missing, and naming the region's
  # shape sends the reader after the wrong thing: a `gaussian_family` over a
  # flat `simplex_region` would be told to fix the region, when a
  # full-dimensional one would not certify either.
  claims_family <- any(vapply(
    certify_methods(),
    \(m) S7_inherits(family, m$family),
    logical(1)
  ))
  obstruction <- if (claims_family) bernstein_obstruction(cell) else NULL
  if (!is.null(obstruction)) {
    return(paste0(
      "The Bernstein enclosure cannot bound ",
      class_name(family),
      " expectations over this ",
      geometry,
      ", because ",
      obstruction,
      ".\n",
      "The enclosure pushes the Bernstein basis on the standard simplex onto ",
      "the region's vertices, which needs one vertex per coordinate, all of ",
      "them in that simplex. Only certification is affected: the region ",
      "charts, projects and fits like any other, and `sup_lb()` still ",
      "searches it."
    ))
  }
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
      as.vector(crossprod(
        score(family, theta, enumerate_space(family@sample_space)),
        weight
      ))
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
#' Cheap by comparison, and defined wherever the search is. An unbounded part
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
#'     simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' X <- likelihood(fam(c(0.4, 0.35, 0.25)))
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
    parts(null@region),
    function(s) {
      maximise_over(s, obj, n_seeds = n_seeds, n_restarts = n_restarts)
    }
  )
  best <- which.max(vapply(found, function(f) f$value, numeric(1)))
  list(
    sup_lb = found[[best]]$value,
    theta = found[[best]]$theta,
    part = best
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
node_table <- function(result, part) {
  nodes <- result$history
  if (!length(nodes)) {
    return(NULL)
  }
  field <- function(name, template) {
    vapply(
      nodes,
      function(b) {
        if (name %in% names(b)) {
          b[[name]]
        } else {
          template
        }
      },
      template
    )
  }
  data.frame(
    part = as.integer(part),
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
#' of the certificate. One row per leaf of the search, per part.
#'
#' The nodes live at any iteration tile their part exactly, so at `K = 3` the
#' `vertices` column draws the partition of the facet directly, in barycentric
#' coordinates, at every step and not merely at the end.
#' `depth` against `born` shows the shape of the tree: a max-bound queue rule
#' can produce a chain rather than anything balanced, which is visible here and
#' nowhere else.
#'
#' Only bounding methods that use branch and bound populate this. A method that
#' bounds in closed form contributes no rows.
#' @inheritParams certify
#' @return A data frame with `part`, `id`, `parent`, `depth`, `born`,
#'   `retired`, `fate`, `upper`, `volume` and a `vertices` list column, plus the
#'   certificate itself in the `"certificate"` attribute and the per-iteration
#'   bound in `"trace"`.
#' @seealso [certify()]
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' X <- likelihood(fam(c(0.4, 0.35, 0.25)))
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
#' Certification is refused for any (family, cell) combination with no known
#' bounding method.
#'
#' Currently only multinomial families with simplex parts are supported. The
#' upper bound is found via a branch-and-bound algorithm that recursively
#' subdivides the simplex \insertCite{Leroy2012}{ripr}, and bounds each subset
#' via the simplicial Bernstein range enclosure property
#' \insertCite{Garloff1986}{ripr}.
#'
#' Unlike [sup_lb()] this is a proven bound.
#' @param x A [random_variable].
#' @param null A [null_model].
#' @param tol Stop once the bound is within `tol` of the best value found.
#' @param max_nodes Cap on subdivisions per part for branch-and-bound
#'   algorithms.
#' @param max_coefficients Refuse above this many Bernstein coefficients
#'   (for bounding multinomial expectation in simplices).
#' @param .record Also return branch-and-bound nodes. Use [certify_trace()]
#'   rather than this directly.
#' @return A list with `sup_ub`, `sup_lb`, the `random_variable` and `null` it
#'   holds for, the `method` names that produced it (one per distinct part
#'   geometry), and per-part `bounds`, `incumbents`, `iterations`,
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
#'     simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' X <- likelihood(fam(c(0.4, 0.35, 0.25)))
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
  # Get applicable bound method for each part.
  methods <- lapply(
    parts(null@region),
    function(s) certify_method(family, s)
  )
  # Stop if any (family, cell) combination is not implemented. Report the
  # offending parts rather than the null model, and deduplicate: a plurality
  # null has one part per candidate and they share a geometry, so the same
  # message would otherwise repeat K - 1 times.
  unavailable <- vapply(methods, is.null, logical(1L))
  if (any(unavailable)) {
    unimpl_msgs <- vapply(
      parts(null@region)[unavailable],
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

  # Group the parts by resolved method, run each bound_fn once on its group,
  # then put the results back in part order. With one entry in the registry
  # this is a single group; the grouping is what stops that being an assumption.
  control <- list(
    tol = tol,
    max_nodes = max_nodes,
    max_coefficients = max_coefficients
  )
  per_part <- vector("list", n_parts(null@region))
  for (name in method_names) {
    which_parts <- which(
      vapply(methods, function(m) m$name, character(1L)) == name
    )
    method <- methods[[which_parts[[1L]]]]
    results <- method$bound_fn(
      x,
      family,
      parts(null@region)[which_parts],
      control
    )
    check_bound_result(results, name, length(which_parts))
    per_part[which_parts] <- results
  }

  bounds <- vapply(per_part, function(r) r$bound, numeric(1L))
  incumbents <- vapply(per_part, function(r) r$incumbent, numeric(1L))
  iterations <- vapply(per_part, function(r) r$iterations, integer(1L))
  converged <- vapply(per_part, function(r) r$converged, logical(1L))
  budget_hit <- vapply(per_part, function(r) r$budget_hit, logical(1L))
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
    # Records the tables of the node histories for each part. Used by
    # `certify_trace`.
    tables <- lapply(
      seq_along(per_part),
      function(i) node_table(per_part[[i]], i)
    )
    out$record <- do.call(rbind, Filter(Negate(is.null), tables))
    out$traces <- lapply(per_part, function(r) r$trace)
    out$incumbent_traces <- lapply(per_part, function(r) r$incumbent_trace)
  }
  out
}
