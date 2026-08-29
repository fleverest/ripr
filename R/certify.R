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
#' Each entry decides for itself which `(cell, family)` combinations it can
#' bound, and carries the `bound_fn` that does it. A combination no entry
#' claims simply has no implementation yet.
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
#' `bound_fn(x, family, cells, control)` returns one result per element of
#' `cells`, in the same order, each a list with the fields
#' `check_bound_result()` requires.
#'
#' Named `bound_fn` rather than `engine` because `ripr_engine` already means a
#' quadrature rule for expectations under the alternative ([gh_engine()],
#' [mc_engine()]). That is a different object doing a different job, and the two
#' meet in the same conversations often enough for one word to serve both.
#'
#' A `bound_fn` receives a *group* of cells that work for a given `(x, family)`,
#' so that enumerating the sample space, evaluating `x` on it, building the
#' lattice (for multinomial), happens just once for all matching cells.
#' `certify()` groups the cells by resolved method, so a `bound_fn` only
#' ever sees geometries it facilitates, and a null with cells of differing
#' geometries is split across bounding methods rather than refused.
#'
#' # The shared incumbent
#'
#' `control$incumbent` is the largest value attained anywhere on the null so
#' far: `-Inf` for the first group, and the best of the preceding groups'
#' `incumbent` fields after that.
#'
#' A method that ignores the field is correct, just slower. What it costs is
#' that the order the groups run in decides how much work each does, but the
#' final supremum returned should be the same.
#'
#' @return A list of methods, each with `name`, `bound_fn`, a `subject` naming
#'   the method as the subject of a sentence, and a `fit`.
#'
#' `fit(cell, family)` answers three questions at once, for one convex cell (a
#' [convex_region], not a whole [region]) and the family it would be bounded
#' under:
#'
#' \describe{
#'   \item{`TRUE`}{this method can bound that combination.}
#'   \item{a list with `because` and `remedy`}{this method is *about* that
#'   combination but cannot proceed; the elements describe why and what to do.}
#'   \item{`NULL`}{not this method's business, so it has nothing to say.}
#' }
#' @keywords internal
#' @noRd
certify_methods <- function() {
  list(
    list(
      name = "point",
      subject = "Evaluation at a point",
      fit = point_fit,
      bound_fn = point_bound
    ),
    list(
      name = "bernstein",
      subject = "The Bernstein enclosure",
      fit = bernstein_fit,
      bound_fn = bernstein_bound
    )
  )
}


#' Can a supremum over this cell be found by evaluation alone?
#'
#' A [point_region] is a single parameter, so the supremum over it is the
#' expectation at that parameter and there is nothing more to enclose. That
#' makes the method independent of the family in a way no other bound is: it
#' does not require a polynomial form, there are no vertices and no subdivision
#' is required.
#'
#' But the expectation has to be computable *precisely*, which means summing
#' over an enumerable sample space (for now; TODO?), and not through monte carlo
#' or quadrature.
#'
#' Also the parameter must belong to the families parameter space.
#' @keywords internal
#' @noRd
point_fit <- function(cell, family) {
  if (!S7_inherits(cell, point_region)) {
    return(NULL)
  }
  if (!is_finite_space(family@sample_space)) {
    return(point_unenumerable(family))
  }
  if (!contains(family@parameter_space, cell@theta)) {
    return(NULL)
  }
  TRUE
}


#' The one refusal point evaluation owns
#' @keywords internal
#' @noRd
point_unenumerable <- function(family) {
  list(
    because = paste0(
      "its expectation is an integral over a `",
      class_name(family@sample_space),
      "` rather than a sum over an enumerable one"
    ),
    remedy = paste0(
      "The expectation has to be computed exactly, i.e. not via quadrature or",
      "monte carlo estimates, for certification. Only certification is ",
      "affected: the region charts, projects and fits like any other, and ",
      "`sup_lb()` still searches it."
    )
  )
}


#' Can the Bernstein enclosure bound this cell, and if not, why not?
#'
#' The family gate lives here and only here, which is what entitles everything
#' below it to talk about the standard simplex.
#' @keywords internal
#' @noRd
bernstein_fit <- function(cell, family) {
  if (!S7_inherits(family, multinomial_family)) {
    return(NULL)
  }
  if (bernstein_compatible(cell)) {
    return(TRUE)
  }
  bernstein_obstruction(cell)
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
  if (ncol(edges) == 0L) {
    # A single vertex has no edges, so it is perfect fine
    return(1)
  }
  sv <- svd(edges, nu = 0L, nv = 0L)$d
  if (sv[1L] <= 0) {
    return(0)
  }
  sv[length(sv)] / sv[1L]
}


#' Why the Bernstein enclosure cannot handle this region, or `NULL` if it can
#'
#' Names the condition that fails *and* what you can do about it. There are
#' four, and they are four different kinds of problem: an unbounded region,
#' a region outside the family's parameter space, a simplex of less than full
#' dimension (TODO?), and a simplex too narrow.
#'
#' A bounded region that is not a simplex is not among them, because [cells()]
#' triangulated it before it got here. `NULL` covers that case and anything
#' else unforeseen, and the caller falls back to naming the class.
#' @param space A [convex_region].
#' @return `NULL`, or a list with `because` (the clause completing "... cannot
#'   bound this region, because ...") and `remedy` (the whole of what the
#'   reader should take from it, the caller appending nothing). Three of the
#'   four remedies close by saying that only certification is affected; the
#'   region outside the parameter space does not, because for that one it is
#'   not true.
#' @keywords internal
#' @noRd
bernstein_obstruction <- function(space) {
  if (bernstein_compatible(space)) {
    return(NULL)
  }
  only_cert_affected_msg <- paste0(
    "Only certification is affected: the region ",
    "charts, projects and fits like any other, and `sup_lb()` still ",
    "searches it."
  )
  # Before the class test, because it is the class-independent one: a
  # `halfspace_region` and an `real_region` fail for the same reason,
  # and so would an unbounded `polyhedron_region` belonging to neither.
  if (!is_bounded(space)) {
    return(list(
      because = paste0(
        "it is unbounded, and no finite set of simplices covers an unbounded ",
        "region"
      ),
      remedy = paste0(
        "The Bernstein enclosure only works for bounded polytopes. State the ",
        "null over a bounded region instead. ",
        only_cert_affected_msg
      )
    ))
  }
  # Every bounded cell that reaches here is a `polytope_region`: a
  # `simplex_region` from the fan, or a `point_region`, which is the fan's
  # degenerate output and a `polytope_region` too.
  if (!S7_inherits(space, polytope_region)) {
    return(NULL)
  }
  v <- space@vertices
  # Membership first. The count branch below reads the vertex deficit as a
  # dimension deficit, which is only true once the vertices are known to share
  # the hyperplane `sum(theta) == 1`; a tetrahedron in `R^3` has four
  # affinely independent vertices and is not lower-dimensional at all.
  outside <- paste0(
    "A region reaching outside the standard simplex cannot be bounded by its ",
    "Bernstein coefficients, as the Bernstein basis polynomials may take ",
    "negative values there."
  )
  if (any(v < -1e-12)) {
    return(list(
      because = paste0(
        "its vertices leave the standard simplex: the smallest coordinate is ",
        format(min(v))
      ),
      remedy = outside
    ))
  }
  sums <- colSums(v)
  if (max(abs(sums - 1)) >= 1e-9) {
    return(list(
      because = paste0(
        "its vertices leave the standard simplex: the coordinates of one sum ",
        "to ",
        format(sums[which.max(abs(sums - 1))]),
        " rather than 1"
      ),
      remedy = outside
    ))
  }
  if (ncol(v) != nrow(v)) {
    return(list(
      because = paste0(
        "it has ",
        count_label(ncol(v), "vertex", "vertices"),
        " in ",
        nrow(v),
        " dimensions, so it is a simplex of dimension ",
        ncol(v) - 1L,
        " inside a parameter space of dimension ",
        nrow(v) - 1L
      ),
      remedy = paste0(
        "The enclosure pushes a degree-`n_trials` Bernstein lattice on the ",
        "standard simplex onto a simplex of the same dimension, so it needs ",
        "one vertex per coordinate. A lower-dimensional region would need a ",
        "lattice of its own dimension instead, for which a parametrisation is ",
        "not yet implemented. ",
        only_cert_affected_msg
      )
    ))
  }
  list(
    because = paste0(
      "it is too ill-conditioned to enclose: the reciprocal condition number ",
      "of its edge matrix is ",
      format(simplex_rcond(v))
    ),
    remedy = paste0(
      "The region is a genuine simplex but the enclosure inverts the vertex ",
      "matrix to reparametrise onto it, and an inversion this ill-conditioned ",
      "cannot be trusted to produce a valid bound. This is a numerical limit ",
      "rather than a shape the method excludes. ",
      only_cert_affected_msg
    )
  )
}


#' Check that a bounding method returned what [certify()] needs
#'
#' This enforces a contract between `certify()` and the bounding methods.
#' @keywords internal
#' @noRd
check_bound_result <- function(results, method_name, n_cells) {
  if (!is.list(results) || length(results) != n_cells) {
    stop(
      "The `",
      method_name,
      "` bounding method returned ",
      length(results),
      " results for ",
      n_cells,
      " cells.",
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
#' Each cell is reparametrised onto its own simplex and enclosed separately.
#' @keywords internal
#' @noRd
bernstein_bound <- function(x, family, cells, control) {
  check_bernstein_size(family@n_trials, family@k, control$max_coefficients)

  values <- evaluate_on_space(x, family)
  lattice <- bernstein_lattice(family@n_trials, family@k)
  boxes <- lapply(cells, function(s) {
    list(V = s@vertices, coef = reparametrise_to(values, lattice, s@vertices))
  })
  incumbent <- max(control$incumbent, boxes_best(boxes, lattice)$value)
  lapply(boxes, function(box) {
    result <- certify_sup(
      list(box),
      lattice,
      tol = control$tol,
      max_iter = control$max_nodes,
      shared_incumbent = incumbent
    )
    incumbent <<- max(incumbent, result$incumbent)
    result
  })
}


#' Evaluate a random variable on the whole of an enumerable sample space
#'
#' Raises an error when the random variable is not bounded.
#' @keywords internal
#' @noRd
evaluate_on_space <- function(x, family) {
  values <- x(enumerate_space(family@sample_space))
  if (any(!is.finite(values))) {
    stop(
      "Cannot certify: the variable is not finite everywhere on the sample ",
      "space, so its null expectation is unbounded.",
      call. = FALSE
    )
  }
  values
}


#' Exact evaluation at a single parameter, for any enumerable family
#'
#' \eqn{\sup_{\theta \in \{\theta_0\}} E_\theta[X]}{sup over {theta_0} of
#' E_theta[X]} is \eqn{E_{\theta_0}[X]}{E_theta0[X]}, so there is no
#' enclosure, no subdivision and no budget: the answer is one weighted sum over
#' the sample space, and the search converges before it starts.
#'
#' The sum is evaluated in floating point, like every bound in the package.
#' The mathematics yields a guaranteed bound, but the arithmetic being IEEE
#' double means we do not yield a strictly *proven* bound.
#' @keywords internal
#' @noRd
point_bound <- function(x, family, cells, control) {
  outcomes <- enumerate_space(family@sample_space)
  values <- evaluate_on_space(x, family)
  loglik <- compile_loglik(family, outcomes)
  lapply(cells, function(cell) {
    log_p <- as.vector(loglik(matrix(cell@theta, ncol = 1L)))
    # A zero-probability outcome contributes an exact zero through a `-Inf`
    # log; `exp(-Inf) * x` is a clean 0 for finite `x`.
    value <- sum(exp(log_p) * values)
    list(
      bound = value,
      incumbent = value,
      theta = cell@theta,
      iterations = 0L,
      converged = TRUE,
      budget_hit = FALSE
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
    if (isTRUE(method$fit(cell, family))) {
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
#' @keywords internal
#' @noRd
unimplemented_message <- function(family, cell) {
  geometry <- class_name(cell)
  # A region obstruction is only the reason when some method claims this
  # family. Otherwise the family is what is missing, and naming the region's
  # shape sends the reader after the wrong thing: a `gaussian_family` over a
  # flat `simplex_region` would be told to fix the region, when a
  # full-dimensional one would not certify either.
  for (method in certify_methods()) {
    obstruction <- method$fit(cell, family)
    if (is.list(obstruction)) {
      return(paste0(
        method$subject,
        " cannot bound ",
        class_name(family),
        " expectations over this ",
        geometry,
        ", because ",
        obstruction$because,
        ".\n",
        obstruction$remedy
      ))
    }
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
#' treating it like a certificate on the bound. See [certify()] for global
#' upper bounds where supported.
#'
#' Cheap by comparison, and defined wherever the search is. An unbounded part
#' or a continuous parameter space can still be searched, though the lower bound
#' may sit far below the supremum.
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
#' One run is one cell, and `id` restarts at 1 in each, so `(cell, id)` is what
#' identifies a node and `parent` is to be matched within a cell. `part` is
#' along for reporting: it is `cell_part[cell]`, constant down a cell's rows,
#' and a part with several cells contributes several trees rather than one.
#'
#' Coefficients are dropped. They are the bulk of a node (10,626 doubles each at
#' `K = 5, n = 20`, against a handful for everything else here) and nothing
#' downstream of a finished run evaluates them.
#' @keywords internal
#' @noRd
node_table <- function(result, cell, part) {
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
    cell = as.integer(cell),
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
#' of the certificate. One row per node of the search, per cell.
#'
#' The nodes live at any iteration tile their cell exactly, so at `K = 3` the
#' `vertices` column draws the partition of the facet directly, in barycentric
#' coordinates, at every step and not merely at the end.
#' `depth` against `born` shows the shape of the tree: a max-bound queue rule
#' can produce a chain rather than anything balanced, which is visible here and
#' nowhere else.
#'
#' A cell is one branch-and-bound run and `id` restarts at 1 in each, so group
#' by `cell` before reading `id` or matching `parent`. `part` says which of the
#' declared parts a cell came from, and a triangulated part contributes one
#' tree per cell rather than one between them; `traces` and `incumbent_traces`
#' are per cell on the same terms.
#'
#' Only bounding methods that use branch and bound populate this. A method that
#' bounds in closed form contributes no rows.
#' @inheritParams certify
#' @return A data frame with `part`, `cell`, `id`, `parent`, `depth`, `born`,
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
#' An upper bound on
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
#' bounding method.A part that is a [polytope_region()] is triangulated first,
#' so any bounded null will certify for a multinomial expectation. What is
#' currently refused is a part that is unbounded (e.g. a [polyhedron_region()]
#' with rays or lineality, of which [halfspace_region()] and
#' [real_region()] are instances of).
#'
#' Only two methods are currently implemented. A [point_region()] is certified
#' by evaluation, for any family whose sample space can be enumerated (and thus
#' the expectation computed precisely) the supremum over a single parameter is
#' the expectation at that point, so there is nothing to enclose. Anything
#' larger needs a specific bounding method, and the only one currently
#' implemented is Bernstein branch-and-bound for multinomial families over
#' simplices, found via a branch-and-bound algorithm that recursively subdivides
#' the simplex \insertCite{Leroy2012}{ripr} and bounds each subset via the
#' simplicial Bernstein range enclosure property \insertCite{Garloff1986}{ripr}.
#'
#' Unlike [sup_lb()] this is a bound on the global supremum rather than a local
#' optimum.
#'
#' ## Numerical limitations
#'
#' The derived bounds are mathematically guaranteed, but evaluated in floating
#' point arithmetic. The geometry underneath (triangulation, set algebra,
#' emptiness) is exact via GMP rationals, so the cells genuinely tile the null
#' and share facets with no gaps or overlaps. Still, the bounding arithmetic
#' itself (Bernstein coefficients, point evaluation) uses ordinary IEEE double
#' precision, and no accounting is made for its rounding, so a certificate here
#' is a mathematical bound computed in floating point, not a formally proven
#' one.
#' @param x A [random_variable].
#' @param null A [null_model].
#' @param tol Stop once the bound is within `tol` of the best value found.
#' @param max_nodes Cap on subdivisions *per cell* for branch-and-bound
#'   algorithms, so a part triangulated into several cells is allowed
#'   `max_nodes` in each. The per-part `iterations` reported below is the total
#'   actually spent, which is what to read the cap against.
#' @param max_coefficients Refuse above this many Bernstein coefficients
#'   (for bounding multinomial expectation in simplices).
#' @param .record Also return branch-and-bound nodes. Use [certify_trace()]
#'   rather than this directly.
#' @return A list with `sup_ub`, `sup_lb`, the `random_variable` and `null` it
#'   holds for, the `method` names that produced it (one per distinct cell
#'   geometry), and per-part `bounds`, `incumbents`, `iterations`,
#'   `converged` and `budget_hit`, each reduced over the part's cells: the
#'   largest bound and incumbent, the total iterations, converged only if every
#'   cell did and `budget_hit` if any did. `budget_hit` flags when a search
#'   stopped at `max_nodes` with the gap still open, so its bound is valid but
#'   likely loose.
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
  # A bound is derived on one cell at a time.
  cells <- null@cells
  cell_part <- null@cell_part
  methods <- lapply(cells, function(s) certify_method(family, s))
  # Stop if any (family, cell) combination is not implemented. Report the
  # offending cells rather than the null model, and deduplicate: a plurality
  # null has one cell per candidate and they share a geometry, so the same
  # message would otherwise repeat K - 1 times.
  unavailable <- vapply(methods, is.null, logical(1L))
  if (any(unavailable)) {
    unimpl_msgs <- vapply(
      cells[unavailable],
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

  # Group the cells by resolved method, run each bound_fn once on its group,
  # then put the results back in cell order.
  control <- list(
    tol = tol,
    max_nodes = max_nodes,
    max_coefficients = max_coefficients,
    incumbent = -Inf
  )
  per_cell <- vector("list", length(cells))
  for (name in method_names) {
    which_cells <- which(
      vapply(methods, function(m) m$name, character(1L)) == name
    )
    method <- methods[[which_cells[[1L]]]]
    results <- method$bound_fn(x, family, cells[which_cells], control)
    check_bound_result(results, name, length(which_cells))
    per_cell[which_cells] <- results
    # Carry the incumbent from this group into onto the next group. The
    # incumbent from the previous run bounds the supremum from below for all
    # future steps, so the next groups may use it to prune more aggressively.
    control$incumbent <- max(
      control$incumbent,
      vapply(results, function(r) r$incumbent, numeric(1L))
    )
  }

  # Reduce the cells back to the parts the caller declared. A part is the union
  # of its cells, so its supremum is the largest of theirs and so is the best
  # value attained in any of them; the work spent is the total, and a part has
  # converged only if all of its cells did.
  by_part <- split(
    seq_along(cells),
    factor(
      cell_part,
      seq_len(n_parts(
        null@region
      ))
    )
  )
  reduce <- function(field, combine, template) {
    per <- vapply(per_cell, function(r) r[[field]], template)
    vapply(by_part, function(i) combine(per[i]), template)
  }
  bounds <- unname(reduce("bound", max, numeric(1L)))
  incumbents <- unname(reduce("incumbent", max, numeric(1L)))
  iterations <- unname(reduce("iterations", sum, integer(1L)))
  converged <- unname(reduce("converged", all, logical(1L)))
  budget_hit <- unname(reduce("budget_hit", any, logical(1L)))
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
    # Records the tables of the node histories for each cell. Used by
    # `certify_trace`. These stay per cell rather than being reduced: the tree
    # is what the recording is for, and two cells' trees do not combine into
    # one. The `part` column carries the reduction's grouping.
    tables <- lapply(
      seq_along(per_cell),
      function(i) node_table(per_cell[[i]], cell = i, part = cell_part[[i]])
    )
    out$record <- do.call(rbind, Filter(Negate(is.null), tables))
    out$traces <- lapply(per_cell, function(r) r$trace)
    out$incumbent_traces <- lapply(per_cell, function(r) r$incumbent_trace)
  }
  out
}
