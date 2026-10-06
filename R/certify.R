#' @include random_variable.R null.R bernstein.R
NULL

# Bounds on the expectation of a random variable under a null. `sup_lb()`
# searches for the supremum via local gradient ascent (yielding a lower bound),
# while `certify()` encloses the supremum over regions. Per cell, `certify()`
# uses evaluation for a `point_region` and the Bernstein enclosure otherwise
# (multinomial family, simplex in the standard one).

#' Why a cell cannot be certified under a family, or `NULL` if it can
#'
#' `NULL` also fixes the method: `"point"` for a [point_region], `"bernstein"`
#' otherwise. Messages say what is missing (a condition, or an unimplemented
#' method), never that the expectation is unbounded.
#' @keywords internal
#' @noRd
certify_obstruction <- function(cell, family) {
  is_point <- S7_inherits(cell, point_region)
  enumerable <- is_finite_space(family@sample_space)
  if (is_point && enumerable && contains(family@parameter_space, cell@theta)) {
    return(NULL)
  }
  why <- NULL
  if (is_point && !enumerable) {
    subject <- "Evaluation at a point"
    why <- point_unenumerable(family)
  } else if (S7_inherits(family, multinomial_family)) {
    # Blame the region only for a family the enclosure claims; otherwise the
    # family is what is missing.
    subject <- "The Bernstein enclosure"
    why <- bernstein_obstruction(cell)
    if (is.null(why) && !is_point) {
      return(NULL)
    }
  }
  geometry <- class_name(cell)
  if (!is.null(why)) {
    return(paste0(
      subject,
      " cannot bound ",
      class_name(family),
      " expectations over this ",
      geometry,
      ", because ",
      why$because,
      ".\n",
      why$remedy
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


#' Why evaluation at a point cannot certify under this family
#'
#' The expectation must be exact, i.e. a sum over an enumerable sample space,
#' not quadrature or Monte Carlo. A `theta` outside the parameter space is left
#' to the family to explain. Returns `list(because, remedy)`.
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
      "The expectation has to be computed exactly, i.e. not via quadrature or ",
      "Monte Carlo estimates, for certification. Only certification is ",
      "affected: the region charts, projects and fits like any other, and ",
      "`sup_lb()` still searches it."
    )
  )
}


#' Why the Bernstein enclosure cannot handle this region, or `NULL` if it can
#'
#' `NULL` exactly when `reparametrise_to()` will accept the vertices: a simplex
#' inside the standard simplex. Returns `list(because, remedy)`, where `because`
#' completes "... cannot bound this region, because ...", `remedy` is complete
#' and says only certification is affected except when the region leaves the
#' parameter space.
#' @keywords internal
#' @noRd
bernstein_obstruction <- function(space) {
  only_cert_affected_msg <- paste0(
    "Only certification is affected: the region ",
    "charts, projects and fits like any other, and `sup_lb()` still ",
    "searches it."
  )
  # Before the class tests: it applies to any unbounded region.
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
  not_simplex <- list(
    because = paste0("it is a `", class_name(space), "` rather than a simplex"),
    remedy = paste0(
      "The enclosure reparametrises onto a simplex's vertices, and `cells()` ",
      "triangulates every bounded polytope into simplices before it gets ",
      "here, so this region is one it could not triangulate. ",
      only_cert_affected_msg
    )
  )
  # Bounded cells arrive triangulated: `simplex_region` or `point_region`.
  if (!S7_inherits(space, polytope_region)) {
    return(not_simplex)
  }
  v <- space@vertices
  departure <- simplex_departure(v)
  if (!is.null(departure)) {
    return(list(
      because = paste0("its vertices leave the standard simplex: ", departure),
      remedy = paste0(
        "A region reaching outside the standard simplex cannot be bounded by ",
        "its Bernstein coefficients, as the Bernstein basis polynomials may ",
        "take negative values there."
      )
    ))
  }
  if (ncol(v) < 2L) {
    return(list(
      because = "it has a single coordinate, so its simplex is one point",
      remedy = paste0(
        "A one-category multinomial has a single outcome, and its expectation ",
        "is that outcome's value at every parameter; there is nothing to bound."
      )
    ))
  }
  if (!S7_inherits(space, simplex_region)) {
    return(not_simplex)
  }
  NULL
}


#' Bernstein enclosure over simplices, for multinomial expectations
#'
#' The multinomial pmf is the Bernstein basis, so `x` on the sample space is
#' already the coefficient vector. All `cells` share one lattice and one
#' evaluation of `x`. `incumbent` (the best value attained on the null so far)
#' is raised as the runs go and each prunes against it.
#' @keywords internal
#' @noRd
bernstein_bound <- function(
  x,
  family,
  cells,
  tol,
  max_splits,
  max_coefficients,
  incumbent = -Inf
) {
  check_bernstein_size(family@n_trials, family@k, max_coefficients)

  values <- evaluate_on_space(x, family)
  lattice <- bernstein_lattice(family@n_trials, family@k)
  boxes <- lapply(cells, function(s) {
    V <- pad_vertices(s@vertices, lattice$K)
    list(V = V, coef = reparametrise_to(values, lattice, V))
  })
  incumbent <- max(incumbent, boxes_best(boxes, lattice)$value)
  lapply(boxes, function(box) {
    result <- certify_sup(
      box,
      lattice,
      tol = tol,
      max_iter = max_splits,
      shared_incumbent = incumbent
    )
    incumbent <<- max(incumbent, result$incumbent)
    result
  })
}


#' `x` on the whole enumerable sample space; errors unless finite everywhere
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
#' The supremum over `{theta_0}` is `E_theta0[X]`: one weighted sum, no
#' enclosure. Returns the `certify_sup()` fields [certify()] reduces, per cell.
#' @keywords internal
#' @noRd
point_bound <- function(x, family, cells) {
  outcomes <- enumerate_space(family@sample_space)
  values <- evaluate_on_space(x, family)
  loglik <- compile_loglik(family, outcomes)
  lapply(cells, function(cell) {
    log_p <- as.vector(loglik(matrix(cell@theta, nrow = 1L)))
    # `exp(-Inf) * x` is a clean 0 for finite `x`.
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


#' The value(s) that `certify()`'s `incumbent_at` seeds the search with
#' @keywords internal
#' @noRd
seed_incumbents <- function(x, null, at) {
  regions <- parts(null@region)
  if (is.null(at)) {
    return(rep(-Inf, length(regions)))
  }
  if (is.numeric(at) && is.null(dim(at))) {
    at <- matrix(at, nrow = 1L)
  }
  d <- space_dim(null@region)
  if (!is.matrix(at) || !is.numeric(at) || ncol(at) != d || anyNA(at)) {
    stop(
      "`incumbent_at` must be a parameter vector of length ",
      d,
      ", or a matrix of them with one per row.",
      call. = FALSE
    )
  }
  points <- lapply(seq_len(nrow(at)), function(i) at[i, ])
  inside <- vapply(
    points,
    function(theta) vapply(regions, contains, logical(1L), theta = theta),
    logical(length(regions))
  )
  inside <- matrix(inside, nrow = length(regions))
  outside <- which(colSums(inside) == 0L)
  if (length(outside)) {
    stop(
      "`incumbent_at` must lie in the null, but row ",
      paste(outside, collapse = ", "),
      " lies in no part of it.",
      call. = FALSE
    )
  }
  values <- vapply(
    point_bound(
      x,
      null@family,
      lapply(points, \(theta) point_region(theta = theta))
    ),
    function(r) r$incumbent,
    numeric(1L)
  )
  if (any(!is.finite(values))) {
    stop(
      "`incumbent_at` row ",
      paste(which(!is.finite(values)), collapse = ", "),
      " gives a non-finite expectation, so it cannot seed the search.",
      call. = FALSE
    )
  }
  vapply(
    seq_along(regions),
    function(p) max(values[inside[p, ]], -Inf),
    numeric(1L)
  )
}


#' `E_theta[X]` and its gradient as an [objective()], under the rule `spec`
#' resolves to at `P_theta`; on log scale when `x` has a log form.
#' @keywords internal
#' @noRd
expectation_objective <- function(family, x, spec) {
  log_scale <- has_log_form(x)
  terms <- memoise_last(function(theta) {
    engine <- resolve_engine(spec, family(theta), family)
    nodes <- engine@nodes
    if (log_scale) {
      log_terms <- engine@log_w + x@log_f(nodes)
      value <- matrixStats::logSumExp(log_terms)
      weight <- if (is.finite(value)) {
        exp(log_terms - value)
      } else {
        numeric(length(log_terms))
      }
    } else {
      weight <- exp(engine@log_w) * x(nodes)
      value <- sum(weight)
    }
    list(value = value, weight = weight, nodes = nodes)
  })

  objective(
    value = function(theta) terms(theta)$value,
    grad = function(theta) {
      at <- terms(theta)
      as.vector(crossprod(score(family, theta, at$nodes), at$weight))
    }
  )
}


# --- Result objects -----------------------------------------------------------

#' A certified upper bound on the largest null expectation
#'
#' A certificate returned by [certify()]: a certified bound
#' \eqn{\sup_{\theta \in \Theta_0} E_\theta[X] \le}{sup_theta E_theta[X] <=}
#' `sup_ub`, with the random variable, the null and how the bound was reached.
#' Passed to [e_variable()], this certificate can be used to yield an e-variable
#' for `null`. `sup_lb` is the largest value attained, so (`sup_lb`, `sup_ub`)
#' encloses the supremum.
#'
#' Per-part properties have one entry per declared part of the null's region,
#' reduced over that part's cells.
#'
#' @param sup_ub The certified upper bound on the supremum.
#' @param sup_lb The largest expected value found anywhere on the null, a lower
#'   bound on the same supremum.
#' @param random_variable The [random_variable] the bound is for.
#' @param null The [null_model] the bound holds over.
#' @param method Names of the bounding methods that produced it, one per
#'   distinct cell geometry.
#' @param bounds,incumbents Per-part upper bounds and attained values.
#' @param iterations Per-part iterations spent, in total over each part's cells.
#' @param converged Per-part: did the search end with no cell left that could
#'   raise `sup_ub` by more than `tol`, rather than on the budget?
#' @param budget_hit Per-part: did any cell stop at `max_splits` with its gap
#'   still open? The bound is still valid, but likely loose.
#' @return A `ripr_certificate`.
#' @seealso [certify()], [e_variable()], and [ripr_search] for what the
#'   searched lower bound [sup_lb()] returns.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' null <- null_model(
#'   fam,
#'   simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1)))
#' )
#' cert <- certify(likelihood(fam(c(0.4, 0.35, 0.25))), null)
#' cert
#' cert@sup_ub
#' cert@bounds
#' @export
ripr_certificate <- new_class(
  "ripr_certificate",
  properties = list(
    sup_ub = class_numeric,
    sup_lb = class_numeric,
    random_variable = random_variable,
    null = null_model,
    method = class_character,
    bounds = class_numeric,
    incumbents = class_numeric,
    iterations = class_integer,
    converged = class_logical,
    budget_hit = class_logical
  ),
  validator = function(self) {
    if (length(self@sup_ub) != 1L || length(self@sup_lb) != 1L) {
      return("`sup_ub` and `sup_lb` must each be a single number")
    }
    n <- length(self@bounds)
    per_part <- list(
      self@incumbents,
      self@iterations,
      self@converged,
      self@budget_hit
    )
    if (any(lengths(per_part) != n)) {
      return("the per-part properties must all have one entry per part")
    }
    NULL
  }
)


#' @description `print()` gives both bounds, the per-part account for up to
#'   eight parts, and what [e_variable()] will make of the certificate.
#' @rdname ripr_certificate
#' @usage NULL
#' @export
method(print, ripr_certificate) <- function(x, ...) {
  cat("<", class_name(x), ">\n", sep = "")
  cat("  X = ", format(x@random_variable), "\n", sep = "")
  cat("  under ", format(x@null), "\n", sep = "")
  cat(
    "  sup E[X] <= ",
    format(x@sup_ub, digits = 7L),
    "  (certified via ",
    paste(x@method, collapse = " and "),
    ")\n",
    sep = ""
  )
  cat(
    "  sup E[X] >= ",
    format(x@sup_lb, digits = 7L),
    "  (interval width: ",
    format(x@sup_ub - x@sup_lb, digits = 3L),
    ")\n",
    sep = ""
  )
  n <- length(x@bounds)
  if (n > 1L && n <= 8L) {
    table <- data.frame(
      bound = signif(x@bounds, 7L),
      attained = signif(x@incumbents, 7L),
      iterations = x@iterations,
      converged = x@converged,
      row.names = paste0("  part ", seq_len(n))
    )
    print(table)
  }
  if (any(x@budget_hit)) {
    cat(
      "  node budget reached in ",
      ngettext(sum(x@budget_hit), "part ", "parts "),
      toString(which(x@budget_hit)),
      ": the bound holds but is likely loose\n",
      sep = ""
    )
  }
  cat("  e_variable(): ", e_variable_label(x), "\n", sep = "")
  invisible(x)
}


#' @description `format()` gives the certified bound on one line.
#' @rdname ripr_certificate
#' @usage NULL
#' @export
method(format, ripr_certificate) <- function(x, ...) {
  sprintf(
    "%s: sup E[%s] <= %s over %s",
    class_name(x),
    format(x@random_variable),
    format(x@sup_ub, digits = 7L),
    parts_label(length(x@bounds))
  )
}


#' What `e_variable()` makes of a certificate, in words
#' @keywords internal
#' @noRd
e_variable_label <- function(certificate) {
  if (certificate@sup_ub <= 1) {
    "X unchanged, already an e-variable"
  } else {
    paste0("X / ", format(certificate@sup_ub, digits = 7L))
  }
}


#' A searched lower bound on the largest null expectation
#'
#' An estimate of the supremum, returned by [sup_lb()]. It is the largest
#' expected value found through a multi-start gradient-ascent over the null, and
#' the location that the local optimum was found. It yields a **lower** bound,
#' so `sup_lb > 1` shows `X` is *not* an e-variable, but nothing here can prove
#' that it is one; use [certify()] for that.
#'
#' @param sup_lb The largest expectation found.
#' @param log_sup_lb Its logarithm, computed directly when the random variable
#'   has a log form, so that it remains finite where `sup_lb` underflows.
#' @param theta The parameter attaining it.
#' @param part The part of the null's region `theta` lies in.
#' @param random_variable The [random_variable] searched over.
#' @param null The [null_model] searched.
#' @return A `ripr_search`.
#' @seealso [sup_lb()], [ripr_certificate]
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' null <- null_model(
#'   fam,
#'   simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1)))
#' )
#' found <- sup_lb(likelihood(fam(c(0.4, 0.35, 0.25))), null, n_seeds = 20L)
#' found
#' found@theta
#' @export
ripr_search <- new_class(
  "ripr_search",
  properties = list(
    sup_lb = class_numeric,
    log_sup_lb = class_numeric,
    theta = class_numeric,
    part = class_integer,
    random_variable = random_variable,
    null = null_model
  )
)


#' @description `print()` gives the value found and where it was attained.
#' @rdname ripr_search
#' @usage NULL
#' @export
method(print, ripr_search) <- function(x, ...) {
  cat("<", class_name(x), ">\n", sep = "")
  cat("  X = ", format(x@random_variable), "\n", sep = "")
  cat("  under ", format(x@null), "\n", sep = "")
  cat(
    "  sup E[X] >= ",
    format(x@sup_lb, digits = 7L),
    "  (searched, not certified)\n",
    sep = ""
  )
  cat(
    "  attained at theta = ",
    theta_label(x@theta),
    ", in part ",
    x@part,
    "\n",
    sep = ""
  )
  invisible(x)
}


#' @description `format()` gives the value found on one line.
#' @rdname ripr_search
#' @usage NULL
#' @export
method(format, ripr_search) <- function(x, ...) {
  sprintf(
    "%s: sup E[%s] >= %s at theta = %s",
    class_name(x),
    format(x@random_variable),
    format(x@sup_lb, digits = 7L),
    theta_label(x@theta)
  )
}


#' An e-variable from a certified bound
#'
#' Returns `X / sup_ub`, which has expectation at most 1 under every
#' distribution in `H0`. If `sup_ub <= 1`, `X` is returned unchanged. The
#' result keeps `X`'s log form.
#'
#' `X` must be non-negative (e.g. a likelihood ratio), but this is not checked.
#' A lower bound via [ripr_search] from [sup_lb()] is illegal!
#' @param x A [ripr_certificate], as returned by [certify()].
#' @param ... Unused, for methods.
#' @return A [random_variable] with expectation at most 1 under the null.
#' @seealso [certify()]
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = rbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' # Against a uniform null point that is not the RIPr, the likelihood ratio
#' # is not an e-variable as it stands, and the certificate says by how much.
#' Q <- fam(c(0.4, 0.35, 0.25))
#' X <- likelihood(Q) / likelihood(fam(c(1, 1, 1) / 3))
#' cert <- certify(X, plurality, tol = 1e-6)
#' cert@sup_ub
#'
#' E <- e_variable(cert)
#' E
#' E(c(2, 1, 1))
#' @export
e_variable <- new_generic("e_variable", "x", function(x, ...) {
  S7::S7_dispatch()
})


method(e_variable, ripr_certificate) <- function(x, ...) {
  bound <- x@sup_ub
  if (is.na(bound) || !is.finite(bound)) {
    stop(
      "the certified bound is ",
      format(bound),
      ", so no rescaling of `X` makes it an e-variable.",
      call. = FALSE
    )
  }
  if (bound <= 1) {
    return(x@random_variable)
  }
  x@random_variable / bound
}


#' @rdname e_variable
#' @usage NULL
method(e_variable, ripr_search) <- function(x, ...) {
  stop(
    "a `ripr_search` from `sup_lb()` is a lower bound on the null ",
    "expectation, and dividing by it guarantees nothing. Use `certify()` for ",
    "an upper bound to rescale by.",
    call. = FALSE
  )
}


#' Estimate the largest null expectation of a random variable
#'
#' Multi-start local ascent, giving a **lower** bound on the supremum: a larger
#' value may exist where it did not look. Use it for diagnosis; see [certify()]
#' for a global upper bound. Works on unbounded parts and continuous spaces,
#' though the bound may then sit far below the supremum.
#'
#' A variable with a log form (e.g. from [likelihood()], `*`, `/`, `+`) is
#' integrated in log space.
#' @param x A [random_variable].
#' @param null A [null_model].
#' @param engine An engine spec for the expectation under `P_theta`, e.g.
#'   [exact_engine()], [gh_engine()] or [mc_engine()].
#' @param n_seeds,n_restarts Resolution of the search.
#' @return A [ripr_search], with properties `sup_lb`, its logarithm
#'   `log_sup_lb`, the `theta` attaining it and the `part` that `theta` lies
#'   in, alongside the `random_variable` and `null` searched.
#' @seealso [certify()], [ripr_search]
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = rbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' X <- likelihood(fam(c(0.4, 0.35, 0.25)))
#' sup_lb(X, plurality)
#'
#' # A continuous sample space is searched the same way, under a rule that
#' # integrates over it. Every evaluation is a quadrature rule here, so the
#' # search is run at a lower resolution than the enumerable one above.
#' gaussian <- gaussian_family(d = 2L)
#' halfspace <- null_model(
#'   gaussian,
#'   halfspace_region(normal = c(1, -1), offset = 0)
#' )
#' Y <- likelihood(gaussian(c(0, 0)))
#' sup_lb(
#'   Y,
#'   halfspace,
#'   engine = gh_engine(n_nodes = 9L),
#'   n_seeds = 20L,
#'   n_restarts = 2L
#' )
#' @export
sup_lb <- function(
  x,
  null,
  engine = exact_engine(),
  n_seeds = 200L,
  n_restarts = 25L
) {
  if (!S7_inherits(x, random_variable)) {
    stop("`x` must be a `random_variable`.", call. = FALSE)
  }

  obj <- expectation_objective(null@family, x, engine)

  found <- lapply(
    parts(null@region),
    function(s) {
      maximise_over(s, obj, n_seeds = n_seeds, n_restarts = n_restarts)
    }
  )
  values <- vapply(found, function(f) f$value, numeric(1))
  best <- which.max(values)
  if (length(best) == 0L) {
    # Every value `NaN`, e.g. a ratio underflowing to `0 / 0`.
    stop(
      "the expectation is `NaN` at every part, so there is nothing to report. ",
      "A likelihood ratio evaluated where both of its densities have ",
      "underflowed is `0 / 0`.",
      call. = FALSE
    )
  }
  value <- values[[best]]
  log_scale <- has_log_form(x)
  ripr_search(
    sup_lb = if (log_scale) exp(value) else value,
    log_sup_lb = if (log_scale) value else suppressWarnings(log(value)),
    theta = as.numeric(found[[best]]$theta),
    part = as.integer(best),
    random_variable = x,
    null = null
  )
}


#' Flatten a run's branch-and-bound nodes into a table, without coefficients
#'
#' `id` restarts per cell, so `(cell, id)` identifies a node and `parent`
#' matches within a cell. `part` is `cell_part[cell]`, for reporting.
#' @keywords internal
#' @noRd
node_table <- function(result, cell, part) {
  nodes <- result$history
  if (!length(nodes)) {
    return(NULL)
  }
  field <- function(name, template) vapply(nodes, `[[`, template, name)
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
#' Same computation as [certify()], returning the branch-and-bound tree: one
#' row per node, per cell. The nodes live at any iteration tile their cell, so
#' at `K = 3` the `vertices` column (`(K, K)` matrices, one vertex per row)
#' draws the partition at every step.
#'
#' Each cell is a separate run and `id` restarts at 1 in each, so group by
#' `cell` before matching `parent`. A triangulated part contributes one tree per
#' cell. Point cells, bounded in closed form, contribute no rows.
#' @inheritParams certify
#' @return A data frame with `part`, `cell`, `id`, `parent`, `depth`, `born`,
#'   `retired`, `fate`, `upper`, `volume` and a `vertices` list column, plus the
#'   certificate itself in the `"certificate"` attribute and the per-iteration
#'   bound in `"trace"`. The certificate is the same [ripr_certificate] that
#'   [certify()] returns.
#' @seealso [certify()]
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = rbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
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
  max_splits = 20000L,
  max_coefficients = 1024^2,
  incumbent_at = NULL
) {
  run <- certify_run(
    x,
    null,
    tol = tol,
    max_splits = max_splits,
    max_coefficients = max_coefficients,
    incumbent_at = incumbent_at,
    record = TRUE
  )
  nodes <- run$record
  attr(nodes, "trace") <- run$traces
  attr(nodes, "incumbent_trace") <- run$incumbent_traces
  attr(nodes, "certificate") <- run$certificate
  nodes
}


#' Certify an upper bound on the largest null expectation
#'
#' A global upper bound on
#' \eqn{\sup_{\theta \in \Theta_0} E_\theta[X]}{sup_theta E_theta[X]}, where
#' available. For non-negative `X`, `X / bound` is then an e-variable for `H0`
#' (or `X` itself if the bound is at most 1); [e_variable()] makes that choice.
#'
#' Two methods are implemented. A [point_region()] is certified by exact
#' evaluation, for any family with an enumerable sample space. Otherwise,
#' multinomial families over simplices use branch and bound
#' \insertCite{Leroy2012}{ripr} on the simplicial Bernstein range enclosure
#' \insertCite{Garloff1986}{ripr}, subdividing by de Casteljau's algorithm
#' \insertCite{PrautzschBoehmPaluszny2002}{ripr}. Bounded polytopes are
#' triangulated first, so any bounded null certifies for a multinomial.
#' Anything else, including an unbounded part such as a [halfspace_region()]
#' or [real_region()], is refused with a message saying why.
#'
#' ## Speeding things up: seeding the incumbent
#'
#' Branch and bound algorithms prune nodes once their bounds cannot beat the
#' best value attained so far (the incumbent). So seeding with an incumbent
#' that is close to the supremum early can save a lot of unnecessary branching.
#'
#' A good value to seed from is the Frank--Wolfe oracle value; it is a direct
#' result of a local optimisation that seeks the supremum. However, to guarantee
#' that the incumbent is attained at a point in the null that is actively being
#' certified, we accept incumbents (via `incumbent_at`) only by the location
#' that attains it. So rather than passing the incumbent value directly, you
#' would pass a point \eqn{\theta}{theta} in the null for which
#' \eqn{E_\theta[X]}{E_theta[X]} is the incumbent; i.e. the argmax instead of
#' the max.
#'
#' A fit's Frank--Wolfe oracle measured its gap (the `gap_after_theta` of
#' the fit's last trace row), pass it to certify via `incumbent_at`:
#' \eqn{E_\theta[X]}{E_theta[X]} is evaluated there exactly before the search
#' starts, and is effective as a principled starting point for the incumbent.
#'
#' ## Numerical limitations
#'
#' The geometry (triangulation, set algebra) is exact in GMP rationals, so the
#' cells tile the null exactly. The bounding arithmetic is IEEE double with no
#' accounting for rounding, so a certificate is a mathematical bound computed in
#' floating point, not a formally proven one.
#' @param x A [random_variable].
#' @param null A [null_model].
#' @param tol Stop once the bound is within `tol` of the best value found.
#' @param max_splits Cap on branch-and-bound subdivisions *per cell*; a
#'   triangulated part gets `max_splits` in each of its cells.
#' @param max_coefficients Refuse above this many Bernstein coefficients.
#' @param incumbent_at (Optional) Points in the null at which to evaluate
#'   \eqn{E_\theta[X]} before the search, to seed it with the largest value:
#'   a parameter vector, or a matrix of them with one per row. Each point must
#'   lie in the null (checked via [contains()]). See "Speeding things up" above.
#' @return A [ripr_certificate]. Where `budget_hit` is set, the bound is valid
#'   but likely loose.
#' @seealso [e_variable()], [sup_lb()], [certify_trace()]
#' @references
#' \insertAllCited{}
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = rbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' X <- likelihood(fam(c(0.4, 0.35, 0.25)))
#' cert <- certify(X, plurality, tol = 1e-6)
#' cert
#' c(upper = cert@sup_ub, attained = cert@sup_lb)
#'
#' # The bound is below 1 here, so X comes back unchanged.
#' E <- e_variable(cert)
#' E
#'
#' # Seeding the search with the approximate supremum from `sup_lb`.
#' found <- sup_lb(X, plurality, n_seeds = 20L)
#' seeded <- certify(X, plurality, tol = 1e-6, incumbent_at = found@theta)
#' c(cert@sup_ub, seeded@sup_ub)
#' c(sum(cert@iterations), sum(seeded@iterations))
#' @export
certify <- function(
  x,
  null,
  tol = 1e-6,
  max_splits = 20000L,
  max_coefficients = 1024^2,
  incumbent_at = NULL
) {
  certify_run(
    x,
    null,
    tol = tol,
    max_splits = max_splits,
    max_coefficients = max_coefficients,
    incumbent_at = incumbent_at
  )$certificate
}


#' The work behind `certify()` and `certify_trace()`; with `record`, also the
#' per-cell node table and traces
#' @keywords internal
#' @noRd
certify_run <- function(
  x,
  null,
  tol,
  max_splits,
  max_coefficients,
  incumbent_at = NULL,
  record = FALSE
) {
  if (!S7_inherits(x, random_variable)) {
    stop("`x` must be a `random_variable`.", call. = FALSE)
  }
  rlang::check_number_decimal(tol, min = 0)
  rlang::check_number_whole(max_splits, min = 1, max = 2147483647)
  rlang::check_number_whole(max_coefficients, min = 1, max = 2147483647)

  family <- null@family
  cells <- null@cells
  cell_part <- null@cell_part
  # Deduplicated: a plurality null's cells share a geometry.
  obstructions <- unlist(lapply(cells, certify_obstruction, family = family))
  if (length(obstructions)) {
    stop(
      "Cannot certify:\n",
      paste(unique(obstructions), collapse = "\n"),
      call. = FALSE
    )
  }
  is_point <- vapply(cells, S7_inherits, logical(1L), point_region)
  methods <- unique(ifelse(is_point, "point", "bernstein"))

  # Points first: they are cheap and give the Bernstein runs an incumbent.
  per_cell <- vector("list", length(cells))
  incumbent <- -Inf
  if (any(is_point)) {
    per_cell[is_point] <- point_bound(x, family, cells[is_point])
    incumbent <- max(vapply(per_cell[is_point], function(r) r$incumbent, 0))
  }
  # Then the caller's seeds, by the same exact evaluation.
  seeded <- seed_incumbents(x, null, incumbent_at)
  incumbent <- max(incumbent, seeded)
  if (!all(is_point)) {
    per_cell[!is_point] <- bernstein_bound(
      x,
      family,
      cells[!is_point],
      tol = tol,
      max_splits = max_splits,
      max_coefficients = max_coefficients,
      incumbent = incumbent
    )
  }

  # Reduce cells to declared parts: max bound and incumbent, total iterations,
  # converged if all cells did, budget_hit if any did.
  by_part <- split(
    seq_along(cells),
    factor(cell_part, seq_len(n_parts(null@region)))
  )
  reduce <- function(field, combine, template) {
    per <- vapply(per_cell, function(r) r[[field]], template)
    unname(vapply(by_part, function(i) combine(per[i]), template))
  }
  bounds <- reduce("bound", max, numeric(1L))
  incumbents <- pmax(reduce("incumbent", max, numeric(1L)), seeded)
  out <- list(
    certificate = ripr_certificate(
      sup_ub = max(bounds),
      sup_lb = max(incumbents),
      random_variable = x,
      null = null,
      method = methods,
      bounds = bounds,
      incumbents = incumbents,
      iterations = reduce("iterations", sum, integer(1L)),
      converged = reduce("converged", all, logical(1L)),
      budget_hit = reduce("budget_hit", any, logical(1L))
    )
  )
  if (record) {
    # Kept per cell: two cells' trees don't combine into one.
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
