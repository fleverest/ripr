#' @include mixing.R distribution.R multinomial.R region.R
#' @include polytope_region.R gauss.R
NULL

# Dirichlet priors, untruncated and truncated. The induced mixture over counts
# is a ratio of integrals over the prior's support `A`:
#
#   P_W(x) = C(x) I(alpha + x) / I(alpha),   I(b) = int_A prod_j theta_j^(b_j-1)
#
# so the Beta normaliser and truncation constant cancel and the prior mass of
# `A` is never computed. On the full simplex `I` is the multivariate Beta
# function; a truncation replaces it with a quadrature rule.

# --- Concentrations -----------------------------------------------------------

#' Validate Dirichlet concentrations. Only the quadrature in
#' [truncated_dirichlet()] needs integers; [dirichlet()]'s closed form is exact.
#' @keywords internal
#' @noRd
as_concentration <- function(alpha, integer = FALSE) {
  alpha <- as.numeric(alpha)
  if (length(alpha) < 2L) {
    stop(
      "`alpha` must have at least 2 entries; a Dirichlet over one category ",
      "is the point mass at 1.",
      call. = FALSE
    )
  }
  if (anyNA(alpha) || any(!is.finite(alpha)) || any(alpha <= 0)) {
    stop(
      "every entry of `alpha` must be a finite positive number",
      call. = FALSE
    )
  }
  if (integer) {
    if (any(abs(alpha - round(alpha)) > 1e-9)) {
      stop(
        "every entry of `alpha` must be a positive integer; got (",
        paste(signif(alpha, 4L), collapse = ", "),
        ").\n",
        "With integer concentrations and integer counts the exponent vector ",
        "`alpha + x - 1` is non-negative, so the integrand is a polynomial and ",
        "a quadrature rule of matching degree integrates it exactly. A ",
        "non-integer concentration parameter makes the Dirichlet distribution ",
        "singular on the boundary faces of the simplex, and a quadrature ",
        "approach becomes only approximate, which is not yet supported.",
        call. = FALSE
      )
    }
    alpha <- as.integer(round(alpha))
  }
  alpha
}


#' The mode when interior (every concentration above 1), else the mean.
#' @keywords internal
#' @noRd
dirichlet_centre <- function(alpha) {
  if (all(alpha > 1)) {
    (alpha - 1) / (sum(alpha) - length(alpha))
  } else {
    alpha / sum(alpha)
  }
}


#' `n` draws from `Dir(alpha)` as an `(n, K)` matrix, by row. Each gamma
#' variate is drawn on the log scale as `log Gamma(alpha + 1) + log(U) / alpha`,
#' so small `alpha` cannot underflow to zero.
#' @keywords internal
#' @noRd
dirichlet_draws <- function(alpha, n) {
  k <- length(alpha)
  by_row <- \(x) matrix(x, nrow = n, ncol = k, byrow = TRUE)
  log_g <- by_row(log(stats::rgamma(n * k, shape = alpha + 1))) +
    by_row(log(stats::runif(n * k)) / alpha)
  exp(log_g - row_logsumexp(log_g))
}


# --- The shared parent --------------------------------------------------------

#' Dirichlet priors over the simplex
#'
#' A Dirichlet law \eqn{W = \mathrm{Dir}(\alpha)}{W = Dir(alpha)} over the
#' probability simplex, possibly truncated to a [region]. Paired with a
#' [multinomial_family()] it induces a continuous mixture. `dirichlet_dist` is
#' the abstract parent of [dirichlet()] and [truncated_dirichlet()].
#'
#' @param alpha Length-`K` vector of positive concentrations, `K >= 2`.
#'   [truncated_dirichlet()] requires whole numbers.
#' @param region A [region] of the probability simplex in `R^K`, bounded and
#'   full-dimensional. [truncated_dirichlet()] only.
#' @param degree_slack Raise the quadrature rule's degree by this much. The
#'   default of `0` is already exact; this exists for testing.
#'   [truncated_dirichlet()] only.
#' @param max_nodes Refuse a rule with more nodes than this.
#'   [truncated_dirichlet()] only.
#' @return A `dirichlet` or a `truncated_dirichlet`; both are
#'   `dirichlet_dist` objects.
#' @examples
#' # `dirichlet_dist` is abstract; the two constructors subclass it:
#' S7::S7_inherits(dirichlet(alpha = c(4, 3, 2)), dirichlet_dist)
#'
#' fam <- multinomial_family(n_trials = 6L, k = 3L)
#' Q <- fam(dirichlet(alpha = c(4, 3, 2)))
#' sum(exp(log_density(Q, enumerate_space(fam@sample_space))))
#'
#' # A uniform prior gives each of the `choose(n + K - 1, K - 1)` outcomes
#' # equal mass.
#' flat <- fam(dirichlet(alpha = c(1, 1, 1)))
#' unique(round(exp(log_density(flat, enumerate_space(fam@sample_space))), 12))
#' 1 / choose(6 + 3 - 1, 3 - 1)
#'
#' # The plurality null
#' plurality <- simplex_region(vertices = rbind(
#'   c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1)
#' )) | simplex_region(vertices = rbind(
#'   c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)
#' ))
#' # and the region where candidate 1 wins outright
#' alt <- fam@parameter_space - plurality
#'
#' W1 <- truncated_dirichlet(alpha = c(4, 3, 2), region = alt)
#'
#' # Still a proper distribution, on the outcomes that region can produce.
#' sum(exp(log_density(fam(W1), enumerate_space(fam@sample_space))))
#' @rdname dirichlet
#' @export
dirichlet_dist <- new_class(
  "dirichlet_dist",
  parent = continuous_dist,
  abstract = TRUE,
  properties = list(
    alpha = class_numeric,
    sample_space = new_property(
      space,
      getter = function(self) {
        simplex_region(vertices = diag(length(self@alpha)))
      }
    )
  )
)


#' @rdname dirichlet
#' @usage NULL
method(format, dirichlet_dist) <- function(x, ...) {
  sprintf(
    "%s: alpha %s",
    class_name(x),
    theta_label(x@alpha)
  )
}


#' Log of the unnormalised Dirichlet integral over a region
#'
#' \eqn{\log \int_A \prod_j \theta_j^{\beta_j - 1} \mathrm{d}\theta}{
#' log int_A prod_j theta_j^(beta_j - 1) dtheta} for every row of `shape`.
#' Unnormalised because callers take a ratio of two, so common factors cancel.
#' @param mixing A [dirichlet_dist].
#' @param shape `(M, K)` matrix of Dirichlet shape vectors, one per row:
#'   `alpha + x` for the numerator, `alpha` for the normaliser.
#' @return Length-`M` numeric vector.
#' @keywords internal
log_region_integral <- new_generic(
  "log_region_integral",
  "mixing",
  function(mixing, shape) S7::S7_dispatch()
)


#' @rdname dirichlet
#' @usage NULL
method(mixture_log_density, list(dirichlet_dist, multinomial_family)) <-
  function(mixing, family, x) {
    k <- length(mixing@alpha)
    if (k != family@k) {
      stop(
        "`alpha` has ",
        k,
        " entries but the family has ",
        family@k,
        " categories; a Dirichlet prior needs one concentration per category.",
        call. = FALSE
      )
    }
    x <- as_row_matrix(x)
    if (ncol(x) != k) {
      stop(
        "outcomes must have ",
        k,
        " columns, one per category; got ",
        ncol(x),
        ".",
        call. = FALSE
      )
    }
    shape <- rbind(add_by_col(x, mixing@alpha), mixing@alpha)
    log_i <- log_region_integral(mixing, shape)
    denom <- length(log_i)
    log_multinom_coef(x, family@n_trials) + log_i[-denom] - log_i[denom]
  }


# --- The untruncated case -----------------------------------------------------

#' @section Over the full probability simplex:
#' [dirichlet()] is \eqn{\mathrm{Dir}(\alpha)}{Dir(alpha)} over the whole
#' simplex. With a [multinomial_family()] it induces the Dirichlet-multinomial,
#' evaluated exactly in closed form; `alpha` may be any positive reals.
#'
#' @rdname dirichlet
#' @order 2
#' @export
dirichlet <- new_class(
  "dirichlet",
  parent = dirichlet_dist,
  constructor = function(alpha) {
    new_object(S7_object(), alpha = as_concentration(alpha))
  }
)


#' @description Over the whole simplex the integral is the multivariate Beta
#'   function, so this is exact and ignores the degree `shape` implies.
#' @rdname log_region_integral
#' @usage NULL
method(log_region_integral, dirichlet) <- function(mixing, shape) {
  rowSums(lgamma(shape)) - lgamma(rowSums(shape))
}


#' @rdname draw
#' @usage NULL
method(draw, dirichlet) <- function(dist, n) {
  dirichlet_draws(dist@alpha, n)
}


#' @noRd
method(reference_point, dirichlet) <- function(x) {
  dirichlet_centre(x@alpha)
}


# --- Gauss-Jacobi quadrature on a simplex -------------------------------------
#
# The integrand is a monomial, so a rule exact to its degree is exact. The
# collapsed-coordinate (Duffy) map from the (K-1)-cube,
#
#   lambda_1 = u_1,  lambda_2 = (1 - u_1) u_2,  ...,  lambda_K = prod_i (1 - u_i)
#
# has Jacobian prod_{i < K-1} (1 - u_i)^(K-1-i); each factor is a Jacobi weight,
# so the rule is a tensor product of 1-D Gauss-Jacobi rules. Gauss weights are
# positive (log space needs that) and nodes strictly interior (keeping
# `log(theta)` finite).

#' Points per direction for a rule exact to `degree` (`q` points: `2q - 1`).
#' @keywords internal
#' @noRd
rule_points <- function(degree) max(1L, as.integer(ceiling((degree + 1) / 2)))


#' A rule on the reference `(k-1)`-simplex exact to total degree `degree`.
#' @return `list(lambda = (M, k) barycentric nodes, log_w = length-M)`.
#' @keywords internal
#' @noRd
reference_simplex_rule <- function(k, degree) {
  d <- k - 1L
  q <- rule_points(degree)
  # Direction `i` carries the Jacobian factor `(1 - u)^(d - i)`; the last one
  # carries none, so it is plain Gauss-Legendre.
  grid <- tensor_rule(
    lapply(seq_len(d), function(i) gauss_jacobi_01(q, a = d - i))
  )

  u <- grid$nodes
  lambda <- matrix(0, nrow = nrow(u), ncol = k)
  remainder <- rep(1, nrow(u))
  for (i in seq_len(d)) {
    lambda[, i] <- remainder * u[, i]
    remainder <- remainder * (1 - u[, i])
  }
  lambda[, k] <- remainder
  list(lambda = lambda, log_w = grid$log_w)
}


#' Refuse a rule with more than `max_nodes` nodes (`n_cells * q^(K-1)`).
#' `n_trials` is only reported.
#' @keywords internal
#' @noRd
check_quadrature_size <- function(mixing, degree, n_trials) {
  k <- length(mixing@alpha)
  q <- rule_points(degree)
  n_cells <- length(mixing@cells)
  total <- n_cells * q^(k - 1L)
  if (total <= mixing@max_nodes) {
    return(invisible(total))
  }
  count <- function(x) format(x, big.mark = ",", scientific = FALSE)
  stop(
    "the quadrature rule would need ",
    count(total),
    " nodes (",
    n_cells,
    " cells x ",
    q,
    "^",
    k - 1L,
    " per cell), above `max_nodes` (",
    count(mixing@max_nodes),
    ").\n",
    "The rule is exact to degree ",
    degree,
    ", for sum(alpha) = ",
    sum(mixing@alpha),
    ", n_trials = ",
    n_trials,
    " and k = ",
    k,
    ", which needs ",
    q,
    " points in each of the k - 1 = ",
    k - 1L,
    " directions.\n",
    "Try with fewer trials or fewer categories, or raise `max_nodes`.",
    call. = FALSE
  )
}


# --- The truncated case -------------------------------------------------------

#' @section Truncated to a region of the simplex:
#' [truncated_dirichlet()] is \eqn{\mathrm{Dir}(\alpha)}{Dir(alpha)} conditioned
#' on \eqn{\theta \in A}{theta in A} for a [region] `A` of the simplex, e.g. the
#' complement of the null, giving a mixing measure supported on the alternative.
#'
#' `region` is split into disjoint simplices by `cells(disjoin(region))`. Cells
#' with fewer than `K` vertices have measure zero and are dropped with a
#' warning; a region with no full-dimensional cell is refused.
#'
#' The integrand is a monomial of degree `|alpha| + n - K`, so a Gauss-Jacobi
#' rule of that degree on each cell is exact. This needs integer `alpha`:
#' non-integer concentrations make the integrand singular on the boundary.
#'
#' @rdname dirichlet
#' @order 3
#' @export
truncated_dirichlet <- new_class(
  "truncated_dirichlet",
  parent = dirichlet_dist,
  properties = list(
    region = region,
    sample_space = new_property(space, getter = function(self) self@region),
    cells = class_list,
    degree_slack = class_numeric,
    max_nodes = class_numeric
  ),
  constructor = function(alpha, region, degree_slack = 0L, max_nodes = 1e6) {
    alpha <- as_concentration(alpha, integer = TRUE)
    k <- length(alpha)
    region <- as_region(region)
    rlang::check_number_whole(degree_slack)
    rlang::check_number_decimal(max_nodes, min = 1)

    if (is_empty(region)) {
      stop(
        "`region` is empty, so there is nothing to truncate to.",
        call. = FALSE
      )
    }
    if (space_dim(region) != k) {
      stop(
        "`region` lives in ",
        space_dim(region),
        " dimensions but `alpha` has ",
        k,
        " entries; a Dirichlet over `K` categories truncates to a region of ",
        "R^K.",
        call. = FALSE
      )
    }
    if (!is_bounded(region)) {
      stop(
        "`region` is unbounded, so it is not a subset of the probability ",
        "simplex. Intersect it with the family's `parameter_space` first.",
        call. = FALSE
      )
    }
    cells <- cells(disjoin(region))
    for (cell in cells) {
      check_simplex_cell(cell, k)
    }
    full <- vapply(cells, is_full_cell, logical(1), k = k)
    if (!any(full)) {
      stop(
        "every cell spans fewer than ",
        k,
        " vertices, so the whole region has measure zero. A truncated ",
        "Dirichlet needs a full-dimensional region.",
        call. = FALSE
      )
    }
    if (!all(full)) {
      warning(degenerate_cell_warning(sum(!full), length(full)))
      # Reset the support so containment checks ignore a dropped sliver, which
      # may lie a rounding error outside the simplex.
      region <- union_region(cells[full])
    }
    cells <- cells[full]
    new_object(
      S7_object(),
      alpha = alpha,
      region = region,
      cells = cells,
      degree_slack = as.integer(degree_slack),
      max_nodes = max_nodes
    )
  }
)


#' @description `print()` on a truncated Dirichlet shows the concentration and
#'   the region it is truncated to.
#' @rdname dirichlet
#' @usage NULL
method(print, truncated_dirichlet) <- function(x, ...) {
  cat("<truncated_dirichlet>\n")
  cat("  alpha  ", theta_label(x@alpha), "\n", sep = "")
  cat("  region ", format(x@region), "\n", sep = "")
  invisible(x)
}


#' Stop unless a cell is a polytope within the simplex.
#' @keywords internal
#' @noRd
check_simplex_cell <- function(cell, k, tol = 1e-9) {
  if (!S7_inherits(cell, polytope_region)) {
    stop(
      "every cell of `region` must be a bounded polytope; got a `",
      class_name(cell),
      "`.",
      call. = FALSE
    )
  }
  outside <- simplex_departure(cell@vertices, neg_tol = tol, sum_tol = tol)
  if (!is.null(outside)) {
    stop(
      "`region` is not a subset of the probability simplex: it has a vertex ",
      "with a negative coordinate or with coordinates not summing to 1. A ",
      "Dirichlet places no mass outside the simplex.",
      call. = FALSE
    )
  }
  invisible(NULL)
}


#' Is a cell full-dimensional within the simplex?
#' @keywords internal
#' @noRd
is_full_cell <- function(cell, k) nrow(cell@vertices) == k


#' Classed so a caller expecting degenerate cells can silence just this one.
#' @keywords internal
#' @noRd
degenerate_cell_warning <- function(n_dropped, n_total) {
  structure(
    class = c("ripr_degenerate_warning", "warning", "condition"),
    list(
      message = paste0(
        n_dropped,
        " of ",
        n_total,
        " cells that `region` decomposes into span too few vertices, so they ",
        "have Dirichlet-measure zero in the simplex. They have been dropped. ",
        "Only ",
        n_total - n_dropped,
        " cells remain."
      ),
      call = NULL
    )
  )
}


#' The pooled rule over a truncated Dirichlet's cells: reference nodes mapped by
#' `theta = lambda V`, weights scaled by `abs(det(V))` per cell.
#' @return `list(log_nodes = (M, k), log_omega = length-M)`.
#' @keywords internal
#' @noRd
truncated_rule <- function(mixing, degree) {
  k <- length(mixing@alpha)
  # Recover `n` from `degree = |alpha| + n - k`, for the error message only.
  n_trials <- degree + k - sum(mixing@alpha)
  degree <- max(0L, degree + mixing@degree_slack)
  check_quadrature_size(mixing, degree, n_trials)

  reference <- reference_simplex_rule(k, degree)
  nodes <- vector("list", length(mixing@cells))
  log_omega <- vector("list", length(mixing@cells))
  for (i in seq_along(mixing@cells)) {
    vertices <- mixing@cells[[i]]@vertices
    nodes[[i]] <- reference$lambda %*% vertices
    log_omega[[i]] <- reference$log_w + log(abs(det(vertices)))
  }
  nodes <- do.call(rbind, nodes)

  # A node on a face `theta_j = 0` would give `0 * -Inf = NaN` downstream.
  if (any(nodes <= 0)) {
    stop(
      "a quadrature node landed on the boundary of the simplex, where the ",
      "log density is undefined.",
      call. = FALSE
    )
  }
  list(log_nodes = log(nodes), log_omega = do.call(c, log_omega))
}


#' @description A quadrature sum over the cells, exact to degree
#'   `max(rowSums(shape)) - K`.
#' @rdname log_region_integral
#' @usage NULL
method(log_region_integral, truncated_dirichlet) <- function(mixing, shape) {
  rule <- truncated_rule(mixing, max(rowSums(shape)) - length(mixing@alpha))
  col_logsumexp(tcrossprod(rule$log_nodes, shape - 1) + rule$log_omega)
}


#' @description Rejection sampling from `Dir(alpha)`, keeping the proposals the
#'   region contains.
#' @rdname draw
#' @usage NULL
method(draw, truncated_dirichlet) <- function(dist, n) {
  mixing <- dist
  k <- length(mixing@alpha)
  out <- matrix(NA_real_, nrow = n, ncol = k)
  filled <- 0L
  proposed <- 0L
  accepted <- 0L
  cap <- 5000 + 500 * n

  while (filled < n) {
    batch <- max(256L, 2L * (n - filled))
    proposal <- dirichlet_draws(mixing@alpha, batch)
    keep <- apply(proposal, 1L, function(theta) contains(mixing@region, theta))
    proposed <- proposed + batch
    accepted <- accepted + sum(keep)

    take <- min(sum(keep), n - filled)
    if (take > 0L) {
      out[filled + seq_len(take), ] <- proposal[
        which(keep)[seq_len(take)],
        ,
        drop = FALSE
      ]
      filled <- filled + take
    }
    rate <- accepted / proposed
    too_rare <- proposed >= 2000L && rate < 0.01
    if (filled < n && (too_rare || proposed >= cap)) {
      stop(
        "rejection sampling from this truncated Dirichlet accepted ",
        accepted,
        " of ",
        proposed,
        " proposals (",
        signif(100 * rate, 3),
        "%), too few to draw ",
        n,
        " in reasonable time. The region holds very little of `Dir(alpha)`; ",
        "concentrate `alpha` towards it, or enlarge the region.",
        call. = FALSE
      )
    }
  }
  out
}


#' @description The untruncated centre if the region contains it, else its
#'   projection onto the nearest cell.
#' @noRd
method(reference_point, truncated_dirichlet) <- function(x) {
  centre <- dirichlet_centre(x@alpha)
  if (contains(x@region, centre)) {
    return(centre)
  }
  candidates <- lapply(x@cells, function(cell) project(cell, centre))
  distances <- vapply(candidates, function(p) sum((p - centre)^2), numeric(1))
  candidates[[which.min(distances)]]
}
