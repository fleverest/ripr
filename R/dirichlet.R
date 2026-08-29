#' @include mixing_measure.R distribution.R multinomial.R region.R
#' @include polytope_region.R
NULL

# Dirichlet priors, both truncated and non-truncated versions
#
# The mixture a Dirichlet prior induces over multinomial counts is a ratio of
# two integrals of the same shape:
#
#   P_W(x) = C(x) I(alpha + x) / I(alpha),   I(b) = int_A prod_j theta_j^(b_j-1)
#
# over the prior's support `A`. The Beta normaliser and the truncation constant
# both sit in `I(alpha)` and cancel, so nothing here ever computes the mass the
# prior places on the region. A = full simplex gives `I = B`, the multivariate
# Beta function, and P_W is a Dirichlet-multinomial; a truncation replaces it
# with a quadrature rule.

# --- Concentrations -----------------------------------------------------------

#' Validate Dirichlet concentration parameters
#'
#' The closed form [dirichlet_mixing()] evaluates is exact for real ones, so
#' only the quadrature in [truncated_dirichlet()] needs the integer restriction.
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


#' A reference point for a Dirichlet distribution
#'
#' The mode when it is interior, which needs every concentration above 1, and
#' the mean otherwise. Only used to seed an optimiser.
#' @keywords internal
#' @noRd
dirichlet_centre <- function(alpha) {
  if (all(alpha > 1)) {
    (alpha - 1) / (sum(alpha) - length(alpha))
  } else {
    alpha / sum(alpha)
  }
}


#' `(K, n)` draws from `Dir(alpha)`
#' @keywords internal
#' @noRd
dirichlet_draws <- function(alpha, n) {
  k <- length(alpha)
  g <- matrix(stats::rgamma(n * k, shape = alpha), nrow = k, ncol = n)
  div_by_col(g, colSums(g))
}


# --- The shared parent --------------------------------------------------------

#' Dirichlet priors over the simplex
#'
#' The abstract parent of [dirichlet_mixing()] and [truncated_dirichlet()]:
#' a Dirichlet law \eqn{W = \mathrm{Dir}(\alpha)}{W = Dir(alpha)} over the
#' probability simplex, possibly truncated to a particular [region]. Paired with
#' a [multinomial_family()] it induces a continuous mixture.
#'
#' Concentrations must be positive. [truncated_dirichlet()] narrows that to
#' positive **integers**, because the exactness of Gauss-Jacobi quadrature
#' requires it; [dirichlet_mixing()] evaluates a closed form and takes any
#' positive reals.
#'
#' @param alpha Length-`K` vector of positive concentrations, `K >= 2`. The
#'   property both subclasses share; `dirichlet_type` is abstract and is
#'   not constructed directly. [truncated_dirichlet()] narrows this to whole
#'   numbers; [dirichlet_mixing()] does not.
#' @param region A [region] of the probability simplex in `R^K`, bounded and
#'   full-dimensional. [truncated_dirichlet()] only.
#' @param degree_slack Internal: Raise the rule's degree by this much. The
#'   default of `0` is already exact, so this is only useful for testing that it
#'   is working as expected. [truncated_dirichlet()] only.
#' @param max_nodes Refuse a rule with more nodes than this.
#'   [truncated_dirichlet()] only.
#' @return A `dirichlet_mixing` or a `truncated_dirichlet`; both are
#'   `dirichlet_type` objects.
#' @examples
#' # `dirichlet_type` is abstract; the two constructors subclass it:
#' S7::S7_inherits(dirichlet_mixing(alpha = c(4, 3, 2)), dirichlet_type)
#'
#' fam <- multinomial_family(n_trials = 6L, k = 3L)
#' Q <- fam(dirichlet_mixing(alpha = c(4, 3, 2)))
#' sum(exp(log_density(Q, enumerate_space(fam@sample_space))))
#'
#' # A uniform prior gives every outcome the same mass -- the Bose-Einstein
#' # count, `choose(n + K - 1, K - 1)` outcomes each of equal probability.
#' flat <- fam(dirichlet_mixing(alpha = c(1, 1, 1)))
#' unique(round(exp(log_density(flat, enumerate_space(fam@sample_space))), 12))
#' 1 / choose(6 + 3 - 1, 3 - 1)
#'
#' # The plurality null, and the region where candidate 1 wins outright.
#' plurality <- union(simplex_region(vertices = cbind(
#'   c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1)
#' )), simplex_region(vertices = cbind(
#'   c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)
#' )))
#' alt <- setdiff(fam@parameter_space, plurality)
#'
#' W1 <- truncated_dirichlet(alpha = c(4, 3, 2), region = alt)
#'
#' # Still a proper distribution, on the outcomes that region can produce.
#' sum(exp(log_density(fam(W1), enumerate_space(fam@sample_space))))
#' @export
dirichlet_type <- new_class(
  "dirichlet_type",
  parent = continuous_mixing,
  abstract = TRUE,
  properties = list(alpha = class_numeric)
)


#' Log of the unnormalised Dirichlet integral over a region
#'
#' \eqn{\log \int_A \prod_j \theta_j^{\beta_j - 1} \mathrm{d}\theta}{
#' log int_A prod_j theta_j^(beta_j - 1) dtheta}, evaluated at every column of
#' `shape` at once. This term is what differentiates the [dirichlet_mixing()]
#' from a [truncated_dirichlet()].
#'
#' No normalisation is applied here because call sites require a ratio of two
#' such integrals at a time, so any common factor cancels.
#'
#' The degree of the integrand is `colSums(shape) - K`, so `shape` has all the
#' information that a quadrature implementation needs.
#' @param mixing A [dirichlet_type].
#' @param shape `(K, M)` matrix of Dirichlet shape vectors, one per column:
#'   `alpha + x` for the numerator, `alpha` for the normaliser.
#' @return Length-`M` numeric vector.
#' @keywords internal
log_region_integral <- new_generic(
  "log_region_integral",
  "mixing",
  function(mixing, shape) S7::S7_dispatch()
)


#' @rdname dirichlet_type
#' @usage NULL
method(induced_log_density, list(dirichlet_type, multinomial_family)) <-
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
    x <- as_outcome_matrix(x)
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
    shape <- cbind(t(x) + mixing@alpha, mixing@alpha)
    log_i <- log_region_integral(mixing, shape)
    denom <- length(log_i)
    log_multinom_coef(x, family@n_trials) + log_i[-denom] - log_i[denom]
  }


# --- The untruncated case -----------------------------------------------------

#' @section Over the full probability simplex:
#' \eqn{W = \mathrm{Dir}(\alpha)}{W = Dir(alpha)} over the whole simplex. Paired
#' with a [multinomial_family()] it induces a Dirichlet-multinomial distribution
#' over count data.
#'
#' It has a closed form, so evaluating the density does not require quadrature
#' rules and is therefore exact.
#'
#' @rdname dirichlet_type
#' @order 2
#' @export
dirichlet_mixing <- new_class(
  "dirichlet_mixing",
  parent = dirichlet_type,
  constructor = function(alpha) {
    new_object(S7_object(), alpha = as_concentration(alpha))
  }
)


#' @description Over the whole simplex the integral is the multivariate Beta
#'   function, so this is exact and ignores the degree `shape` implies.
#' @rdname log_region_integral
#' @usage NULL
method(log_region_integral, dirichlet_mixing) <- function(mixing, shape) {
  colSums(lgamma(shape)) - lgamma(colSums(shape))
}


method(draw_theta, dirichlet_mixing) <- function(mixing, n) {
  dirichlet_draws(mixing@alpha, n)
}


#' @description The mode when it is interior, i.e. every concentration parameter
#'   is above 1, and the mean otherwise. Both only seed the optimiser's starting
#'   atom, so the fallback costs nothing but a slightly different starting
#'   point.
#' @rdname mode_parameter
#' @usage NULL
method(mode_parameter, dirichlet_mixing) <- function(x) {
  dirichlet_centre(x@alpha)
}


# --- Gauss-Jacobi quadrature on a simplex -------------------------------------
#
# The integrand is a monomial, so a rule exact to its degree gives the exact
# integral. The rule used is the collapsed-coordinate (Duffy) one: map the
# (K-1)-cube to the reference simplex by
#
#   lambda_1 = u_1,  lambda_2 = (1 - u_1) u_2,  ...,
#   lambda_K = prod_i (1 - u_i)
#
# whose Jacobian is prod_{i < K-1} (1 - u_i)^(K-1-i). Each factor is a Jacobi
# weight function, so a Gauss-Jacobi rule in direction `i` absorbs it exactly
# and the whole thing is a tensor product of one-dimensional rules. Gauss
# weights are positive, which log space requires, and Gauss nodes are strictly
# interior, which is what keeps `log(theta)` finite below.

#' Gauss-Jacobi nodes and weights on `[0, 1]` for the weight `(1-u)^a u^b`
#'
#' Golub-Welsch: the nodes are the eigenvalues of the symmetric tridiagonal
#' Jacobi matrix built from the three-term recurrence, and the weights come from
#' the first component of each orthonormalised eigenvector, scaled by the
#' weight function's total mass.
#'
#' Returned on `[0, 1]` rather than `[-1, 1]`, which is where the collapsed
#' coordinates live; the affine change of variable rescales the weights by
#' `2^-(a+b+1)`. Weights sum to `int_0^1 (1-u)^a u^b du`.
#' @param n Number of nodes.
#' @param a,b Jacobi exponents, both `>= 0`.
#' @return `list(nodes = , weights = )`, nodes strictly inside `(0, 1)`.
#' @keywords internal
#' @noRd
gauss_jacobi_01 <- function(n, a, b = 0) {
  ab <- a + b
  k <- seq_len(n) - 1L
  # Monic recurrence coefficients (Gautschi). The n = 0 diagonal entry is the
  # special case where the general formula divides by zero at a + b = 0.
  diagonal <- ifelse(
    k == 0L,
    (b - a) / (ab + 2),
    (b^2 - a^2) / ((2 * k + ab) * (2 * k + ab + 2))
  )
  # `mass` is beta_0, the integral of the weight function over [-1, 1].
  mass <- exp(
    (ab + 1) * log(2) + lgamma(a + 1) + lgamma(b + 1) - lgamma(ab + 2)
  )

  jacobi <- matrix(0, n, n)
  jacobi[cbind(seq_len(n), seq_len(n))] <- diagonal
  if (n > 1L) {
    j <- seq_len(n - 1L)
    beta_j <- ifelse(
      j == 1L,
      4 * (a + 1) * (b + 1) / ((ab + 2)^2 * (ab + 3)),
      4 *
        (j + a) *
        (j + b) *
        j *
        (j + ab) /
        ((2 * j + ab)^2 * (2 * j + ab + 1) * (2 * j + ab - 1))
    )
    off <- sqrt(beta_j)
    jacobi[cbind(j, j + 1L)] <- off
    jacobi[cbind(j + 1L, j)] <- off
  }
  ev <- eigen(jacobi, symmetric = TRUE)
  ord <- order(ev$values)
  list(
    nodes = (1 + ev$values[ord]) / 2,
    weights = mass * ev$vectors[1L, ord]^2 / 2^(ab + 1)
  )
}


#' Points per direction for a rule exact to a given polynomial degree
#'
#' An `q`-point Gauss rule is exact to degree `2q - 1`.
#' @keywords internal
#' @noRd
rule_points <- function(degree) max(1L, as.integer(ceiling((degree + 1) / 2)))


#' A degree-exact rule on the reference `(K-1)`-simplex
#'
#' Barycentric nodes and log weights for
#' \eqn{\int_T f \,\mathrm{d}\lambda}{int_T f dlambda} over the reference
#' simplex, exact for every polynomial of total degree at most `degree`.
#' @param k Number of barycentric coordinates.
#' @param degree Polynomial degree the rule must integrate exactly.
#' @return `list(lambda = (M, k) barycentric nodes, log_w = length-M)`.
#' @keywords internal
#' @noRd
reference_simplex_rule <- function(k, degree) {
  d <- k - 1L
  q <- rule_points(degree)
  # Direction `i` carries the Jacobian factor `(1 - u)^(d - i)`; the last one
  # carries none, so it is plain Gauss-Legendre.
  rules <- lapply(seq_len(d), function(i) gauss_jacobi_01(q, a = d - i))

  grid <- as.matrix(expand.grid(rep(list(seq_len(q)), d)))
  u <- matrix(0, nrow = nrow(grid), ncol = d)
  log_w <- numeric(nrow(grid))
  for (i in seq_len(d)) {
    u[, i] <- rules[[i]]$nodes[grid[, i]]
    log_w <- log_w + log(rules[[i]]$weights[grid[, i]])
  }

  lambda <- matrix(0, nrow = nrow(grid), ncol = k)
  remainder <- rep(1, nrow(grid))
  for (i in seq_len(d)) {
    lambda[, i] <- remainder * u[, i]
    remainder <- remainder * (1 - u[, i])
  }
  lambda[, k] <- remainder
  list(lambda = lambda, log_w = log_w)
}


#' Refuse a quadrature rule too large to build
#'
#' A resource limit rather than a correctness one, as [certify()]'s
#' `max_coefficients` is: raising it costs time and memory and nothing else.
#' The node count is `n_cells * q^(K-1)`, so it is the category count that
#' drives it.
#' @keywords internal
#' @noRd
check_quadrature_size <- function(
  n_cells,
  q,
  k,
  alpha,
  n_trials,
  degree,
  max_nodes
) {
  per_cell <- q^(k - 1L)
  total <- n_cells * per_cell
  if (total <= max_nodes) {
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
    count(max_nodes),
    ").\n",
    "The rule is exact to degree ",
    degree,
    ", for sum(alpha) = ",
    sum(alpha),
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
#' \eqn{W = \mathrm{Dir}(\alpha)}{W = Dir(alpha)} conditioned on
#' \eqn{\theta \in A}{theta in A} for a [region] `A` of the probability simplex.
#' This is how a mixing measure supported strictly on the alternative may be
#' written down: e.g. take the complement of the null region and truncate a
#' Dirichlet prior to it.
#'
#' @section The truncation region:
#' `region` is decomposed by `cells(disjoin(region))` into disjoint simplices,
#' which is what lets the integral be a sum over cells without double-counting
#' an overlap for arbitrary unions. Every cell must be full-dimensional within
#' the simplex: a cell with fewer than `K` vertices has Lebesgue measure zero,
#' so construction refuses a region with only lower-dimension cells. A region
#' with *some* of them keeps the rest and warns, since dropping a set of
#' measure zero does not change the resulting integral.
#'
#' @section Quadrature:
#' The integrand is a monomial of total degree `|alpha| + n - K`, the same for
#' every outcome since counts always sum to `n`. A collapsed-coordinate
#' Gauss-Jacobi rule of that degree on each cell therefore evaluates it exactly
#' and agrees with [dirichlet_mixing()] over the whole simplex to floating point
#' rounding.
#'
#' Non-integer concentrations are refused here, unlike in [dirichlet_mixing()]:
#' they make the integrand singular on the boundary faces, so the quadrature
#' rule would only be approximate.
#'
#' @rdname dirichlet_type
#' @order 3
#' @export
truncated_dirichlet <- new_class(
  "truncated_dirichlet",
  parent = dirichlet_type,
  properties = list(
    region = region,
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
    # Lower-dimensional cells integrate to zero, so they are dropped.
    # Only a region with nothing left is an error.
    full <- vapply(cells, is_full_cell, logical(1), k = k)
    if (!any(full)) {
      stop(
        "every cell spans fewer than ",
        k,
        " vertices, so the whole region has measure zero. A  truncated ",
        "Dirichlet needs a full-dimensional region.",
        call. = FALSE
      )
    }
    if (!all(full)) {
      warning(degenerate_cell_warning(sum(!full), length(full)))
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


#' Does a cell lie within the simplex?
#' @keywords internal
#' @noRd
check_simplex_cell <- function(cell, k, tol = 1e-9) {
  if (!S7_inherits(cell, polytope_region)) {
    stop(
      "every cell of `region` must be a bounded polytope; got a `",
      attr(S7_class(cell), "name"),
      "`.",
      call. = FALSE
    )
  }
  vertices <- cell@vertices
  # Check containment
  in_simplex <- all(vertices >= -tol) && all(abs(colSums(vertices) - 1) <= tol)
  if (!in_simplex) {
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
is_full_cell <- function(cell, k) ncol(cell@vertices) == k


#' The warning a dropped cell earns
#'
#' Classed, as `slice_warning()` is, so a caller who already knows their region
#' carries degenerate pieces can silence just this one.
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


#' The pooled quadrature rule over a truncated Dirichlet's cells
#'
#' Nodes are the reference simplex's, mapped into each cell by
#' `theta = V lambda`; the Jacobian is the constant `abs(det(V))` per cell,
#' which multiplies that cell's weights.
#'
#' The degree depends on the family's `n_trials`, which the measure meets only
#' through a density call. The rule is a deterministic function of
#' `(alpha, region, degree)`, so repeated calls share it and the
#' self-normalisation identity survives across them as well as within one.
#' @return `list(log_nodes = (M, k), log_omega = length-M)`.
#' @keywords internal
#' @noRd
truncated_rule <- function(mixing, degree) {
  k <- length(mixing@alpha)
  # `degree` is the exact one, `|alpha| + n - k`, so the trial count it came
  # from is recoverable for the error message below.
  n_trials <- degree + k - sum(mixing@alpha)
  degree <- max(0L, degree + mixing@degree_slack)
  q <- rule_points(degree)
  check_quadrature_size(
    length(mixing@cells),
    q,
    k,
    mixing@alpha,
    n_trials,
    degree,
    mixing@max_nodes
  )

  reference <- reference_simplex_rule(k, degree)
  nodes <- vector("list", length(mixing@cells))
  log_omega <- vector("list", length(mixing@cells))
  for (i in seq_along(mixing@cells)) {
    vertices <- mixing@cells[[i]]@vertices
    nodes[[i]] <- reference$lambda %*% t(vertices)
    log_omega[[i]] <- reference$log_w + log(abs(det(vertices)))
  }
  nodes <- do.call(rbind, nodes)

  # Guaranteed by the full-dimensionality check at construction rather than
  # assumed: a node on a face `theta_j = 0` would meet a zero exponent as
  # `0 * -Inf` below and give `NaN`.
  if (any(nodes <= 0)) {
    stop(
      "a quadrature node landed on the boundary of the simplex, where the ",
      "log density is undefined.",
      call. = FALSE
    )
  }
  list(log_nodes = log(nodes), log_omega = do.call(c, log_omega))
}


#' @description A quadrature sum over the cells, exact because the integrand is
#' a monomial of degree `max(colSums(shape)) - K`.
#' @rdname log_region_integral
#' @usage NULL
method(log_region_integral, truncated_dirichlet) <- function(mixing, shape) {
  rule <- truncated_rule(mixing, max(colSums(shape)) - length(mixing@alpha))
  # `truncated_rule()` has already established that no node coordinate is zero,
  # so no `0 * -Inf` arises.
  col_logsumexp(rule$log_nodes %*% (shape - 1) + rule$log_omega)
}


#' @description Rejection sampling from `Dir(alpha)`, keeping the proposals the
#'   region contains.
#' @rdname draw_theta
#' @usage NULL
method(draw_theta, truncated_dirichlet) <- function(mixing, n) {
  k <- length(mixing@alpha)
  out <- matrix(NA_real_, nrow = k, ncol = n)
  filled <- 0L
  proposed <- 0L
  accepted <- 0L
  # At the 1% floor below, `n` draws need about `100 n` proposals, so the cap
  # is only reached by a region that passed the rate check and was unlucky.
  cap <- 5000 + 500 * n

  while (filled < n) {
    batch <- max(256L, 2L * (n - filled))
    proposal <- dirichlet_draws(mixing@alpha, batch)
    keep <- apply(proposal, 2L, function(theta) contains(mixing@region, theta))
    proposed <- proposed + batch
    accepted <- accepted + sum(keep)

    take <- min(sum(keep), n - filled)
    if (take > 0L) {
      out[, filled + seq_len(take)] <- proposal[,
        which(keep)[seq_len(take)],
        drop = FALSE
      ]
      filled <- filled + take
    }
    rate <- accepted / proposed
    if (filled < n && (rate < 0.01 || proposed >= cap)) {
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


#' @description The untruncated mode or mean when the region contains it, and
#'   otherwise its projection onto whichever cell of the region is nearest.
#' @rdname mode_parameter
#' @usage NULL
method(mode_parameter, truncated_dirichlet) <- function(x) {
  centre <- dirichlet_centre(x@alpha)
  if (contains(x@region, centre)) {
    return(centre)
  }
  candidates <- lapply(x@cells, function(cell) project(cell, centre))
  distances <- vapply(candidates, function(p) sum((p - centre)^2), numeric(1))
  candidates[[which.min(distances)]]
}
