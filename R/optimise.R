#' @include region.R
NULL

# --- Optimisation over a convex_region's chart ---------------------------------

#' Bundle an objective for [maximise_over()]
#'
#' @param value Function of a parameter vector returning a scalar.
#' @param grad Function of a parameter vector returning the gradient, length `d`.
#' @param value_batch Optional function of a `(d, N)` matrix returning `N`
#'   values. Used to score the seed grid; defaults to applying `value` per
#'   column, which is correct but slower.
#' @return A list to pass to [maximise_over()].
#' @keywords internal
objective <- function(value, grad, value_batch = NULL) {
  if (is.null(value_batch)) {
    value_batch <- function(theta_mat) {
      vapply(seq_len(ncol(theta_mat)), \(i) value(theta_mat[, i]), numeric(1))
    }
  }
  list(value = value, grad = grad, value_batch = value_batch)
}


#' Maximise an objective over a parameter space
#'
#' Multi-start SLSQP in the space's own [chart()]: seed coordinates are scored
#' in batch, the best `n_restarts` are refined under the chart's declared
#' constraints, and the best refinement wins. Written once and shared by all
#' geometries, since only the charts should differ.
#'
#' **The result is a lower bound on the true supremum, not the supremum.** The
#' objective is generally non-convex, so this search is a heuristic. Restarts
#' will converge to local maxima. Anything needing an upper bound must obtain
#' it some other way, e.g. via certification.
#'
#' `seeds` should always include the current mixture's atoms. Without them the
#' returned value can fall below `max_j G(theta_j)`, which the mixture already
#' guarantees is at least 1, and a duality gap computed from it would come out
#' spuriously negative.
#'
#' @param space A [convex_region].
#' @param obj An [objective()].
#' @param seeds Optional `(d, m)` matrix of parameter-space points to seed from.
#' @param n_seeds Random seeds drawn from the chart.
#' @param n_restarts How many of the best seeds to refine.
#' @return `list(theta = , value = )` with `theta` in the space.
#' @keywords internal
maximise_over <- function(
  space,
  obj,
  seeds = NULL,
  n_seeds = 200L,
  n_restarts = 25L
) {
  ch <- chart(space)

  # A zero-dimensional chart has nothing to search: the space is a point.
  if (ch$n_par == 0L) {
    theta <- ch$to_theta(numeric(0))
    return(list(theta = theta, value = obj$value(theta)))
  }

  # slsqp() calls fn and gr separately at the same point, so cache the pair.
  last_u <- NULL
  last <- NULL
  fn_gr <- function(u) {
    if (!is.null(last_u) && identical(u, last_u)) {
      return(last)
    }
    theta <- ch$to_theta(u)
    last_u <<- u
    last <<- list(
      value = -obj$value(theta),
      gradient = -as.numeric(obj$grad(theta) %*% ch$jacobian(u))
    )
    last
  }

  refine <- function(u0, fallback) {
    res <- tryCatch(
      {
        fit <- nloptr::slsqp(
          u0,
          fn = \(u) fn_gr(u)$value,
          gr = \(u) fn_gr(u)$gradient,
          lower = ch$lower,
          heq = ch$heq,
          heqjac = ch$heqjac,
          control = list(xtol_rel = 1e-8, maxeval = 1000L)
        )
        list(par = fit$par, value = fit$value)
      },
      error = function(e) list(par = u0, value = fallback)
    )
    if (!is.finite(res$value)) {
      res <- list(par = u0, value = fallback)
    }
    res
  }

  starts <- ch$seed(n_seeds)
  if (!is.null(seeds)) {
    seeds <- as.matrix(seeds)
    coords <- vapply(
      seq_len(ncol(seeds)),
      \(i) ch$from_theta(seeds[, i]),
      numeric(ch$n_par)
    )
    starts <- cbind(matrix(coords, nrow = ch$n_par), starts)
  }

  scores <- obj$value_batch(ch$to_theta_batch(starts))
  top <- order(scores, decreasing = TRUE)[seq_len(min(
    n_restarts,
    length(scores)
  ))]

  best <- list(par = starts[, top[1L]], value = Inf)
  for (i in top) {
    fallback <- if (is.finite(scores[i])) -scores[i] else Inf
    res <- refine(starts[, i], fallback = fallback)
    if (is.finite(res$value) && res$value < best$value) best <- res
  }
  # SLSQP constraints hold only to its own tolerance, so the final iterate can
  # lie just outside the space, and its objective value marginally above the
  # supremum over the space. Here we just project the result back in and
  # re-evaluate, so the value is attained at a point of the space.
  theta <- project(space, ch$to_theta(best$par))
  value <- obj$value(theta)
  # The search shouldn't return something worse than its best seed.
  # SLSQP could stop a tiny step away from where it started and get something
  # worse because the projection re-evaluates. Without this the refined value
  # could be below one already attained by the seed.
  seed_theta <- project(space, ch$to_theta(starts[, top[1L]]))
  seed_value <- obj$value(seed_theta)
  if (seed_value > value) {
    theta <- seed_theta
    value <- seed_value
  }
  list(theta = theta, value = value)
}
