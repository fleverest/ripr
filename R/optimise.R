#' @include region.R
NULL

# --- Optimisation over a convex_region's chart ---------------------------------

#' Bundle an objective for [maximise_over()]
#'
#' @param value Function of a parameter vector returning a scalar.
#' @param grad Function of a parameter vector returning the gradient, length `d`.
#' @param value_batch Optional function of an `(N, d)` matrix, one parameter
#'   per row, returning `N` values. Used to score the seed grid; defaults to
#'   applying `value` per row.
#' @return A list to pass to [maximise_over()].
#' @keywords internal
objective <- function(value, grad, value_batch = NULL) {
  if (is.null(value_batch)) {
    value_batch <- function(theta_mat) {
      vapply(seq_len(nrow(theta_mat)), \(i) value(theta_mat[i, ]), numeric(1))
    }
  }
  list(value = value, grad = grad, value_batch = value_batch)
}


#' Cache only the last call of a one-argument function (by `identical()`)
#'
#' An optimiser asks for the value and then the gradient at the same point.
#' @keywords internal
#' @noRd
memoise_last <- function(f) {
  last_arg <- NULL
  last <- NULL
  function(a) {
    if (!is.null(last_arg) && identical(a, last_arg)) {
      return(last)
    }
    last <<- f(a)
    last_arg <<- a
    last
  }
}


#' Maximise an objective over a parameter space
#'
#' Multi-start SLSQP in the space's own [chart()]: seeds are scored in batch,
#' the best `n_restarts` are refined under the chart's constraints, and the
#' best refinement wins.
#'
#' **The result is a lower bound on the supremum, not the supremum**: the
#' objective is generally non-convex. Upper bounds need certification.
#'
#' `seeds` should include the current mixture's atoms; otherwise the value can
#' fall below `max_j G(theta_j)` (at least 1) and a duality gap computed from
#' it can come out spuriously negative.
#'
#' @param space A [convex_region].
#' @param obj An [objective()].
#' @param seeds Optional `(m, d)` matrix of parameter-space points to seed
#'   from, one per row.
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
  fn_gr <- memoise_last(function(u) {
    theta <- ch$to_theta(u)
    list(
      value = -obj$value(theta),
      gradient = -as.numeric(obj$grad(theta) %*% ch$jacobian(u))
    )
  })

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
    seeds <- as_row_matrix(seeds)
    coords <- lapply(seq_len(nrow(seeds)), \(i) ch$from_theta(seeds[i, ]))
    starts <- rbind(do.call(rbind, coords), starts)
  }

  scores <- obj$value_batch(ch$to_theta_batch(starts))
  top <- order(scores, decreasing = TRUE)[seq_len(min(
    n_restarts,
    length(scores)
  ))]

  best <- list(par = starts[top[1L], ], value = Inf)
  for (i in top) {
    fallback <- if (is.finite(scores[i])) -scores[i] else Inf
    res <- refine(starts[i, ], fallback = fallback)
    if (is.finite(res$value) && res$value < best$value) best <- res
  }
  # SLSQP constraints hold only to tolerance, so project back into the space
  # and re-evaluate: the value must be attained at a point of the space.
  theta <- project(space, ch$to_theta(best$par))
  value <- obj$value(theta)
  # Projection can leave the refined value below the best seed's; never
  # return worse than that.
  seed_theta <- project(space, ch$to_theta(starts[top[1L], ]))
  seed_value <- obj$value(seed_theta)
  if (seed_value > value) {
    theta <- seed_theta
    value <- seed_value
  }
  list(theta = theta, value = value)
}
