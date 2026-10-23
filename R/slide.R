# Sliding: descend the KL over every weight and atom position jointly
# (sfw_step()).

#' One sliding descent: every weight and every atom position at once
#'
#' The sliding step of \insertCite{Denoyelle2020}{ripr}:
#' minimise KL jointly over all weights and all atom positions, from the
#' current mixture, each atom in its own part's chart coordinates. SLSQP
#' handles the constraints directly: weights nonnegative and summing to one,
#' and each chart's bounds and equality. The gradient is exact: `-G` in the
#' weights and `-w_c E_Q[grad p_c / P]` in atom `c`, through the chart's
#' Jacobian.
#'
#' Returns the descended `atoms` and `weights`, never with higher KL than the
#' start: a failed or worsening descent returns the start unchanged.
#' @keywords internal
#' @noRd
slide_sweep <- function(state, ld, dg, max_eval) {
  engine <- state@engine
  atoms <- state@mixing@atoms
  w <- state@mixing@weights
  n <- nrow(atoms)
  prts <- parts(state@null@region)[state@part]
  charts <- lapply(prts, chart)
  n_par <- vapply(charts, \(ch) as.integer(ch$n_par), integer(1))
  idx <- split(
    n + seq_len(sum(n_par)),
    factor(rep(seq_len(n), n_par), levels = seq_len(n))
  )
  with_heq <- which(!vapply(charts, \(ch) is.null(ch$heq), logical(1)))
  # Gets the atom locations (in theta-space) as a matrix from the joint chart-
  # space representation
  theta_of <- function(x) {
    rows <- lapply(seq_len(n), \(c_i) charts[[c_i]]$to_theta(x[idx[[c_i]]]))
    matrix(unlist(rows), nrow = n, byrow = TRUE)
  }

  # SLSQP calls fn and gr separately at the same point, so cache the pair.
  fg <- memoise_last(function(x) {
    wx <- x[seq_len(n)]
    theta <- theta_of(x)
    ld_all <- ld(theta)
    log_p <- mixture_log_p(ld_all, wx)
    value <- kl_at(engine, log_p)
    if (!is.finite(value)) {
      return(list(value = Inf, gradient = numeric(length(x))))
    }
    g_theta <- -wx * dg(theta, engine@log_w - log_p)
    g_u <- lapply(seq_len(n), function(c_i) {
      as.vector(g_theta[c_i, ] %*% charts[[c_i]]$jacobian(x[idx[[c_i]]]))
    })
    list(
      value = value,
      gradient = c(-atom_g(ld_all, log_p, engine), unlist(g_u))
    )
  })
  heq <- function(x) {
    c(
      sum(x[seq_len(n)]) - 1,
      unlist(lapply(with_heq, \(c_i) charts[[c_i]]$heq(x[idx[[c_i]]])))
    )
  }
  heqjac <- function(x) {
    jac <- matrix(0, nrow = 1L + length(with_heq), ncol = length(x))
    jac[1L, seq_len(n)] <- 1
    for (r in seq_along(with_heq)) {
      c_i <- with_heq[[r]]
      jac[r + 1L, idx[[c_i]]] <- charts[[c_i]]$heqjac(x[idx[[c_i]]])
    }
    jac
  }

  # Joint chart-space representation is c(weights, chart_vals)
  x0 <- c(
    w,
    unlist(lapply(seq_len(n), \(c_i) charts[[c_i]]$from_theta(atoms[c_i, ])))
  )
  start <- fg(x0)$value
  unchanged <- list(atoms = atoms, weights = w)
  fit <- tryCatch(
    nloptr::slsqp(
      x0,
      fn = \(x) fg(x)$value,
      gr = \(x) fg(x)$gradient,
      lower = c(numeric(n), unlist(lapply(charts, `[[`, "lower"))),
      heq = heq,
      heqjac = heqjac,
      control = list(xtol_rel = 1e-10, maxeval = max_eval)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(unchanged)
  }

  # SLSQP meets its constraints only to tolerance: snap weights back onto the
  # simplex (exact zeros below rounding, so they drop) and atoms into their
  # parts.
  w_new <- fit$par[seq_len(n)]
  w_new[w_new < 1e-12 * max(w_new)] <- 0
  w_new <- normalise_weights(w_new)
  theta <- theta_of(fit$par)
  for (c_i in seq_len(n)) {
    theta[c_i, ] <- project(prts[[c_i]], theta[c_i, ])
  }
  kl_new <- kl_at(engine, mixture_log_p(ld(theta), w_new))
  if (!is.finite(kl_new) || kl_new > start) {
    return(unchanged)
  }
  list(atoms = theta, weights = w_new)
}


#' Slide a state: descend, merge coincident atoms, re-solve the weights
#'
#' The descent stops once KL stops changing in floating point, which can leave
#' the weights slightly off optimal for the atoms it found, so they are
#' re-solved with the atoms fixed by `correct_state()`, which first merges
#' atoms that slid onto the same point.
#' @keywords internal
#' @noRd
slide_state <- function(state, ld, dg, max_eval) {
  slid <- slide_sweep(state, ld, dg, max_eval = max_eval)
  correct_state(
    set_mixture(state, atoms = slid$atoms, weights = slid$weights),
    ld
  )$state
}
