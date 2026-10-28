# Oracles and step rules
#
# Oracles propose where a new atom should go. Step rules decide whether or not
# to add it and how to update the weights accordingly. EM sweeps also move
# atoms but do not change the support size by adding or removing atoms.
#
# Oracles and steps work on `(ld_all, w)`: the log density of every atom at the
# quadrature nodes, and the weights. Atoms do not move, so the caller builds
# `ld_all`; `log_p` is rebuilt after each weight update. The candidate atom is
# always the last column of `ld_all`, entering with weight zero, which is where
# `add_atom()` puts it.

# --- Step verbs ---------------------------------------------------------------

#' Search every cell; return the best `theta`, `value` and `part`, and with
#' `keep_maxima = TRUE` every cell's refined local maxima (`maxima`: `theta`,
#' `value`, `part`)
#'
#' `seeds` defaults to the current atoms: `sum_c w_c G(theta_c) = 1` forces
#' `max_c G(theta_c) >= 1`, so including them keeps the maximum at least one.
#' @keywords internal
#' @noRd
search_null <- function(
  state,
  obj,
  seeds = state@mixing@atoms,
  keep_maxima = FALSE
) {
  ctl <- state@control
  found <- lapply(
    state@null@cells,
    \(s) {
      maximise_over(
        s,
        obj,
        seeds = seeds,
        n_seeds = ctl$n_seeds,
        n_restarts = ctl$n_restarts,
        xtol = ctl$search_xtol,
        max_eval = ctl$search_max_eval,
        keep_maxima = keep_maxima
      )
    }
  )
  best <- which.max(vapply(found, \(f) f$value, numeric(1)))
  out <- list(
    theta = found[[best]]$theta,
    value = found[[best]]$value,
    part = state@null@cell_part[[best]]
  )
  if (!keep_maxima) {
    return(out)
  }
  out$maxima <- list(
    theta = do.call(rbind, lapply(found, \(f) f$maxima$theta)),
    value = unlist(lapply(found, \(f) f$maxima$value)),
    part = unlist(lapply(seq_along(found), function(i) {
      rep(state@null@cell_part[[i]], length(found[[i]]$maxima$value))
    }))
  )
  out
}


#' Write a planned step back to the state
#'
#' The candidate joins only if the step weighted it. Incumbents driven towards
#' zero stay, since only an oracle can grow the support back.
#' @keywords internal
#' @noRd
commit_step <- function(state, theta, part, planned) {
  stepped <- if (planned$uses_candidate) {
    merge_state(add_atom(state, theta, part, planned$weights))
  } else {
    set_mixture(state, weights = planned$weights[-length(planned$weights)])
  }
  drop_empty(stepped)
}


# --- Multi-point insertion ----------------------------------------------------

#' The oracle's atom, then every other distinct local maximiser of `G` above one
#'
#' The multi-point insertion of the primal-dual active point (PDAP) method
#' \insertCite{PieperWalter2021}{ripr}. The inserted points come from the
#' refined multi-start search, `n_restarts` starts per cell, and
#' `found$maxima` keeps their end points. Those with `G > 1` are candidates,
#' best first. Candidates within `tol` (`near_tol` in [ripr_control()]; in
#' every coordinate) of each other are the same local maximum reached twice,
#' and merge into the better one, in any part: a maximum on a face two parts
#' share is found from both. A candidate within `tol` of a current atom, in
#' any part, is dropped. Returns a list of `list(theta, part)`, the oracle's
#' atom first.
#' @keywords internal
#' @noRd
insert_candidates <- function(
  state,
  found,
  tol = state@control$near_tol
) {
  oracle <- list(theta = as.vector(found$theta), part = as.integer(found$part))
  maxima <- found$maxima
  if (is.null(maxima)) {
    return(list(oracle))
  }
  up <- which(is.finite(maxima$value) & maxima$value > 1)
  up <- up[order(maxima$value[up], decreasing = TRUE)]
  # Across parts: a maximum on a face two parts share is found from both.
  merged <- merge_coincident(
    rbind(oracle$theta, maxima$theta[up, , drop = FALSE], deparse.level = 0),
    c(found$value, maxima$value[up]),
    c(oracle$part, as.integer(maxima$part[up])),
    atol = tol,
    rtol = tol,
    by_part = FALSE
  )
  prts <- parts(state@null@region)
  candidate <- function(i) {
    p <- merged$part[[i]]
    list(theta = as.vector(project(prts[[p]], merged$atoms[i, ])), part = p)
  }
  atoms <- state@mixing@atoms
  kept <- list(candidate(1L))
  for (i in seq_len(nrow(merged$atoms))[-1L]) {
    theta <- merged$atoms[i, ]
    same_as_atom <- vapply(
      seq_len(nrow(atoms)),
      \(c_i) coincident(atoms[c_i, ], theta, tol, tol),
      logical(1)
    )
    if (!any(same_as_atom)) {
      kept[[length(kept) + 1L]] <- candidate(i)
    }
  }
  kept
}


#' Insert every candidate at weight zero, then solve the weights
#'
#' `correct_state()` over the current atoms and every candidate gives the
#' candidates their mass: each has `G > 1`, so the solve's first QP frees it.
#' No Frank--Wolfe step seeds the solve, since its result would be replaced.
#' Atoms left at zero are dropped. Returns `list(state, log_p, row)`, as
#' `step_towards()` does; the row's `part` is the oracle's atom's, if it kept
#' any weight.
#' @keywords internal
#' @noRd
step_towards_all <- function(state, candidates, ld) {
  found <- candidates[[1L]]
  inserted <- set_mixture(
    state,
    atoms = rbind(
      state@mixing@atoms,
      do.call(rbind, lapply(candidates, `[[`, "theta")),
      deparse.level = 0
    ),
    weights = c(state@mixing@weights, numeric(length(candidates))),
    part = c(state@part, vapply(candidates, `[[`, integer(1), "part"))
  )
  corrected <- correct_state(inserted, ld)
  # The oracle's atom may have merged with a coincident one.
  ctl <- state@control
  atoms <- corrected$state@mixing@atoms
  oracle_kept <- any(apply(
    atoms,
    1L,
    coincident,
    found$theta,
    ctl$merge_atol,
    ctl$merge_rtol
  ))
  list(
    state = corrected$state,
    log_p = corrected$log_p,
    row = list(
      kl = kl_at(state@engine, corrected$log_p),
      solve_kkt = corrected$kkt,
      solve_converged = corrected$converged,
      oracle_value = NA_real_,
      oracle_theta = found$theta,
      part = if (oracle_kept) found$part else NA_integer_,
      step_size = NA_real_,
      direction = NA_character_
    )
  )
}
