# Merging coincident atoms.

#' Indicates whether two points agree in every coordinate to within
#' `atol + rtol * max(|a|, |b|)`?
#' @keywords internal
#' @noRd
coincident <- function(a, b, atol, rtol) {
  all(abs(a - b) <= atol + rtol * pmax(abs(a), abs(b)))
}


#' Fold atoms of the same part that slid onto the same point into one
#'
#' Two atoms of a part are the same point when every coordinate agrees to
#' within `atol + rtol * max(|a|, |b|)`. The merged atom sits at their
#' weight-averaged position (inside the part, which is convex) with their summed
#' weight, and keeps the earlier atom's part. With `by_part = FALSE` atoms of
#' different parts merge too, as points on a face two parts share; their
#' average may then lie just outside the kept part, so the caller projects it.
#' @keywords internal
#' @noRd
merge_coincident <- function(
  atoms,
  weights,
  part,
  atol = 1e-9,
  rtol = 1e-9,
  by_part = TRUE
) {
  # A merge moves the surviving atom, possibly onto one an earlier pass
  # compared it with, so passes repeat until one merges nothing.
  repeat {
    n <- nrow(atoms)
    alive <- rep(TRUE, n)
    for (a in seq_len(n - 1L)) {
      if (!alive[a]) {
        next
      }
      for (b in (a + 1L):n) {
        if (!alive[b] || (by_part && part[a] != part[b])) {
          next
        }
        if (coincident(atoms[a, ], atoms[b, ], atol, rtol)) {
          total <- weights[a] + weights[b]
          if (total > 0) {
            atoms[a, ] <- (weights[a] * atoms[a, ] + weights[b] * atoms[b, ]) /
              total
          }
          weights[a] <- total
          alive[b] <- FALSE
        }
      }
    }
    atoms <- atoms[alive, , drop = FALSE]
    weights <- weights[alive]
    part <- part[alive]
    if (all(alive)) {
      break
    }
  }
  list(atoms = atoms, weights = weights, part = part)
}


#' `merge_coincident()` on a state's atoms, at the control's tolerances
#'
#' A merged atom sits at the weight-averaged position of the atoms it replaces,
#' so the mixture changes only at second order in their (tiny) separation.
#' @keywords internal
#' @noRd
merge_state <- function(state) {
  ctl <- state@control
  m <- merge_coincident(
    state@mixing@atoms,
    state@mixing@weights,
    state@part,
    atol = ctl$merge_atol,
    rtol = ctl$merge_rtol
  )
  if (nrow(m$atoms) == nrow(state@mixing@atoms)) {
    return(state)
  }
  set_mixture(state, atoms = m$atoms, weights = m$weights, part = m$part)
}
