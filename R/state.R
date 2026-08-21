#' @include null.R quadrature.R mixture.R
NULL

# The state representation of the ripr optimiser: the class, the conversions between
# the per-subnull atom lists and the flat vectors the arithmetic works on, and the
# trace. `steps.R` and `scheduler.R` both depend on this.

#' State of a RIPr fit
#'
#' Atoms and weights are lists indexed by subnull: element `i` holds the atoms
#' belonging to \eqn{\Theta_{0i}}{Theta_0i}. The subnull tag matters because the
#' pieces overlap, so an atom in an intersection would otherwise be ambiguous
#' about which chart to search in, and a parallel index vector would make that
#' an invariant to check rather than one to hold structurally. Weights are
#' normalised across the whole list, not within a subnull.
#'
#' @param atoms List of `(d, n_i)` matrices, one per subnull.
#' @param weights List of numeric vectors matching `atoms`, summing to 1 overall.
#' @param alternative The alternative \eqn{Q}{Q}.
#' @param null A [null_model].
#' @param engine A resolved [quadrature].
#' @param control From [ripr_control()].
#' @param trace Data frame, one row per recorded event.
#' @param snapshots List of recorded mixtures.
#' @param iters Named integer counts of steps taken, one name per verb.
#' @return A `ripr_state`.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_null(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_null(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' Q <- mixture(point_mixing(c(0.4, 0.35, 0.25)), fam)
#' state <- ripr_init(Q, plurality)
#' state
#' @export
ripr_state <- new_class(
  "ripr_state",
  properties = list(
    atoms = class_list,
    weights = class_list,
    alternative = outcome_distribution,
    null = null_model,
    engine = quadrature,
    control = class_list,
    trace = class_any,
    snapshots = class_list,
    iters = class_integer
  ),
  validator = function(self) {
    if (length(self@atoms) != length(self@weights)) {
      return("`weights` must have one element per element of `atoms`")
    }
    if (length(self@atoms) != length(self@null@subnulls)) {
      return("`atoms` must have one element per subnull")
    }
    sizes <- vapply(self@atoms, ncol, integer(1))
    if (!identical(sizes, vapply(self@weights, length, integer(1)))) {
      return("each weight vector must match the column count of its atoms")
    }
    total <- sum(unlist(self@weights))
    if (sum(sizes) > 0L && abs(total - 1) > 1e-8) {
      return("weights must sum to 1 across all subnulls")
    }
    if (is.null(names(self@iters))) {
      return("`iters` must be a named integer vector, one name per verb")
    }
    NULL
  }
)

# --- Flat views ---------------------------------------------------------------

#' All atoms as one `(d, C)` matrix
#'
#' Columns are grouped by subnull, in subnull order, and a new atom is appended
#' within its own subnull's block. So the flat ordering is not chronological: an
#' atom added at iteration 5 may sit before one added at iteration 2. Use
#' [flat_subnull()] to recover which block a column came from.
#' @param state A [ripr_state].
#' @return `(d, C)` numeric matrix.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_null(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_null(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' Q <- mixture(point_mixing(c(0.4, 0.35, 0.25)), fam)
#' state <- ripr_init(Q, plurality)
#' flat_atoms(state)
#' @export
flat_atoms <- function(state) {
  keep <- block_sizes(state) > 0L
  if (!any(keep)) {
    return(matrix(numeric(0), nrow = 0L, ncol = 0L))
  }
  do.call(cbind, state@atoms[keep])
}

#' All weights as one vector, aligned with [flat_atoms()]
#' @param state A [ripr_state].
#' @return Numeric vector summing to 1.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_null(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_null(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' Q <- mixture(point_mixing(c(0.4, 0.35, 0.25)), fam)
#' state <- ripr_init(Q, plurality)
#' flat_weights(state)
#' @export
flat_weights <- function(state) unlist(state@weights, use.names = FALSE)

#' Which subnull each column of [flat_atoms()] belongs to
#' @param state A [ripr_state].
#' @return Integer vector.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_null(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_null(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' Q <- mixture(point_mixing(c(0.4, 0.35, 0.25)), fam)
#' state <- ripr_init(Q, plurality)
#' flat_subnull(state)
#' @export
flat_subnull <- function(state) {
  rep(seq_along(state@atoms), block_sizes(state))
}

#' Number of atoms in each subnull
#' @keywords internal
#' @noRd
block_sizes <- function(state) vapply(state@atoms, ncol, integer(1))

#' Cut a flat vector into blocks of the given sizes
#'
#' Zero-length blocks are kept, so the result always has one element per
#' subnull.
#' @keywords internal
#' @noRd
split_by_sizes <- function(x, sizes) {
  ends <- cumsum(sizes)
  starts <- ends - sizes + 1L
  lapply(
    seq_along(sizes),
    \(i) if (sizes[i] == 0L) numeric(0) else unname(x[starts[i]:ends[i]])
  )
}

#' Redistribute a flat weight vector back into the per-subnull list
#' @keywords internal
#' @noRd
unflatten_weights <- function(state, w) split_by_sizes(w, block_sizes(state))

#' Where an atom appended to `subnull_index` lands in the flat ordering
#'
#' `add_atom` appends within a subnull's own block, so the new column sits at
#' the end of that block rather than at the end of the flat vector. The step
#' layer inserts the candidate at this index from the outset, which is what
#' keeps it from having to re-index anything afterwards.
#' @keywords internal
#' @noRd
insert_index <- function(state, subnull_index) {
  sum(block_sizes(state)[seq_len(subnull_index)]) + 1L
}

#' Replace the weights from a flat vector
#' @keywords internal
#' @noRd
set_weights <- function(state, w) {
  state@weights <- unflatten_weights(state, w)
  state
}

#' Append an atom to a subnull and set every weight at once
#'
#' `weights` is the full post-step flat vector of length `C + 1`, with the new
#' atom's entry at `insert_index`. Atoms and weights are jointly constrained
#' -- the lengths must match -- and S7 validates after every `@<-`, so the
#' assignment has to be a single call.
#' @keywords internal
#' @noRd
add_atom <- function(state, theta, subnull_index, weights) {
  atoms <- state@atoms
  atoms[[subnull_index]] <- cbind(
    atoms[[subnull_index]],
    theta,
    deparse.level = 0
  )
  S7::set_props(
    state,
    atoms = atoms,
    weights = split_by_sizes(weights, vapply(atoms, ncol, integer(1)))
  )
}

#' Redistribute a flat `(d, C)` atom matrix back into the per-subnull list
#' @keywords internal
#' @noRd
unflatten_atoms <- function(state, mat) {
  sizes <- block_sizes(state)
  ends <- cumsum(sizes)
  starts <- ends - sizes + 1L
  lapply(seq_along(sizes), function(i) {
    if (sizes[i] == 0L) {
      matrix(numeric(0), nrow = nrow(mat), ncol = 0L)
    } else {
      mat[, starts[i]:ends[i], drop = FALSE]
    }
  })
}

# --- Core quantities ----------------------------------------------------------

#' Log mixture density at the engine's nodes
#' @keywords internal
#' @noRd
log_p_at_nodes <- function(state, ld = NULL) {
  if (is.null(ld)) {
    ld <- compile_engine(state@engine)
  }
  row_logsumexp(add_by_col(ld(flat_atoms(state)), log(flat_weights(state))))
}

#' `KL(Q || P_W)` under the engine's quadrature rule
#' @keywords internal
#' @noRd
kl_divergence <- function(state, log_p = NULL, ld = NULL) {
  if (is.null(log_p)) {
    log_p <- log_p_at_nodes(state, ld)
  }
  expect_q(state@engine, state@engine@log_q - log_p)
}

# --- Trace ---------------------------------------------------------------------

#' @keywords internal
#' @noRd
empty_trace <- function() {
  data.frame(
    fw = integer(0),
    lb = integer(0),
    em = integer(0),
    weight = integer(0),
    phase = character(0),
    kl = numeric(0),
    gap = numeric(0),
    oracle_value = numeric(0),
    subnull = integer(0),
    step_size = numeric(0),
    direction = character(0),
    support_size = integer(0),
    max_weight = numeric(0),
    stringsAsFactors = FALSE
  )
}

bump <- function(state, which, by = 1L) {
  state@iters[[which]] <- state@iters[[which]] + as.integer(by)
  state
}

#' Append one row to the trace, and a snapshot if the control asks for it
#' @keywords internal
#' @noRd
record <- function(
  state,
  phase,
  kl,
  gap = NA_real_,
  oracle_value = NA_real_,
  subnull = NA_integer_,
  step_size = NA_real_,
  direction = NA_character_
) {
  w <- flat_weights(state)
  row <- data.frame(
    fw = state@iters[["fw"]],
    lb = state@iters[["lb"]],
    em = state@iters[["em"]],
    weight = state@iters[["weight"]],
    phase = phase,
    kl = kl,
    gap = gap,
    oracle_value = oracle_value,
    subnull = as.integer(subnull),
    step_size = step_size,
    direction = direction,
    support_size = length(w),
    max_weight = if (length(w)) max(w) else NA_real_,
    stringsAsFactors = FALSE
  )
  state@trace <- rbind(state@trace, row)

  state
}

#' Record the whole mixture alongside the trace
#'
#' A trace row is a handful of scalars; a snapshot copies every atom and weight,
#' so the two differ in cost by orders of magnitude and in frequency within a
#' single verb call. The verb decides when, via `wants_snapshot` -- `record`
#' cannot, since it does not know whether it sits partway through a `times` loop
#' or at the end of one.
#' @keywords internal
#' @noRd
snapshot_state <- function(state, phase) {
  state@snapshots <- c(
    state@snapshots,
    list(list(
      iters = state@iters,
      phase = phase,
      atoms = state@atoms,
      weights = state@weights
    ))
  )
  state
}


#' Should this iteration of a verb take a snapshot?
#'
#' `"step"` counts calls, so it fires only on the last iteration; `"all"` counts
#' iterations, so it fires on every one.
#' @param last Is this the final iteration of the current call?
#' @keywords internal
#' @noRd
wants_snapshot <- function(state, last) {
  switch(state@control$snapshot, none = FALSE, step = last, all = TRUE)
}
