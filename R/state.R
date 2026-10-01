#' @include null.R quadrature.R distribution.R mixing.R
NULL

# The optimiser's state and trace

#' State of a RIPr fit
#'
#' The mixing measure is a [finite_dist], whose atoms are in the order they
#' joined it (a new atom is appended last). `part` records which part of the
#' null each atom belongs to, since parts may overlap.
#'
#' [ripr_init()] builds the state and the step verbs advance it. The constructor
#' is not public.
#'
#' @param mixing The current mixture, a [finite_dist] over the null.
#' @param part Integer vector with one entry per atom of `mixing`: the part of
#'   the null region holding each atom.
#' @param alternative The alternative \eqn{Q}{Q}.
#' @param null A [null_model].
#' @param engine A resolved `quadrature`.
#' @param control From [ripr_control()].
#' @param trace_rows List of trace rows, one per step, each a named list;
#'   read them through `trace`.
#' @param snapshots List of recorded mixtures, each `list(step, phase, mixing,
#'   part)`: the trace row and verb it was taken at, and the mixture then.
#' @param oracle The linear oracle's verdict on the *current* mixture, or
#'   `NULL` when nothing has measured it: `list(value, theta, part, elapsed)`,
#'   with \eqn{\sup G = 1 + \mathrm{gap}}{sup G = 1 + gap}, its argmax, part,
#'   and search seconds. A cache: any step that changes the mixture clears it,
#'   a `record_gap = TRUE` sweep sets it, and the next [fw_step()] consumes it.
#'   [gap_below()] reads it.
#' @section Trace:
#'   `state@trace` is a read-only data frame with one row per step, the first
#'   being [ripr_init()]'s; see [oracles] for the columns.
#' @return A `ripr_state`.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = rbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' Q <- fam(c(0.4, 0.35, 0.25))
#' state <- ripr_init(Q, plurality)
#' state
#' state@mixing
#' state@part
#'
#' # Nothing has measured the starting mixture yet.
#' is.null(state@oracle)
#' # We may pass `record_gap = TRUE` to `ripr_init()` to sweep the oracle over
#' # the initial state.
#' state <- ripr_init(Q, plurality, record_gap = TRUE)
#' state@oracle$value - 1
#' state@trace$gap_after
#' @keywords internal
ripr_state <- new_class(
  "ripr_state",
  properties = list(
    mixing = finite_dist,
    part = class_integer,
    alternative = distribution,
    null = null_model,
    engine = quadrature,
    control = class_list,
    trace_rows = class_list,
    trace = new_property(
      S7::class_data.frame,
      getter = function(self) trace_frame(self@trace_rows)
    ),
    snapshots = class_list,
    oracle = class_any
  ),
  validator = function(self) {
    if (length(self@part) != n_atoms(self@mixing)) {
      return("`part` must have one element per atom of `mixing`")
    }
    if (!all(self@part %in% seq_len(n_parts(self@null@region)))) {
      return("`part` must index the parts of the null's region")
    }
    oracle_fields <- c("value", "theta", "part", "elapsed")
    if (
      !is.null(self@oracle) &&
        !(is.list(self@oracle) && all(oracle_fields %in% names(self@oracle)))
    ) {
      return("`oracle` must be `NULL` or a list(value, theta, part, elapsed)")
    }
    NULL
  }
)

#' @rdname ripr_state
#' @usage NULL
method(print, ripr_state) <- function(x, ...) {
  tr <- x@trace
  kl <- if (nrow(tr)) format(signif(utils::tail(tr$kl, 1L), 6L)) else NA

  # The cached oracle is the only gap known to describe the current mixture;
  # a recorded one without it belongs to an iterate that has since moved on.
  g_any <- tr$gap_after[!is.na(tr$gap_after)]
  gap <- if (!is.null(x@oracle)) {
    format(signif(x@oracle$value - 1, 3L))
  } else if (length(g_any)) {
    paste(
      format(signif(utils::tail(g_any, 1L), 3L)),
      "at an earlier iterate"
    )
  } else {
    "none recorded"
  }
  counts <- phase_counts(tr)

  cat("<ripr_state>\n")
  cat(
    "  null     ",
    class_name(x@null@family),
    " in ",
    count_label(n_parts(x@null@region), "part"),
    "\n",
    sep = ""
  )
  cat(
    "  iterate  ",
    count_label(n_atoms(x@mixing), "atom"),
    if (!is.na(kl)) paste0(", KL ", kl),
    "\n",
    sep = ""
  )
  cat("  gap      ", gap, "\n", sep = "")
  cat(
    "  steps    ",
    paste(names(counts), counts, collapse = ", "),
    sprintf(" (%.2gs)", sum(tr$elapsed, na.rm = TRUE)),
    "\n",
    sep = ""
  )
  cat(
    "  trace    ",
    count_label(nrow(tr), "row"),
    ", ",
    count_label(length(x@snapshots), "snapshot"),
    "\n",
    sep = ""
  )
  invisible(x)
}


#' @rdname ripr_state
#' @usage NULL
method(format, ripr_state) <- function(x, ...) {
  kl <- last_row(x)$kl
  sprintf(
    "ripr_state: %s, %s",
    count_label(n_atoms(x@mixing), "atom"),
    if (is.null(kl)) "no trace" else paste0("KL ", format(signif(kl, 6L)))
  )
}


# --- The mixture --------------------------------------------------------------

#' Replace the mixture and parts in one call
#'
#' Atoms, weights and parts must agree in length, so they change together.
#' @keywords internal
#' @noRd
set_mixture <- function(
  state,
  atoms = state@mixing@atoms,
  weights = state@mixing@weights,
  part = state@part
) {
  S7::set_props(
    state,
    mixing = finite_dist(atoms = atoms, weights = weights),
    part = part
  )
}

#' Append an atom; `weights` is the full length-`C + 1` vector, new atom last
#' @keywords internal
#' @noRd
add_atom <- function(state, theta, part_index, weights) {
  set_mixture(
    state,
    atoms = rbind(state@mixing@atoms, theta, deparse.level = 0),
    weights = weights,
    part = c(state@part, as.integer(part_index))
  )
}

#' Remove every atom the step left with no weight
#'
#' Not [prune()]: that would renormalise, and knows nothing of `part`.
#' @keywords internal
#' @noRd
drop_empty <- function(state) {
  w <- state@mixing@weights
  keep <- w > 0
  if (all(keep)) {
    return(state)
  }
  set_mixture(
    state,
    atoms = state@mixing@atoms[keep, , drop = FALSE],
    weights = w[keep],
    part = state@part[keep]
  )
}

# --- Core quantities ----------------------------------------------------------

#' Log density of a mixture with weights `w` over the atoms of `ld_all`
#' @keywords internal
#' @noRd
mixture_log_p <- function(ld_all, w) {
  row_logsumexp(add_by_col(ld_all, log(pmax(w, 0))))
}

#' Log mixture density at the engine's nodes
#' @keywords internal
#' @noRd
log_p_at_nodes <- function(state, ld) {
  mixture_log_p(ld(state@mixing@atoms), state@mixing@weights)
}

#' `KL(Q || P)` via the engine's quadrature, `P` given by log density at nodes
#' @keywords internal
#' @noRd
kl_at <- function(engine, log_p) expect_q(engine, engine@log_q - log_p)


#' Returns a function giving wall-clock seconds since the call
#' @keywords internal
#' @noRd
stopwatch <- function() {
  started <- proc.time()[["elapsed"]]
  function() proc.time()[["elapsed"]] - started
}

# --- Trace ---------------------------------------------------------------------

#' One trace row: the single place the columns are listed, in print order
#'
#' The columns are documented under "What a trace row records" in `?oracles`.
#' `w` is the produced mixture's weights.
#' @keywords internal
#' @noRd
trace_row <- function(
  step,
  w,
  phase,
  kl,
  oracle_value = NA_real_,
  oracle_theta = NULL,
  part = NA_integer_,
  step_size = NA_real_,
  direction = NA_character_,
  elapsed = NA_real_
) {
  list(
    step = as.integer(step),
    phase = phase,
    kl = kl,
    gap_after = NA_real_,
    gap_after_theta = theta_cell(NULL),
    gap_after_part = NA_integer_,
    gap_after_elapsed = NA_real_,
    oracle_value = oracle_value,
    oracle_theta = theta_cell(oracle_theta),
    part = as.integer(part),
    step_size = step_size,
    direction = direction,
    support_size = length(w),
    max_weight = if (length(w)) max(w) else NA_real_,
    elapsed = elapsed
  )
}

#' A trace row of missing values: each column's name and type
#' @keywords internal
#' @noRd
trace_prototype <- function() {
  trace_row(0L, numeric(0), phase = NA_character_, kl = NA_real_)
}

#' Bind trace rows into a data frame, one pass per column
#' @keywords internal
#' @noRd
trace_frame <- function(rows) {
  proto <- trace_prototype()
  cols <- lapply(names(proto), function(nm) {
    if (nm %in% theta_columns) {
      lapply(rows, .subset2, nm)
    } else {
      vapply(rows, .subset2, proto[[nm]], nm)
    }
  })
  names(cols) <- names(proto)
  structure(
    cols,
    class = "data.frame",
    row.names = .set_row_names(length(rows))
  )
}

#' The list columns of a trace
#' @keywords internal
#' @noRd
theta_columns <- c("gap_after_theta", "oracle_theta")

#' The last trace row, or `NULL`, without building the data frame
#' @keywords internal
#' @noRd
last_row <- function(state) {
  rows <- state@trace_rows
  if (length(rows)) rows[[length(rows)]] else NULL
}

#' Steps taken so far by verb (`fw`, `lb`, `em`, `weight`), from `phase`;
#' the `init` row is not counted
#' @keywords internal
#' @noRd
phase_counts <- function(trace) {
  vapply(
    c(fw = "fw", lb = "lb", em = "em", weight = "weight"),
    \(p) sum(trace$phase == p),
    integer(1)
  )
}

#' Append one row to the trace, describing the current mixture
#'
#' Arguments are those of `trace_row()` after `step` and `w`. A new row means a
#' new mixture, so this clears the cached oracle.
#' @keywords internal
#' @noRd
record <- function(state, ...) {
  rows <- state@trace_rows
  n <- length(rows)
  rows[[n + 1L]] <- trace_row(n, state@mixing@weights, ...)
  # Neither property can invalidate the state; skip the validator.
  S7::set_props(state, trace_rows = rows, oracle = NULL, .check = FALSE)
}

#' One cell of a `theta` list column, `NA` when there is nothing to say
#' @keywords internal
#' @noRd
theta_cell <- function(theta) {
  if (is.null(theta)) NA else theta
}

#' Write an oracle result into the last row's `gap_after` columns
#'
#' The last row produced the current mixture, so this is its Frank--Wolfe gap.
#' Never overwrites a recorded gap.
#' @keywords internal
#' @noRd
fill_gap <- function(state, gap, theta, part, elapsed) {
  rows <- state@trace_rows
  i <- length(rows)
  if (!i || !is.na(rows[[i]]$gap_after)) {
    return(state)
  }
  # Element by element: `[<-` with a list would spread a list-valued `theta`.
  rows[[i]][["gap_after"]] <- gap
  rows[[i]]["gap_after_theta"] <- list(theta_cell(theta))
  rows[[i]][["gap_after_part"]] <- as.integer(part)
  rows[[i]][["gap_after_elapsed"]] <- elapsed
  S7::set_props(state, trace_rows = rows, .check = FALSE)
}


#' Record the whole mixture alongside the trace
#'
#' The verb decides when, via `wants_snapshot()`.
#' @keywords internal
#' @noRd
snapshot_state <- function(state, phase) {
  state@snapshots <- c(
    state@snapshots,
    list(list(
      step = last_row(state)$step,
      phase = phase,
      mixing = state@mixing,
      part = state@part
    ))
  )
  state
}


#' Should this iteration of a verb take a snapshot? `last`: final iteration
#' of the call
#' @keywords internal
#' @noRd
wants_snapshot <- function(state, last) {
  switch(state@control$snapshot, none = FALSE, step = last, all = TRUE)
}
