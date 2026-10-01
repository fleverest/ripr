#' @include control.R mixing.R state.R steps.R
NULL

# Fitting the RIPr mixture. There is no fixed pipeline: the algorithm is the
# sequence of verbs called.
#
#     fit <- ripr_init(Q, H0) |>
#       fw_step(10) |>
#       em_step(100) |>
#       ripr_finish()
#
# Every verb takes `times` and `until`, and records one trace row per iteration.

#' Begin a RIPr fit
#'
#' @param alternative The alternative \eqn{Q}{Q}, an [distribution].
#' @param null A [null_model].
#' @param engine An engine spec, e.g. [exact_engine()].
#' @param atoms Optional list of `(n_i, d)` matrices, one atom per row and one
#'   matrix per part of the null region (`nrow = 0` for an empty part). `NULL`
#'   places one atom per part by projecting the alternative's reference point.
#' @param weights Optional list matching `atoms`, one weight per atom; the
#'   weights over all parts together must sum to 1. Defaults to uniform.
#' @param record_gap Sweep the Frank--Wolfe oracle over the starting mixture,
#'   filling the init row's `gap_after` columns and caching the result for the
#'   first [fw_step()]. Off by default: it costs about one [fw_step()].
#' @param control From [ripr_control()].
#' @return A [ripr_state] with no iterations run.
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
#' ripr_init(Q, plurality)
#' @export
ripr_init <- function(
  alternative,
  null,
  engine = exact_engine(),
  atoms = NULL,
  weights = NULL,
  record_gap = FALSE,
  control = ripr_control()
) {
  clock <- stopwatch()
  resolved <- resolve_engine(engine, alternative, null@family)

  if (is.null(atoms)) {
    # Seed near where W_1 puts its mass. Q need not be a mixture, hence the
    # family fallback.
    ref <- if (S7_inherits(alternative, mixture)) {
      reference_point(alternative@mixing)
    } else {
      reference_point(null@family)
    }
    atoms <- lapply(
      parts(null@region),
      \(s) matrix(init_point(s, ref), nrow = 1L)
    )
  }
  if (length(atoms) != n_parts(null@region)) {
    stop(
      "`atoms` must be a list with one element per part (",
      n_parts(null@region),
      "), not ",
      length(atoms),
      ".",
      call. = FALSE
    )
  }
  atoms <- lapply(atoms, as_row_matrix)
  prts <- parts(null@region)
  for (i in seq_along(atoms)) {
    inside <- apply(atoms[[i]], 1L, \(theta) contains(prts[[i]], theta))
    if (!all(as.logical(inside))) {
      stop(
        "every atom must lie in its own part; an atom given for part ",
        i,
        " does not.",
        call. = FALSE
      )
    }
  }

  sizes <- vapply(atoms, nrow, integer(1))
  if (sum(sizes) == 0L) {
    stop("at least one part must carry an atom.", call. = FALSE)
  }
  weights <- if (is.null(weights)) {
    rep(1 / sum(sizes), sum(sizes))
  } else {
    unlist(weights, use.names = FALSE)
  }

  state <- ripr_state(
    mixing = finite_dist(
      # An empty part may be any `nrow = 0` matrix, whatever its column count.
      atoms = do.call(rbind, unname(atoms[sizes > 0L])),
      weights = weights
    ),
    part = rep(seq_along(atoms), sizes),
    alternative = alternative,
    null = null,
    engine = resolved,
    control = control
  )
  ld <- compile_engine(resolved)
  log_p <- log_p_at_nodes(state, ld)
  # Exclude the optional `record_gap` sweep from the init row's time.
  state <- record(
    state,
    phase = "init",
    kl = kl_at(state@engine, log_p),
    elapsed = clock()
  )
  if (record_gap) {
    state <- search_gap(state, log_p, ld)
  }
  if (wants_snapshot(state, last = TRUE)) {
    state <- snapshot_state(state, "init")
  }
  state
}


# --- Stopping predicates ------------------------------------------------------

#' Predicates for early stopping
#'
#' The `until` argument of each step verb accepts these predicates: functions
#' `function(state, candidate)` tested before each step, where `candidate` is
#' the state the pending step would produce (its trace row already appended).
#' If satisfied, the step is not taken and the verb returns `state`. So
#' `em_step(1000, until = kl_flat(1e-12))` means "at most a thousand EM sweeps,
#' or until a sweep would decrease KL by less than 1e-12". Your own predicate
#' must accept both arguments, e.g.
#' `function(state, candidate) nrow(candidate@trace) > 20`.
#'
#' [gap_below()] and [support_gap_below()] estimate the Frank--Wolfe gap of
#' `state` over the whole null and over the current support respectively. They
#' ignore `candidate`, so can be asked of one state as `gap_below(tol)(state)`.
#'
#' [kl_flat()] is a per-step difference and bounds nothing: under slow linear
#' convergence a small `dKL` means converged *or* crawling. Use it as a budget.
#'
#' @param tol Threshold.
#' @return A function of `(state, candidate)`, two [ripr_state]s, returning
#'   `TRUE` or `FALSE`.
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = rbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' Q <- fam(c(0.4, 0.35, 0.25))
#' state <- ripr_init(Q, plurality) |>
#'   fw_step(times = 25L, until = gap_below(1e-2))
#' gap_below(1e-2)(state)
#'
#' # Would one more EM sweep leave the KL flat?
#' kl_flat(1e-6)(state, em_step(state))
#'
#' # A predicate of your own sees both states: here, stop before the support
#' # would grow past a dozen atoms.
#' too_many_atoms <- function(state, candidate) {
#'   n_atoms(candidate@mixing) > 12L
#' }
#' fw_step(state, 10L, until = too_many_atoms)
#' @name predicates
NULL

#' @describeIn predicates Checks whether the pending step would change the KL
#'   by less than `tol`, comparing the last trace row of `candidate` with that
#'   of `state`. It is recommended for [em_step()] in particular.
#' @export
kl_flat <- function(tol = 1e-12) {
  function(state, candidate) {
    abs(last_row(candidate)$kl - last_row(state)$kl) < tol
  }
}

#' @describeIn predicates Checks whether the Frank--Wolfe gap of `state` is
#'   below `tol`, reading the oracle cached in `state@oracle`. Inside
#'   [fw_step()] the cache is always there; elsewhere set `record_gap = TRUE`
#'   on the previous verb.
#' @export
gap_below <- function(tol = 1e-8) {
  function(state, candidate) {
    if (is.null(state@oracle)) {
      stop(
        "no gap recorded for the current mixture. Set `record_gap = TRUE` ",
        "on the previous verb, or use this predicate as `until` in `fw_step()`.",
        call. = FALSE
      )
    }
    state@oracle$value - 1 < tol
  }
}

#' @describeIn predicates Checks whether the Frank--Wolfe gap on the current
#'   support, `max_c G(theta_c) - 1`, is below `tol`: whether the *weights* are
#'   optimal for the current atoms, not whether the fit is optimal.
#' @export
support_gap_below <- function(tol = 1e-8) {
  function(state, candidate) {
    ld_all <- compile_engine(state@engine)(state@mixing@atoms)
    log_p <- mixture_log_p(ld_all, state@mixing@weights)
    max(atom_g(ld_all, log_p, state@engine)) - 1 < tol
  }
}


# --- Step verbs ---------------------------------------------------------------

#' Run a verb's loop, recording and snapshotting as it goes
#'
#' `advance(state, ld)` returns `list(state, log_p, row)` for one uncommitted
#' step; the candidate is adopted unless `until(state, candidate)`.
#'
#' For `uses_oracle` verbs the linear oracle is searched before the step if not
#' cached, so `advance` only reads `state@oracle`; this also fills the last
#' row's `gap_after` and lets `gap_below()` answer inside the loop. The clock
#' spans `advance` alone, plus the consumed search's `elapsed`.
#' @keywords internal
#' @noRd
run_steps <- function(
  state,
  times,
  until,
  phase,
  record_gap,
  advance,
  uses_oracle = FALSE
) {
  rlang::check_number_whole(times, min = 1, max = 2147483647)
  check_until(until)
  ld <- compile_engine(state@engine)

  for (i in seq_len(times)) {
    if (uses_oracle && is.null(state@oracle)) {
      state <- search_gap(state, log_p_at_nodes(state, ld), ld)
    }
    clock <- stopwatch()
    stepped <- advance(state, ld)
    elapsed <- clock() + if (uses_oracle) state@oracle$elapsed else 0
    candidate <- do.call(
      record,
      c(list(stepped$state, phase = phase, elapsed = elapsed), stepped$row)
    )
    if (!is.null(until) && isTRUE(until(state, candidate))) {
      break
    }
    state <- candidate
    if (record_gap) {
      state <- search_gap(state, stepped$log_p, ld)
    }
    if (wants_snapshot(state, i == times)) {
      state <- snapshot_state(state, phase)
    }
  }
  state
}


#' Refuse an `until` that cannot take two states
#' @keywords internal
#' @noRd
check_until <- function(until) {
  if (is.null(until)) {
    return(invisible())
  }
  args <- if (is.function(until)) formals(args(until)) else NULL
  if (!is.function(until) || (!"..." %in% names(args) && length(args) < 2L)) {
    stop(
      "`until` must be a function of `(state, candidate)`; see ?predicates.",
      call. = FALSE
    )
  }
  invisible()
}


#' Search the linear oracle over the current mixture and cache the result
#'
#' Sets `state@oracle` and fills the last row's `gap_after`. The only place the
#' linear oracle runs during a fit. `log_p` is the current mixture's.
#' @keywords internal
#' @noRd
search_gap <- function(state, log_p, ld) {
  clock <- stopwatch()
  found <- search_null(state, linear_oracle(state, log_p, ld))
  elapsed <- clock()
  state <- fill_gap(state, found$value - 1, found$theta, found$part, elapsed)
  state@oracle <- list(
    value = found$value,
    theta = found$theta,
    part = found$part,
    elapsed = elapsed
  )
  state
}


#' Step towards an oracle's `found` `theta`/`part`, as `run_steps()` expects
#' of `advance`; `...` goes to `plan_step()`
#' @keywords internal
#' @noRd
step_towards <- function(state, found, oracle_value, log_p, ld, ...) {
  planned <- plan_step(state, log_p, ld, ...)(found$theta)
  list(
    state = commit_step(state, found$theta, found$part, planned),
    log_p = planned$log_p,
    row = list(
      kl = planned$kl,
      oracle_value = oracle_value,
      oracle_theta = found$theta,
      part = if (planned$uses_candidate) found$part else NA_integer_,
      step_size = planned$gamma,
      direction = planned$direction
    )
  )
}


#' Frank--Wolfe step
#'
#' Maximises \eqn{G(\theta)}{G(theta)} over the null, then moves the iterate
#' towards the maximiser; see [oracles]. The search is skipped if
#' `state@oracle` already holds it (e.g. after `record_gap = TRUE`), so each
#' step costs one search.
#'
#' `correct = TRUE` gives fully-corrective Frank--Wolfe: after each step every
#' weight is re-solved to optimality over the current atoms.
#'
#' @param state A [ripr_state].
#' @param times Steps to take.
#' @param variant Which Frank--Wolfe algorithm to run. `"standard"` moves only
#'   towards the oracle's atom. `"away-step"` may instead move away from the
#'   worst current atom, whichever maximises
#'   \eqn{\langle -\nabla f, d\rangle}{<-grad f, d>}. `"pairwise"` moves mass
#'   from that worst atom directly to the new one.
#' @param size `"line-search"`, or `"fixed"` for the open-loop schedule.
#' @param correct Re-solve the weights with atoms fixed after each step, to
#'   `fc_tol` and `fc_max_iter` from [ripr_control()]. Can be expensive; the
#'   `size` step only seeds the solve.
#' @param record_gap Sweep the Frank--Wolfe oracle over the produced mixture
#'   after each step, filling `gap_after`. Earlier rows fill anyway from the
#'   next step's search, so `TRUE` buys the final row's gap for about one extra
#'   step.
#' @param until Optional predicate `function(state, candidate)` tested before
#'   each step; see [predicates].
#' @return The updated [ripr_state].
#' @seealso [oracles], [predicates]
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = rbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' Q <- fam(c(0.4, 0.35, 0.25))
#' state <- ripr_init(Q, plurality) |> fw_step(times = 10L)
#' state@trace$kl
#'
#' # Where each step put its atom.
#' do.call(rbind, state@trace$oracle_theta[state@trace$phase == "fw"])
#' @export
fw_step <- function(
  state,
  times = 1L,
  variant = c("standard", "away-step", "pairwise"),
  size = c("line-search", "fixed"),
  correct = FALSE,
  record_gap = FALSE,
  until = NULL
) {
  directions <- variant_directions(rlang::arg_match(variant))
  size <- rlang::arg_match(size)

  run_steps(
    state,
    times,
    until,
    "fw",
    record_gap = record_gap,
    uses_oracle = TRUE,
    advance = function(state, ld) {
      step_towards(
        state,
        state@oracle,
        # Already recorded as the previous row's `1 + gap_after`.
        oracle_value = NA_real_,
        log_p_at_nodes(state, ld),
        ld,
        directions = directions,
        size = size,
        gamma_fixed = schedule_gamma(state),
        correct = correct
      )
    }
  )
}


#' Li--Barron greedy step
#'
#' Scores each candidate by the KL left *after* its weight is chosen, so the
#' step rule is an inner optimisation. Considerably more expensive than
#' [fw_step()]; see [oracles].
#'
#' @inheritParams fw_step
#' @param record_gap Sweep the Frank--Wolfe oracle over the produced mixture,
#'   filling `gap_after` and caching it for [gap_below()] and [fw_step()].
#'   `FALSE` by default.
#' @return The updated [ripr_state].
#' @seealso [oracles], [predicates]
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = rbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' Q <- fam(c(0.4, 0.35, 0.25))
#' state <- ripr_init(Q, plurality) |> lb_step(times = 5L)
#' state@trace$kl
#' @export
lb_step <- function(
  state,
  times = 1L,
  size = c("line-search", "fixed"),
  correct = FALSE,
  record_gap = FALSE,
  until = NULL
) {
  size <- rlang::arg_match(size)

  run_steps(
    state,
    times,
    until,
    "lb",
    record_gap,
    function(state, ld) {
      log_p <- log_p_at_nodes(state, ld)
      gamma_fixed <- schedule_gamma(state)
      obj <- nonlinear_oracle(
        state,
        log_p,
        ld,
        size = size,
        gamma_fixed = gamma_fixed,
        correct = correct
      )
      found <- search_null(state, obj)
      step_towards(
        state,
        found,
        oracle_value = found$value,
        log_p,
        ld,
        size = size,
        gamma_fixed = gamma_fixed,
        correct = correct
      )
    }
  )
}


#' EM sweep
#'
#' Each sweep updates every weight, then moves every atom within its own part.
#' The support never changes size; only [fw_step()] or [lb_step()] add atoms.
#'
#' @inheritParams fw_step
#' @param record_gap Sweep the Frank--Wolfe oracle over the produced mixture,
#'   filling `gap_after` and caching it for [gap_below()] and [fw_step()].
#'   Off by default: a full oracle sweep per row.
#' @return The updated [ripr_state].
#' @seealso [oracles], [predicates]
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = rbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' Q <- fam(c(0.4, 0.35, 0.25))
#' state <- ripr_init(Q, plurality) |> fw_step(times = 3L) |> em_step(times = 10L)
#' state@trace$kl
#' @export
em_step <- function(state, times = 1L, record_gap = FALSE, until = NULL) {
  run_steps(
    state,
    times,
    until,
    "em",
    record_gap,
    function(state, ld) {
      stepped <- em_sweep(state, ld)
      log_p <- log_p_at_nodes(stepped, ld)
      list(
        state = stepped,
        log_p = log_p,
        row = list(kl = kl_at(stepped@engine, log_p))
      )
    }
  )
}


#' Weight correction step
#'
#' \eqn{w_c \leftarrow w_c G(\theta_c)}{w_c <- w_c G(theta_c)} with the atoms
#' held fixed: the exact M-step for the weights, monotone in KL. Repeated, it
#' converges only slowly to the optimal weights; to solve for optimal weights in
#' one go, the current API supports only `fw_step(correct = TRUE)` or
#' `ripr_finish(reoptimise = TRUE)`.
#'
#' Stop with `until = support_gap_below(tol)`, making `times` a budget.
#' Convergence slows as atoms crowd together.
#'
#' @inheritParams em_step
#' @return The updated [ripr_state].
#' @seealso [predicates]
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = rbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' Q <- fam(c(0.4, 0.35, 0.25))
#' state <- ripr_init(Q, plurality) |>
#'   fw_step(times = 5L) |>
#'   weight_step(times = 20L, until = support_gap_below(1e-6))
#' state@trace$kl
#' @export
weight_step <- function(state, times = 1L, record_gap = FALSE, until = NULL) {
  run_steps(
    state,
    times,
    until,
    "weight",
    record_gap,
    function(state, ld) {
      ld_all <- ld(state@mixing@atoms)
      w <- state@mixing@weights
      sweep <- weight_sweep(
        ld_all,
        w,
        mixture_log_p(ld_all, w),
        engine = state@engine
      )
      stepped <- set_mixture(state, weights = sweep$weights)
      new_log_p <- log_p_at_nodes(stepped, ld)
      list(
        state = stepped,
        log_p = new_log_p,
        row = list(
          kl = kl_at(stepped@engine, new_log_p),
          # The pre-sweep support gap; no `oracle_theta`, since its location
          # is already an atom.
          oracle_value = sweep$residual + 1
        )
      )
    }
  )
}


# --- Finishing ----------------------------------------------------------------

#' A finished RIPr fit
#'
#' [ripr_finish()] returns the fitted mixing measure, the mixture it induces,
#' fit diagnostics, and the state it came from. The constructor is internal.
#'
#' The trace, snapshots and pre-finish mixture stay on `fit@state`. Finishing
#' may reorder, reweight or drop atoms, so use `W0` and `part`, not
#' `fit@state@mixing`.
#'
#' @param W0 The fitted mixing measure, a [finite_dist] over the null.
#' @param P_star The mixture `W0` induces through the null's family, the
#'   approximate RIPr.
#' @param part Integer vector, the part of the null holding each atom of `W0`.
#'   Atoms are grouped in part order, so `part` is non-decreasing.
#' @param kl \eqn{KL(Q \| P^*)}{KL(Q || P*)} of the returned mixture, under
#'   the state's quadrature rule.
#' @param gap_fit The last Frank--Wolfe gap recorded during fitting, which may
#'   describe an earlier iterate than the returned mixture; `NA` if none was.
#' @param gap_final The Frank--Wolfe gap over the returned mixture, `NA` unless
#'   [ripr_finish()] was called with `record_gap = TRUE`.
#' @param rounds Number of refinement rounds [ripr_finish()] ran; `0` when
#'   neither `reoptimise` nor `identify` was on.
#' @param state The [ripr_state] the fit was finished from.
#' @return A `ripr_fit`.
#' @seealso [ripr_finish()]
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = rbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' Q <- fam(c(0.4, 0.35, 0.25))
#' fit <- ripr_init(Q, plurality) |>
#'   fw_step(times = 10L) |>
#'   ripr_finish(record_gap = TRUE)
#' fit
#' fit@W0
#' fit@part
#' utils::tail(fit@state@trace$kl)
ripr_fit <- new_class(
  "ripr_fit",
  properties = list(
    W0 = finite_dist,
    P_star = mixture,
    part = class_integer,
    kl = class_numeric,
    gap_fit = class_numeric,
    gap_final = class_numeric,
    rounds = class_integer,
    state = ripr_state
  ),
  validator = function(self) {
    if (length(self@part) != nrow(self@W0@atoms)) {
      return("`part` must have one entry per atom of `W0`")
    }
    scalars <- list(
      kl = self@kl,
      gap_fit = self@gap_fit,
      gap_final = self@gap_final,
      rounds = self@rounds
    )
    for (nm in names(scalars)) {
      if (length(scalars[[nm]]) != 1L) {
        return(sprintf("`%s` must be a single number", nm))
      }
    }
    NULL
  }
)


#' @description `print()` shows the returned mixture, its KL, the gaps that
#'   were measured, and how much history the state carries.
#' @rdname ripr_fit
#' @usage NULL
method(print, ripr_fit) <- function(x, ...) {
  tr <- x@state@trace
  gaps <- c(
    if (!is.na(x@gap_final)) {
      paste(format(signif(x@gap_final, 3L)), "final")
    },
    if (!is.na(x@gap_fit)) {
      paste(format(signif(x@gap_fit, 3L)), "while fitting")
    }
  )

  cat("<ripr_fit>\n")
  cat("  mixture  ", format(x@P_star), "\n", sep = "")
  cat("  KL       ", format(signif(x@kl, 6L)), "\n", sep = "")
  cat(
    "  gap      ",
    if (length(gaps)) paste(gaps, collapse = ", ") else "none recorded",
    "\n",
    sep = ""
  )
  cat(
    "  refined  ",
    if (x@rounds) count_label(x@rounds, "round") else "no",
    "\n",
    sep = ""
  )
  cat(
    "  trace    ",
    nrow(tr),
    ngettext(nrow(tr), " row", " rows"),
    ", ",
    length(x@state@snapshots),
    ngettext(length(x@state@snapshots), " snapshot", " snapshots"),
    "\n",
    sep = ""
  )
  invisible(x)
}


#' @description `format()` gives the atom count and KL on one line.
#' @rdname ripr_fit
#' @usage NULL
method(format, ripr_fit) <- function(x, ...) {
  sprintf(
    "ripr_fit: %s, KL %s",
    count_label(nrow(x@W0@atoms), "atom"),
    format(signif(x@kl, 6L))
  )
}


#' Turn a fitted state into a mixture
#'
#' Converts a fitted state into a mixture. By default the weights come across
#' as they stand and nothing is dropped.
#'
#' *Experimental*: two optional refinements, both off by default. `reoptimise`
#' re-solves the weights over the current atoms (to `fc_tol`/`fc_max_iter`).
#' `identify` zeroes, lightest first, each atom whose removal does not increase
#' KL. With both on they alternate until a round removes nothing. They are sound
#' only because the atoms have stopped moving; `identify` never increases KL but
#' does not *certify* that a removed atom is zero at the optimum.
#'
#' Atoms with weight at or below `threshold` are then dropped: by default only
#' those `identify` zeroed. A positive value may raise KL.
#'
#' Refining lowers KL but can *raise* the gap, since the weights are optimised
#' over the current support and dropping atoms leaves more of the null
#' uncovered. Choose according to whether you want the fit or a certificate
#' resting on the gap.
#'
#' @param state A [ripr_state].
#' @param threshold Drop atoms with weight at or below this. Must be below 1.
#' @param reoptimise Re-solve the weights before dropping any. Off by default.
#' @param identify Zero atoms whose removal does not increase KL. Off by
#'   default.
#' @param record_gap Sweep the Frank--Wolfe oracle over the *returned* mixture,
#'   filling `gap_final`. Off by default, since it costs a full oracle sweep.
#' @param max_rounds Cap on refinement rounds.
#' @return A [ripr_fit], whose `W0` is the mixing measure (a [finite_dist]),
#'   `P_star` the mixture it induces, and `kl`, `gap_fit` and `gap_final`
#'   describe it. The state it was finished from, with its trace and
#'   snapshots, is kept as `fit@state`.
#' @references
#'   \insertRef{FercoqGramfortSalmon2015}{ripr}
#' @seealso [ripr_fit]
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = rbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = rbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' Q <- fam(c(0.4, 0.35, 0.25))
#' state <- ripr_init(Q, plurality) |> fw_step(times = 10L)
#' fit <- ripr_finish(state, reoptimise = TRUE, identify = TRUE, record_gap = TRUE)
#' fit
#' fit@kl
#' fit@gap_fit
#' fit@gap_final
#' @export
ripr_finish <- function(
  state,
  threshold = 0,
  reoptimise = FALSE,
  identify = FALSE,
  record_gap = FALSE,
  max_rounds = 10L
) {
  rlang::check_number_decimal(threshold)
  rlang::check_bool(reoptimise)
  rlang::check_bool(identify)
  rlang::check_bool(record_gap)
  if (threshold >= 1) {
    stop(
      "`threshold` must be below 1; every weight is at most 1.",
      call. = FALSE
    )
  }

  engine <- state@engine
  ctl <- state@control
  ld <- compile_engine(engine)
  support <- atoms(state@mixing)
  ld_all <- ld(support)
  w <- weights(state@mixing)

  # Alternating only helps with both on: with one alone a second round would
  # repeat the first.
  rounds <- 0L
  if (reoptimise || identify) {
    for (round in seq_len(if (reoptimise && identify) max_rounds else 1L)) {
      rounds <- round
      if (reoptimise) {
        w <- solve_weights(
          ld_all,
          w,
          engine,
          tol = ctl$fc_tol,
          max_iter = ctl$fc_max_iter
        )
      }
      if (!identify) {
        break
      }
      refined <- identify_support(ld_all, w, engine)
      settled <- identical(which(refined > 0), which(w > 0))
      w <- refined
      if (settled) {
        break
      }
    }
  }

  # Group atoms by part, so `W0` lines up with `part`.
  by_order <- order(state@part)
  w <- w[by_order] / sum(w)
  mixing <- prune(
    finite_dist(atoms = support[by_order, , drop = FALSE], weights = w),
    threshold = threshold
  )
  log_p <- mixture_log_p(ld(mixing@atoms), mixing@weights)
  gaps <- state@trace$gap_after
  gaps <- gaps[!is.na(gaps)]

  ripr_fit(
    W0 = mixing,
    P_star = mixture(engine@family, mixing),
    part = state@part[by_order][w > threshold],
    # Of the returned mixture, not the state's.
    kl = kl_at(engine, log_p),
    gap_fit = if (length(gaps)) utils::tail(gaps, 1L) else NA_real_,
    gap_final = if (record_gap) {
      oracle <- linear_oracle(state, log_p, ld)
      search_null(state, oracle, seeds = mixing@atoms)$value - 1
    } else {
      NA_real_
    },
    rounds = as.integer(rounds),
    state = state
  )
}
