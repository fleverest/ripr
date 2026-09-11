#' @include control.R state.R steps.R
NULL

# Fitting the RIPr mixture
#
# `ripr_init()` builds the starting point, the step verbs advance it, and
# `ripr_finish()` turns it into a mixture object. Everything exported by
# the package for the fitting process is here.
#
# There is no fixed pipeline. Which algorithm runs is decided by which verbs are
# called and in what order, so `ripr_control()` holds no algorithm settings:
#
#     fit <- ripr_init(Q, H0) |>
#       fw_step(10) |>
#       em_step(100) |>
#       ripr_finish()
#
# Every verb takes `times` and an optional `until` predicate, and every verb
# records one trace row per iteration.

#' Begin a RIPr fit
#'
#' @param alternative The alternative \eqn{Q}{Q}, an [distribution].
#' @param null A [null_model].
#' @param engine An engine spec, e.g. [exact_engine()].
#' @param atoms Optional list of `(d, n_i)` matrices, one per part of the
#'   null region. `NULL` places one atom per part by projecting the
#'   alternative's reference point, which is the sensible default and what the
#'   examples use. Empty parts are `ncol = 0` matrices, which the loop handles
#'   without a special case.
#' @param weights Optional list matching `atoms`; defaults to uniform.
#' @param record_gap Sweep the Frank--Wolfe oracle over the starting mixture,
#'   filling the `gap_after` columns on the init row. Off by default: the sweep
#'   costs about as much as a [fw_step()].
#' @param control From [ripr_control()]. Its `snapshot` setting applies here as
#'   it does for other verbs.
#' @return A [ripr_state] with no iterations run.
#' @examples
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
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
  started <- proc.time()[["elapsed"]]
  resolved <- resolve_engine(engine, alternative, null@family)

  if (is.null(atoms)) {
    # The alternative's modal parameter when there is one: the projection lands
    # near where W_1 puts its mass, which is where the RIPr will be. Falls back
    # to the family's canonical point, since Q need not be a mixture at all.
    ref <- if (S7_inherits(alternative, mixture)) {
      reference_point(alternative@mixing)
    } else {
      reference_point(null@family)
    }
    atoms <- lapply(
      parts(null@region),
      \(s) matrix(init_point(s, ref), ncol = 1L)
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
  atoms <- lapply(atoms, as.matrix)

  sizes <- vapply(atoms, ncol, integer(1))
  if (sum(sizes) == 0L) {
    stop("at least one part must carry an atom.", call. = FALSE)
  }
  if (is.null(weights)) {
    weights <- lapply(sizes, \(n) rep(1 / sum(sizes), n))
  }

  state <- ripr_state(
    atoms = atoms,
    weights = weights,
    alternative = alternative,
    null = null,
    engine = resolved,
    control = control,
    trace = empty_trace(),
    snapshots = list(),
    iters = c(fw = 0L, lb = 0L, em = 0L, weight = 0L)
  )
  ld <- compile_engine(resolved)
  log_p <- log_p_at_nodes(state, ld)
  kl <- kl_divergence(state, log_p)
  # The clock stops here: the sweep below is a diagnostic the caller asked for,
  # and its cost belongs to `record_gap` rather than to initialising.
  elapsed <- proc.time()[["elapsed"]] - started
  state <- record(state, phase = "init", kl = kl, elapsed = elapsed)
  if (record_gap) {
    swept_at <- proc.time()[["elapsed"]]
    swept <- linear_gap(state, log_p, ld, flat_atoms(state))
    state <- fill_gap(
      state,
      swept$gap,
      swept$theta,
      swept$part,
      proc.time()[["elapsed"]] - swept_at
    )
  }
  if (wants_snapshot(state, last = TRUE)) {
    state <- snapshot_state(state, "init")
  }
  state
}


# --- Stopping predicates ------------------------------------------------------

#' Predicates for early stopping
#'
#' The `until` argument of each step verb accepts these predicates.
#'
#' Each returns a function of the state, tested before each step: a satisfied
#' predicate means no step is taken. For instance,
#' `em_step(1000, until = kl_flat(1e-12))` means "at most a thousand EM
#' sweeps, or until the KL decreases by less than 1e-12".
#'
#' [gap_below()] and [support_gap_below()] are estimates of the Frank--Wolfe
#' gap, over the whole null and over the current support respectively, so each
#' roughly tracks how much KL is still available at that scope.
#'
#' [kl_flat()] is a per-row difference and bounds nothing. Under linear
#' convergence at rate \eqn{\rho}{rho} it is \eqn{(1-\rho)}{(1 - rho)} times the
#' true suboptimality, and \eqn{\rho}{rho} runs close to 1 here, so a small
#' `dKL` means converged *or* crawling. Use it as a
#' budget.
#'
#' Convergence is a property of the current iterate, not a statement about the
#' total fit. Take another step and a converged state is no longer converged.
#' That is why these are predicates rather than a flag on the state.
#'
#' @param tol Threshold.
#' @return A function of a [ripr_state] returning `TRUE` or `FALSE`.
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' Q <- fam(c(0.4, 0.35, 0.25))
#' state <- ripr_init(Q, plurality) |>
#'   fw_step(times = 25L, until = gap_below(1e-2))
#' gap_below(1e-2)(state)
#' @name predicates
NULL

#' @describeIn predicates Checks whether the pending step would decrease KL by
#'   less than `tol`. It is recommended for [em_step()] in particular. When
#'   called externally on a state, it falls back to the change in KL over the
#'   last two steps.
#' @export
kl_flat <- function(tol = 1e-12) {
  function(state, stepped = NULL) {
    kl <- state@trace$kl
    if (!is.null(stepped)) {
      return(
        length(kl) >= 1L &&
          abs(utils::tail(stepped@trace$kl, 1L) - utils::tail(kl, 1L)) < tol
      )
    }
    length(kl) >= 2L && abs(diff(utils::tail(kl, 2L))) < tol
  }
}

#' @describeIn predicates Checks whether the Frank--Wolfe gap that the current
#'   mixture attains is below `tol`. Inside [fw_step()] the row is filled by
#'   the step's own oracle before the check; elsewhere it needs a recorded gap,
#'   so make sure to set `record_gap = TRUE`.
#' @export
gap_below <- function(tol = 1e-8) {
  function(state) {
    tr <- state@trace
    fresh <- tr$fw == state@iters[["fw"]] &
      tr$lb == state@iters[["lb"]] &
      tr$em == state@iters[["em"]] &
      tr$weight == state@iters[["weight"]]
    g <- tr$gap_after[fresh & !is.na(tr$gap_after)]
    if (!length(g)) {
      stop(
        "no gap recorded for the current mixture. Use this predicate as ",
        "`until` in a step verb, set `record_gap = TRUE` on `em_step()`, ",
        "`lb_step()` or `weight_step()`, or measure the final mixture with ",
        "`ripr_finish(record_gap = TRUE)`.",
        call. = FALSE
      )
    }
    utils::tail(g, 1L) < tol
  }
}

#' @describeIn predicates Checks whether the Frank--Wolfe gap on the current
#'   support, `max_c G(theta_c) - 1` exceeds `tol`.
#'   It tests that the *weights* are optimal for the atoms in the current
#'   support, rather than the current fit being optimal for the full null.
#' @export
support_gap_below <- function(tol = 1e-8) {
  function(state) {
    ld <- compile_engine(state@engine)
    ld_all <- ld(flat_atoms(state))
    log_p <- log_p_at_nodes(state, ld)
    g <- exp(col_logsumexp(ld_all - log_p + state@engine@log_w))
    max(g) - 1 < tol
  }
}


# --- Step verbs ---------------------------------------------------------------

#' Run a verb's loop, recording and snapshotting as it goes
#'
#' `advance(state, ld)` computes one step without committing it, returning the
#' stepped state, that state's log-density at the nodes, the row for `record`,
#' and optionally `before`: the input state after any trace write the step's
#' work produced. The loop adopts `before`, builds the state the step would
#' make, and asks `until(state, candidate)`: a break keeps the `before` write
#' and discards the candidate; otherwise the candidate is adopted.
#'
#' The clock spans `advance` alone: snapshots, `record_gap` sweeps and the
#' `until` predicate are excluded. An advance that consumed a read-back oracle
#' reports the stored search time as `reused`, which is added to its row's
#' `elapsed`.
#' @keywords internal
#' @noRd
run_steps <- function(
  state,
  times,
  until,
  phase,
  record_gap,
  advance
) {
  rlang::check_number_whole(times, min = 1, max = 2147483647)
  ld <- compile_engine(state@engine)

  for (i in seq_len(times)) {
    started <- proc.time()[["elapsed"]]
    stepped <- advance(state, ld)
    # A read-back adds the stored search time, so the row prices the search
    # its step consumed wherever that search physically ran.
    elapsed <- proc.time()[["elapsed"]] -
      started +
      (if (is.null(stepped$reused)) 0 else stepped$reused)
    if (!is.null(stepped$before)) {
      state <- stepped$before
    }
    candidate <- do.call(
      record,
      c(
        list(bump(stepped$state, phase), phase = phase, elapsed = elapsed),
        stepped$row
      )
    )
    if (!is.null(until) && isTRUE(call_until(until, state, candidate))) {
      break
    }
    state <- candidate
    if (record_gap) {
      swept_at <- proc.time()[["elapsed"]]
      swept <- linear_gap(state, stepped$log_p, ld, flat_atoms(state))
      state <- fill_gap(
        state,
        swept$gap,
        swept$theta,
        swept$part,
        proc.time()[["elapsed"]] - swept_at
      )
    }
    if (wants_snapshot(state, i == times)) {
      state <- snapshot_state(state, phase)
    }
  }
  state
}


#' Ask a predicate about the state and the step waiting to be taken
#'
#' A predicate with two or more formals also receives the state the step
#' would make; one with a single formal is asked about the current state
#' alone.
#' @keywords internal
#' @noRd
call_until <- function(until, state, candidate) {
  if (length(formals(until)) >= 2L) {
    until(state, candidate)
  } else {
    until(state)
  }
}


#' The linear oracle at the current mixture, read back or searched
#'
#' Reuses an oracle already recorded for the current mixture; otherwise
#' searches and writes the result to the last trace row via `fill_gap()`.
#' @keywords internal
#' @noRd
oracle_at <- function(state, log_p, ld) {
  found <- recorded_oracle(state)
  if (!is.null(found)) {
    return(list(
      state = state,
      found = found,
      reused = if (is.na(found$elapsed)) 0 else found$elapsed
    ))
  }
  started <- proc.time()[["elapsed"]]
  found <- search_null(state, linear_oracle(state, log_p, ld))
  searched <- proc.time()[["elapsed"]] - started
  list(
    state = fill_gap(state, found$value - 1, found$theta, found$part, searched),
    found = found,
    reused = 0
  )
}


#' Frank--Wolfe step
#'
#' Maximises \eqn{G(\theta)}{G(theta)} over the null, then moves the iterate
#' towards the maximiser. See [oracles] for technical details.
#'
#' `correct = TRUE` yields the fully-corrective Frank--Wolfe scheme: the step
#' is taken, then every weight is re-solved with the atoms held fixed. It is
#' effectively the same thing as `fw_step(1) |> weight_step(big_num)`, but as
#' a single verb rather than two.
#'
#' @param state A [ripr_state].
#' @param times Steps to take.
#' @param variant Which Frank--Wolfe algorithm to run. `"standard"` moves only
#'   towards the atom the oracle found. `"away-step"` may instead move away from
#'   the worst atom the mixture already carries, taking whichever of the two the
#'   linear model prefers (i.e. the one that maximises \eqn{\langle -\nabla f, d\rangle}{<-grad f, d>})
#'   and then chooses a step length in that direction. `"pairwise"` moves mass
#'   from that worst atom to the new oracle atom directly, leaving every other
#'   weight untouched.
#' @param size `"line-search"`, or `"fixed"` for the open-loop schedule.
#' @param correct Whether to use "fully-corrective" steps that re-solve for the
#'   weights with the atoms held fixed, run to `fc_tol` and `fc_max_iter` from
#'   [ripr_control()]. Can be quite expensive. Note that the step size found by
#'   the method selected via `size` is only used as a seed for the
#'   fully-corrective solve.
#' @param until (Optional) A predicate, tested before each step. It is given
#'   the current state, and when it accepts two arguments the state the step
#'   would make, so that a predicate can decide whether to stop early based
#'   on the change between the two states. See [predicates].
#' @return The updated [ripr_state].
#' @seealso [oracles], [predicates]
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' Q <- fam(c(0.4, 0.35, 0.25))
#' state <- ripr_init(Q, plurality) |> fw_step(times = 10L)
#' state@trace$kl
#'
#' # Where each step put its atom.
#' do.call(cbind, state@trace$oracle_theta[state@trace$phase == "fw"])
#' @export
fw_step <- function(
  state,
  times = 1L,
  variant = c("standard", "away-step", "pairwise"),
  size = c("line-search", "fixed"),
  correct = FALSE,
  until = NULL
) {
  directions <- variant_directions(rlang::arg_match(variant))
  size <- rlang::arg_match(size)

  run_steps(
    state,
    times,
    until,
    "fw",
    record_gap = FALSE,
    advance = function(state, ld) {
      log_p <- log_p_at_nodes(state, ld)
      oracle <- oracle_at(state, log_p, ld)
      state <- oracle$state
      found <- oracle$found
      planned <- plan_step(
        state,
        log_p,
        ld,
        directions = directions,
        size = size,
        gamma_fixed = schedule_gamma(schedule_index(state)),
        correct = correct,
        at = insert_index(state, found$part)
      )(found$theta)

      list(
        before = state,
        reused = oracle$reused,
        state = commit_step(state, found$theta, found$part, planned),
        log_p = planned$log_p,
        row = list(
          kl = planned$kl,
          oracle_value = found$value,
          oracle_theta = found$theta,
          part = if (planned$uses_candidate) found$part else NA_integer_,
          step_size = planned$gamma,
          direction = planned$direction
        )
      )
    }
  )
}


#' Li--Barron greedy step
#'
#' Scores each candidate by the KL it yields *after* the new weights are
#' chosen, so the step rule runs as an inner optimisation. Considerably more
#' expensive than [fw_step()]; see [oracles].
#'
#' @inheritParams fw_step
#' @param record_gap Sweep the Frank--Wolfe oracle over the mixture the step
#'   *produced*, filling the `gap_after` columns. `FALSE` by default.
#' @return The updated [ripr_state].
#' @seealso [oracles], [predicates]
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
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
      obj <- nonlinear_oracle(
        state,
        log_p,
        ld,
        size = size,
        gamma_fixed = schedule_gamma(schedule_index(state)),
        correct = correct
      )
      found <- search_null(state, obj)
      planned <- plan_step(
        state,
        log_p,
        ld,
        size = size,
        gamma_fixed = schedule_gamma(schedule_index(state)),
        correct = correct,
        at = insert_index(state, found$part)
      )(found$theta)

      list(
        state = commit_step(state, found$theta, found$part, planned),
        log_p = planned$log_p,
        row = list(
          kl = planned$kl,
          oracle_value = found$value,
          oracle_theta = found$theta,
          part = if (planned$uses_candidate) found$part else NA_integer_,
          step_size = planned$gamma,
          direction = planned$direction
        )
      )
    }
  )
}


#' EM sweep
#'
#' Each sweep updates every weight, then moves every atom within its own
#' part. The support neither grows nor shrinks: only an oracle method such as
#' [fw_step()] or [lb_step()] can add an atom to the mixture.
#'
#' @inheritParams fw_step
#' @param record_gap Sweep the Frank--Wolfe oracle over the mixture the sweep
#'   *produced*, filling the `gap_after` columns. Off by default, since it
#'   costs a full oracle sweep per row.
#' @return The updated [ripr_state].
#' @seealso [oracles], [predicates]
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
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
        row = list(kl = kl_divergence(stepped, log_p = log_p))
      )
    }
  )
}


#' Weight correction step
#'
#' \eqn{w_c \leftarrow w_c G(\theta_c)}{w_c <- w_c G(theta_c)} with the atoms
#' held fixed: the exact M-step for the weights, guaranteed monotone in KL.
#' Iterated to convergence this is the corrective half of fully-corrective
#' Frank--Wolfe, so `fw_step(1) |> weight_step(big_num)` is one FCFW iteration.
#' `fw_step(1, correct = TRUE)` is the same iteration taken in one verb, which
#' is what to reach for when a row should price the whole of it; this verb is
#' for correcting weights on its own account, with its own rows and its own
#' stopping rule. Note that `lb_step(1) |> weight_step(big_num)` is *not* one
#' Li--Barron step with fully corrective weights incorporated as the inner
#' optimisation.
#'
#' `until = support_gap_below(tol)` is the natural stopping rule, and is what
#' makes `times` a budget rather than a target. Expect to reach it: the rate
#' degrades as atoms crowd together, which is what Frank--Wolfe makes them do.
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
#'     simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
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
      ld_all <- ld(flat_atoms(state))
      log_p <- log_p_at_nodes(state, ld)
      sweep <- weight_sweep(
        ld_all,
        flat_weights(state),
        log_p,
        engine = state@engine
      )
      stepped <- set_weights(state, sweep$weights)
      new_log_p <- log_p_at_nodes(stepped, ld)
      list(
        state = stepped,
        log_p = new_log_p,
        row = list(
          kl = kl_divergence(stepped, log_p = new_log_p),
          # The residual is the gap over the current support, measured
          # before the sweep, so it says what this sweep had left to gain.
          # It is a maximum over the atoms rather than over the null, so no
          # `oracle_theta` goes with it -- the location is already in the
          # mixture.
          oracle_value = sweep$residual + 1
        )
      )
    }
  )
}


# --- Finishing ----------------------------------------------------------------

#' Turn a fitted state into a mixture
#'
#' Converts a fitted state into a mixture. By default that is all it does: the
#' weights come across as they stand and nothing is dropped.
#'
#' *Experimental feature*:
#'
#' Two optional refinements, independently switchable so their effect can be
#' measured. Both are off by default, so what is returned is what was fitted.
#' `reoptimise` re-solves the weights over the current atoms, atoms fixed.
#' `identify` tests each atom in ascending weight order and zeroes it if
#' removing it does not increase KL. With both on they alternate until a round
#' removes nothing, because each makes the other work better: removing an atom
#' frees mass the survivors should absorb, and the removal test renormalises
#' rather than re-optimises, so it understates how good a removal is until the
#' weights have caught up. Usually one or two rounds.
#'
#' Neither is safe earlier in a fit. Nothing can restore a zeroed atom, so run
#' mid-fit `identify` ratchets the support down; it is sound here only because
#' the atoms have stopped moving.
#'
#' `identify` never increases KL -- it compares two feasible points and keeps
#' the better -- but it does not *certify* that a removed atom is zero at the
#' optimum. A certified rule would come from the duality gap, in the manner of
#' the safe screening literature, and needs the dual of this problem deriving
#' first.
#'
#' `prune` then drops atoms at or below its value. At the default of `0` that
#' is nothing, unless `identify` ran, in which case it is exactly the atoms it
#' zeroed. A positive value drops more, which may raise KL.
#'
#' Refining lowers KL and can *raise* the gap: the weight solve optimises over
#' the current support, which need not be where \eqn{\sup G}{sup G} is small,
#' and dropping atoms leaves more of the null uncovered. Measured on a `K = 4`
#' problem, both refinements together took KL from 0.120 to 0.096 and the gap
#' from 0.43 to 1.28. Which matters depends on whether the mixture is wanted for
#' its fit or for a certificate resting on the gap.
#'
#' @param state A [ripr_state].
#' @param prune Drop atoms with weight at or below this. Must be below 1.
#' @param reoptimise Re-solve the weights before pruning. Off by default.
#' @param identify Zero atoms whose removal does not increase KL. Off by
#'   default.
#' @param record_gap Sweep the Frank--Wolfe oracle over the *returned* mixture,
#'   filling `gap_final`. Off by default, since it costs a full oracle sweep.
#' @param tol,max_iter Passed to the weight solve.
#' @param max_rounds Cap on refinement rounds.
#' @return A list with `W0` (a [finite_dist]), `P_star` (a [mixture]), `kl` of
#'   the returned mixture, `gap_fit` (the last Frank--Wolfe gap recorded during
#'   fitting, which may describe an earlier iterate than the returned mixture;
#'   `NA` if none was), `gap_final` (the Frank--Wolfe gap over the final
#'   mixture, `NA` unless `record_gap = TRUE`), `rounds`, `atoms`, `weights`
#'   `part`, `trace` and `snapshots`.
#' @references
#'   \insertRef{FercoqGramfortSalmon2015}{ripr}
#' @examples
#' set.seed(1)
#' fam <- multinomial_family(n_trials = 4L, k = 3L)
#' plurality <- null_model(
#'   fam,
#'   list(
#'     simplex_region(vertices = cbind(c(0.5, 0.5, 0), c(0, 1, 0), c(0, 0, 1))),
#'     simplex_region(vertices = cbind(c(0.5, 0, 0.5), c(0, 1, 0), c(0, 0, 1)))
#'   )
#' )
#' Q <- fam(c(0.4, 0.35, 0.25))
#' state <- ripr_init(Q, plurality) |> fw_step(times = 10L)
#' fit <- ripr_finish(state, reoptimise = TRUE, identify = TRUE, record_gap = TRUE)
#' fit$kl
#' fit$gap_fit
#' fit$gap_final
#' @export
ripr_finish <- function(
  state,
  prune = 0,
  reoptimise = FALSE,
  identify = FALSE,
  record_gap = FALSE,
  tol = 1e-10,
  max_iter = 500L,
  max_rounds = 10L
) {
  rlang::check_number_decimal(prune)
  rlang::check_bool(reoptimise)
  rlang::check_bool(identify)
  rlang::check_bool(record_gap)
  if (prune >= 1) {
    stop("`prune` must be below 1; every weight is at most 1.", call. = FALSE)
  }

  engine <- state@engine
  ld <- compile_engine(engine)
  ld_all <- ld(flat_atoms(state))
  w <- flat_weights(state)

  # One round only when a single refinement is on: with nothing removed there is
  # nothing for a second solve to absorb, and with no solve the removal test
  # sees the same weights every time.
  rounds <- 0L
  if (reoptimise || identify) {
    for (round in seq_len(if (reoptimise && identify) max_rounds else 1L)) {
      rounds <- round
      if (reoptimise) {
        w <- solve_weights(ld_all, w, engine, tol = tol, max_iter = max_iter)
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

  keep <- w > prune
  if (!any(keep)) {
    stop(
      "no atom has weight above `prune` (",
      prune,
      "); the largest is ",
      signif(max(w), 3),
      ".",
      call. = FALSE
    )
  }

  mixing <- finite_dist(
    components = flat_atoms(state)[, keep, drop = FALSE],
    weights = w[keep] / sum(w[keep])
  )
  log_p <- mixture_log_p(ld(mixing@components), mixing@weights)
  gaps <- state@trace$gap_after[!is.na(state@trace$gap_after)]

  list(
    W0 = mixing,
    P_star = mixture(engine@family, mixing),
    # Of what is being returned, not of the state it came from: refining and
    # pruning both change the mixture.
    kl = expect_q(engine, engine@log_q - log_p),
    gap_fit = if (length(gaps)) utils::tail(gaps, 1L) else NA_real_,
    gap_final = if (record_gap) {
      linear_gap(state, log_p, ld, mixing@components)$gap
    } else {
      NA_real_
    },
    rounds = rounds,
    atoms = state@atoms,
    weights = state@weights,
    part = flat_part(state)[keep],
    trace = state@trace,
    snapshots = state@snapshots
  )
}
