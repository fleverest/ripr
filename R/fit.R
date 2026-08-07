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
#' @param alternative The alternative \eqn{Q}{Q}, an [outcome_distribution].
#' @param null A [null_model].
#' @param engine An engine spec, e.g. [exact_engine()].
#' @param atoms Optional list of `(d, n_i)` matrices, one per subnull. `NULL`
#'   places one atom per subnull by projecting the alternative's reference point,
#'   which is the sensible default and what the examples use. Empty subnulls are
#'   `ncol = 0` matrices, which the loop handles without a special case.
#' @param weights Optional list matching `atoms`; defaults to uniform.
#' @param control From [ripr_control()].
#' @return A [ripr_state] with no iterations run.
#' @export
ripr_init <- function(
  alternative,
  null,
  engine = exact_engine(),
  atoms = NULL,
  weights = NULL,
  control = ripr_control()
) {
  resolved <- resolve_engine(engine, alternative, null@family)

  if (is.null(atoms)) {
    # The alternative's modal parameter when there is one: the projection lands
    # near where W_1 puts its mass, which is where the RIPr will be. Falls back
    # to the family's canonical point, since Q need not be a mixture at all.
    ref <- if (S7_inherits(alternative, mixture)) {
      mode_parameter(alternative@mixing)
    } else {
      reference_parameter(null@family)
    }
    atoms <- lapply(null@subnulls, \(s) matrix(init_point(s, ref), ncol = 1L))
  }
  if (length(atoms) != length(null@subnulls)) {
    stop(
      "`atoms` must be a list with one element per subnull (",
      length(null@subnulls),
      "), not ",
      length(atoms),
      ".",
      call. = FALSE
    )
  }
  atoms <- lapply(atoms, as.matrix)

  sizes <- vapply(atoms, ncol, integer(1))
  if (sum(sizes) == 0L) {
    stop("at least one subnull must carry an atom.", call. = FALSE)
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
  record(state, phase = "init", kl = kl_divergence(state))
}


# --- Stopping predicates ------------------------------------------------------

#' Predicates for the `until` argument of a step verb
#'
#' Each returns a function of the state, so a verb stops early when its
#' condition holds: `em_step(1000, until = kl_flat(1e-12))` means "at most a
#' thousand sweeps".
#'
#' The three are not equally trustworthy.
#'
#' [gap_below()] and [support_gap_below()] are Frank--Wolfe gaps, over the whole
#' null and over the current support respectively, so each bounds how much KL is
#' still available at that scope.
#'
#' [gap_below()] is a predicate on the Frank--Wolfe gap as estimated (lower
#' bounded) by the oracle optimiser. It can never be fully trusted.
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
#' @name predicates
NULL

#' @describeIn predicates The change in KL over the last two trace rows.
#' @export
kl_flat <- function(tol = 1e-12) {
  function(state) {
    kl <- state@trace$kl
    length(kl) >= 2L && abs(diff(utils::tail(kl, 2L))) < tol
  }
}

#' @describeIn predicates The Frank--Wolfe gap over the whole null. Errors if no
#'   gap has been recorded at the state's current step counts, since a stale one
#'   would silently answer a question about an earlier iterate.
#' @export
gap_below <- function(tol = 1e-8) {
  function(state) {
    tr <- state@trace
    fresh <- tr$fw == state@iters[["fw"]] &
      tr$lb == state@iters[["lb"]] &
      tr$em == state@iters[["em"]] &
      tr$weight == state@iters[["weight"]]
    g <- tr$gap[fresh & !is.na(tr$gap)]
    if (!length(g)) {
      stop(
        "no gap recorded at the current step. Use `record_gap = TRUE`, or ",
        "`fw_step()`, which always records one.",
        call. = FALSE
      )
    }
    utils::tail(g, 1L) < tol
  }
}

#' @describeIn predicates The Frank--Wolfe gap over the current support,
#'   `max_c G(theta_c) - 1`. Always available, since it needs no oracle sweep.
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
#' The shared skeleton: `advance` takes the state and the compiled log-density
#' and returns the state plus whatever `record` needs.
#' @keywords internal
#' @noRd
run_steps <- function(state, times, until, counter, phase, advance) {
  rlang::check_number_whole(times, min = 1, max = 2147483647)
  ld <- compile_engine(state@engine)

  for (i in seq_len(times)) {
    stepped <- advance(state, ld)
    state <- bump(stepped$state, counter)
    state <- do.call(
      record,
      c(list(state, phase = phase), stepped$row)
    )
    if (wants_snapshot(state, i == times)) {
      state <- snapshot_state(state, phase)
    }
    if (!is.null(until) && isTRUE(until(state))) {
      break
    }
  }
  state
}


#' Frank--Wolfe step
#'
#' Maximises \eqn{G(\theta)}{G(theta)} over the null, then moves the iterate
#' towards the maximiser. See [oracles] for technical details.
#'
#' @param state A [ripr_state].
#' @param times Steps to take.
#' @param directions Any of `"forward"`, `"pairwise"`, `"away"`. More than one
#'   means each is tried and whichever reaches the lowest KL is taken.
#' @param size `"line-search"`, or `"fixed"` for the open-loop schedule.
#' @param until Optional predicate; see [predicates].
#' @return The updated [ripr_state].
#' @seealso [oracles], [predicates]
#' @export
fw_step <- function(
  state,
  times = 1L,
  directions = "forward",
  size = c("line-search", "fixed"),
  until = NULL
) {
  directions <- rlang::arg_match(
    directions,
    c("forward", "pairwise", "away"),
    multiple = TRUE
  )
  size <- rlang::arg_match(size)

  run_steps(state, times, until, "fw", "fw", function(state, ld) {
    log_p <- log_p_at_nodes(state, ld)
    found <- search_null(state, linear_oracle(state, log_p, ld))
    planned <- plan_step(
      state,
      log_p,
      ld,
      directions = directions,
      size = size,
      # Only Frank--Wolfe steps advance the schedule; EM sweeps between two of
      # them must not.
      gamma_fixed = schedule_gamma(length(flat_weights(state))),
      at = insert_index(state, found$subnull)
    )(found$theta)

    list(
      state = commit_step(state, found$theta, found$subnull, planned),
      row = list(
        kl = planned$kl,
        gap = found$value - 1,
        oracle_value = found$value,
        subnull = if (planned$uses_candidate) found$subnull else NA_integer_,
        step_size = planned$gamma,
        direction = planned$direction
      )
    )
  })
}


#' Li--Barron greedy step
#'
#' Scores each candidate by the KL it yields *after* the new weights are
#' chosen, so the step rule runs as an inner optimisation. Considerably more
#' expensive than [fw_step()]; see [oracles].
#'
#' @inheritParams fw_step
#' @param correct Re-solve every weight inside each candidate evaluation, using
#'   `lb_fc_tol` and `lb_fc_max_iter` from [ripr_control()].
#' @param record_gap Also sweep the Frank--Wolfe oracle, purely to record a gap.
#'   `FALSE` by default.
#' @return The updated [ripr_state].
#' @seealso [oracles], [predicates]
#' @export
lb_step <- function(
  state,
  times = 1L,
  directions = "forward",
  size = c("line-search", "fixed"),
  correct = FALSE,
  record_gap = FALSE,
  until = NULL
) {
  directions <- rlang::arg_match(
    directions,
    c("forward", "pairwise", "away"),
    multiple = TRUE
  )
  size <- rlang::arg_match(size)

  run_steps(state, times, until, "lb", "lb", function(state, ld) {
    log_p <- log_p_at_nodes(state, ld)
    obj <- nonlinear_oracle(
      state,
      log_p,
      ld,
      directions = directions,
      size = size,
      gamma_fixed = schedule_gamma(length(flat_weights(state))),
      correct = correct
    )
    found <- search_null(state, obj)
    planned <- plan_step(
      state,
      log_p,
      ld,
      directions = directions,
      size = size,
      gamma_fixed = schedule_gamma(length(flat_weights(state))),
      correct = correct,
      at = insert_index(state, found$subnull)
    )(found$theta)

    list(
      state = commit_step(state, found$theta, found$subnull, planned),
      row = list(
        kl = planned$kl,
        gap = if (record_gap) {
          linear_gap(state, log_p, ld, flat_atoms(state))
        } else {
          NA_real_
        },
        oracle_value = found$value,
        subnull = if (planned$uses_candidate) found$subnull else NA_integer_,
        step_size = planned$gamma,
        direction = planned$direction
      )
    )
  })
}


#' EM sweep
#'
#' Each sweep updates every weight, then moves every atom within its own
#' subnull. The support neither grows nor shrinks: only an oracle can add an
#' atom to the mixture.
#'
#' @inheritParams fw_step
#' @param record_gap Sweep the Frank--Wolfe oracle to record a gap. Off by
#'   default, since it costs a full oracle sweep per row.
#' @return The updated [ripr_state].
#' @seealso [oracles], [predicates]
#' @export
em_step <- function(state, times = 1L, record_gap = FALSE, until = NULL) {
  run_steps(state, times, until, "em", "em", function(state, ld) {
    stepped <- em_sweep(state, ld)
    log_p <- log_p_at_nodes(stepped, ld)
    list(
      state = stepped,
      row = list(
        kl = kl_divergence(stepped, log_p = log_p),
        gap = if (record_gap) {
          linear_gap(stepped, log_p, ld, flat_atoms(stepped))
        } else {
          NA_real_
        }
      )
    )
  })
}


#' Weight correction step
#'
#' \eqn{w_c \leftarrow w_c G(\theta_c)}{w_c <- w_c G(theta_c)} with the atoms
#' held fixed: the exact M-step for the weights, guaranteed monotone in KL.
#' Iterated to convergence this is the corrective half of fully-corrective
#' Frank--Wolfe, so `fw_step(1) |> weight_step(500)` is one FCFW iteration.
#' Note however that `lb_step(1) |> weight_step(500)` is not one Li--Barron
#' step with fully corrective weights incorporated as the inner optimisation.
#'
#' `until = support_gap_below(tol)` is the natural stopping rule, and is what
#' makes `times` a budget rather than a target. Expect to reach it: the rate
#' degrades as atoms crowd together, which is what Frank--Wolfe makes them do.
#'
#' @inheritParams em_step
#' @return The updated [ripr_state].
#' @seealso [predicates]
#' @export
weight_step <- function(state, times = 1L, record_gap = FALSE, until = NULL) {
  run_steps(state, times, until, "weight", "weight", function(state, ld) {
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
      row = list(
        kl = kl_divergence(stepped, log_p = new_log_p),
        # The residual is the gap over the current support, measured before the
        # sweep, so it says what this sweep had left to gain.
        oracle_value = sweep$residual + 1,
        gap = if (record_gap) {
          linear_gap(stepped, new_log_p, ld, flat_atoms(stepped))
        } else {
          NA_real_
        }
      )
    )
  })
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
#'   Turn it on to measure what refining and pruning did to the gap: `gap_fit`
#'   alone cannot show that, as it describes the state on the way in.
#' @param tol,max_iter Passed to the weight solve.
#' @param max_rounds Cap on refinement rounds.
#' @return A list with `W0` (a [finite_mixing]), `P_star` (a [mixture]), `kl` of
#'   the returned mixture, `gap_fit`, `gap_final`, `rounds`, `atoms`, `weights`,
#'   `subnull`, `trace` and `snapshots`.
#' @references
#'   \insertRef{FercoqGramfortSalmon2015}{ripr}
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

  mixing <- finite_mixing(
    components = flat_atoms(state)[, keep, drop = FALSE],
    weights = w[keep] / sum(w[keep])
  )
  log_p <- mixture_log_p(ld(mixing@components), mixing@weights)
  gaps <- state@trace$gap[!is.na(state@trace$gap)]

  list(
    W0 = mixing,
    P_star = mixture(mixing, engine@family),
    # Of what is being returned, not of the state it came from: refining and
    # pruning both change the mixture.
    kl = expect_q(engine, engine@log_q - log_p),
    gap_fit = if (length(gaps)) utils::tail(gaps, 1L) else NA_real_,
    gap_final = if (record_gap) {
      linear_gap(state, log_p, ld, mixing@components)
    } else {
      NA_real_
    },
    rounds = rounds,
    atoms = state@atoms,
    weights = state@weights,
    subnull = flat_subnull(state)[keep],
    trace = state@trace,
    snapshots = state@snapshots
  )
}
