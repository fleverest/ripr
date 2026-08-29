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
#' @param control From [ripr_control()].
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
  control = ripr_control()
) {
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
#'   fw_step(times = 25L, record_gap = TRUE, until = gap_below(1e-2))
#' gap_below(1e-2)(state)
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
        "no gap recorded at the current step. Pass `record_gap = TRUE` to the ",
        "verb that took it.",
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
#' @param record_gap Sweep the Frank--Wolfe oracle over the mixture the step
#'   *produced*, filling `gap` and `gap_theta`. Off by default: the sweep costs
#'   about as much as the step itself. The oracle value the step got for free is
#'   recorded regardless, in `oracle_value` and `oracle_theta`, that effectively
#'   measures the `gap` and `gap_theta` for the previous iterations.
#' @param until Optional predicate; see [predicates].
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
  directions = "forward",
  size = c("line-search", "fixed"),
  record_gap = FALSE,
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
      at = insert_index(state, found$part)
    )(found$theta)

    stepped <- commit_step(state, found$theta, found$part, planned)
    # `planned$log_p` is already the stepped mixture's density, and
    # `commit_step()` does not touch the engine, so the sweep needs no
    # recomputation beyond its own optimisation.
    swept <- if (record_gap) {
      linear_gap(stepped, planned$log_p, ld, flat_atoms(stepped))
    }

    list(
      state = stepped,
      row = list(
        kl = planned$kl,
        gap = if (is.null(swept)) NA_real_ else swept$gap,
        gap_theta = swept$theta,
        oracle_value = found$value,
        oracle_theta = found$theta,
        part = if (planned$uses_candidate) found$part else NA_integer_,
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
#' @param record_gap Sweep the Frank--Wolfe oracle over the mixture the step
#'   *produced*, filling `gap` and `gap_theta`. `FALSE` by default.
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
      at = insert_index(state, found$part)
    )(found$theta)

    stepped <- commit_step(state, found$theta, found$part, planned)
    swept <- if (record_gap) {
      linear_gap(stepped, planned$log_p, ld, flat_atoms(stepped))
    }

    list(
      state = stepped,
      row = list(
        kl = planned$kl,
        gap = if (is.null(swept)) NA_real_ else swept$gap,
        gap_theta = swept$theta,
        oracle_value = found$value,
        oracle_theta = found$theta,
        part = if (planned$uses_candidate) found$part else NA_integer_,
        step_size = planned$gamma,
        direction = planned$direction
      )
    )
  })
}


#' EM sweep
#'
#' Each sweep updates every weight, then moves every atom within its own
#' part. The support neither grows nor shrinks: only an oracle method such as
#' [fw_step()] or [lb_step()] can add an atom to the mixture.
#'
#' @inheritParams fw_step
#' @param record_gap Sweep the Frank--Wolfe oracle over the new mixture that the
#'   sweep *produced*, filling `gap` and `gap_theta`. Off by default, since it
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
  run_steps(state, times, until, "em", "em", function(state, ld) {
    stepped <- em_sweep(state, ld)
    log_p <- log_p_at_nodes(stepped, ld)
    swept <- if (record_gap) {
      linear_gap(stepped, log_p, ld, flat_atoms(stepped))
    }
    list(
      state = stepped,
      row = list(
        kl = kl_divergence(stepped, log_p = log_p),
        gap = if (is.null(swept)) NA_real_ else swept$gap,
        gap_theta = swept$theta
      )
    )
  })
}


#' Weight correction step
#'
#' \eqn{w_c \leftarrow w_c G(\theta_c)}{w_c <- w_c G(theta_c)} with the atoms
#' held fixed: the exact M-step for the weights, guaranteed monotone in KL.
#' Iterated to convergence this is the corrective half of fully-corrective
#' Frank--Wolfe, so `fw_step(1) |> weight_step(big_num)` is one FCFW iteration.
#' Note however that `lb_step(1) |> weight_step(big_num)` is not one Li--Barron
#' step with fully corrective weights incorporated as the inner optimisation.
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
    swept <- if (record_gap) {
      linear_gap(stepped, new_log_p, ld, flat_atoms(stepped))
    }
    list(
      state = stepped,
      row = list(
        kl = kl_divergence(stepped, log_p = new_log_p),
        # The residual is the gap over the current support, measured before the
        # sweep, so it says what this sweep had left to gain. It is a maximum
        # over the atoms rather than over the null, so no `oracle_theta` goes
        # with it -- the location is already in the mixture.
        oracle_value = sweep$residual + 1,
        gap = if (is.null(swept)) NA_real_ else swept$gap,
        gap_theta = swept$theta
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
#' @param tol,max_iter Passed to the weight solve.
#' @param max_rounds Cap on refinement rounds.
#' @return A list with `W0` (a [finite_dist]), `P_star` (a [mixture]), `kl` of
#'   the returned mixture, `gap_fit` (the last Frank--Wolfe gap recorded during
#'   fitting, `NA` if none was), `gap_final` (a fresh Frank--Wolfe gap over the
#'   returned mixture, `NA` unless `record_gap = TRUE`), `rounds`, `atoms`,
#'   `weights`, `part`, `trace` and `snapshots`.
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
  gaps <- state@trace$gap[!is.na(state@trace$gap)]

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
