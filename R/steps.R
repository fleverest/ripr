#' @include null.R quadrature.R numerics.R state.R
NULL

# Oracles and step rules
#
# Oracles propose where a new atom should go, while step rules decide whether to
# add it, and how to update the weights correspondingly. EM steps are different
# in that they may move atom locations as well as changing the weights.
#
# Oracles and steps require `(ld_all, w)`; log-density for each quadrature node
# (e.g. the full support for exact) and the weights over them. The atoms do not
# move here, so `ld_all` is built by the caller and passed in. However, `log_p`
# changes whenever `w` does and is rebuilt after each weight update. EM does not
# as it operates on state directly.
#
# The EM sweeps at the foot of the file change both the weights and the location
# of the atoms, but they do not change the size of the support.
#
# A step rule is a map from a step size to a weight vector, and finding that step
# size means minimising KL, so the step layer evaluates KL tens of times per
# step. However, the corrective solve does not evaluate KL at all. It follows a
# multiplicative update that amounts to an MM step that is monotone by
# construction, so there is nothing to line search and no need to confirm the
# objective fell. It stops on the KKT residual, which the sweep computes anyway.
# `correct_weights` does evaluate KL once per iteration, but only to record it
# in the trace. It takes exactly the same path as `solve_weights` otherwise.
#
# Throughout, the candidate atom is carried as a column of `ld_all` at index
# `new_idx`, entering with weight zero, whether or not the step ends up using
# it. `new_idx` is the position the candidate will occupy in the state's flat
# ordering, so no re-indexing is needed when the weights are written back.

# --- Oracles ------------------------------------------------------------------

# Rd topic for the Li--Barron (greedy) and Frank--Wolfe (linear) oracles.

#' KL Minimisation Oracles
#'
#' The RIPr fit minimises \eqn{KL(Q \| P)}{KL(Q || P)} over mixtures `P` whose
#' components lie in the null. An *oracle* proposes the single component to
#' bring in next, and [fw_step()] and [lb_step()] differ only in how far ahead
#' they look.
#'
#' Write \eqn{P_i}{P_i} for the mixture at the current iterate and
#' \eqn{G(\theta) = E_Q[P_\theta / P_i] = E_\theta[Q / P_i]}{G(theta) = E_Q[P_theta / P_i] = E_theta[Q / P_i]},
#' the same integral written two ways.
#'
#'
#' The Frank--Wolfe \insertCite{Jaggi2013}{ripr} linear oracle [fw_step()] asks
#' how fast KL falls if an infinitesimal amount of mass moves towards
#' \eqn{P_\theta}{P_theta}:
#'
#' \deqn{\left.\frac{d}{d\epsilon} KL\!\left(Q \,\|\, (1-\epsilon) P_i + \epsilon P_\theta\right)\right|_{\epsilon=0} = 1 - G(\theta).}{d/d(eps) KL(Q || (1 - eps) P_i + eps P_theta) =  1 - G(theta) at eps=0.}
#'
#' In the language of safe, anytime-valid inference, this is the location in
#' the null that the current likelihood ratio \eqn{R_i=Q/P_i}{R_i=Q/P_i} "fails
#' to be an e-variable the most", i.e. the point that maximises the expected
#' value of \eqn{R_i}{R_i}.
#'
#' The attained maximum is
#' \eqn{\sup_{\theta \in \Theta_0} G(\theta) = 1 + \mathrm{gap}}{sup G = 1 + gap},
#' which bounds the suboptimality (in KL) via Frank--Wolfe duality, a well known
#' result in the convex optimisation literature \insertCite{Jaggi2013}{ripr}.
#'
#' The implementation of [fw_step()] yields a lower-bound on the gap at no extra
#' cost, since it approximates the Frank--Wolfe oracle. It is to be treated only
#' as a lower bound on the true duality gap, since finding the true supremum is
#' non-convex in general.
#'
#'
#' The Li--Barron \insertCite{LiBarron1999}{ripr} greedy oracle [lb_step()] does
#' not linearise the objective. It asks which atom minimises KL after the weight
#' selection inner optimisation:
#'
#' \deqn{\theta^* = \arg\min_{\theta \in \Theta_0} KL\!\left(Q \,\|\, (1 - w(\theta)) P_i + w(\theta) P_\theta\right),}{theta* = argmin_theta KL(Q || (1 - w(theta)) P_i + w(theta) P_theta),}
#'
#' where \eqn{w(\theta)}{w(theta)} depends on the weight update procedure.
#'
#' The weight selection may be an inner optimisation nested inside every
#' objective evaluation, so evaluating each candidate may cost a line search
#' or, with a corrective solve, a full re-optimisation of every weight.
#' There is no duality result that bounds suboptimality with this oracle, but
#' we cab estimate the Frank--Wolfe gap via a second optimisation sweep with the
#' Frank--Wolfe oracle by setting `record_gap = TRUE` in [lb_step()], without
#' adding the corresponding atom.
#'
#' # Which to use
#'
#' It is advisable to use [fw_step()] unless there is a specific reason not to.
#' It is cheaper per iteration and has similar if not the same convergence
#' guarantees, though no convergence guarantee applies in either case if the
#' oracle can not be trusted---as is the case here.
#'
#' Expect the two to disagree sharply on the reported gap even when they agree
#' on KL, and do not read that as one fitting better. The Frank--Wolfe oracle
#' places its atom exactly at the worst-case \eqn{\theta}{theta}, so the mixture
#' absorbs the very point defining the gap; the Li--Barron oracle places its
#' atom wherever the post-step KL is smallest, leaving that point untouched.
#' Two fits can sit at near-identical KL with gaps orders of magnitude apart.
#'
#'
#' @name oracles
#' @references
#'   \insertAllCited{}
#' @seealso [fw_step()], [lb_step()]
NULL


#' The linear oracle: maximise `G(theta)` over the null
#'
#' See [oracles] for the user-facing account. Implementation notes only here.
#'
#' `value` returns \eqn{G(\theta)}{G(theta)} and `maximise_over()` maximises it,
#' which turns out to be descent on KL (see [oracles] for why). The gradient
#' is \eqn{E_\theta[(Q/P_i)\, s_\theta]}{E_theta[(Q/P_i) s_theta]}, by
#' differentiating under the integral. `value_batch` scores a whole matrix of
#' candidates in one reduction, which is what the multi-start seeding uses.
#' @keywords internal
#' @noRd
linear_oracle <- function(state, log_p, ld) {
  engine <- state@engine
  family <- engine@family
  w <- exp(engine@log_w)

  objective(
    value = function(theta) {
      exp(log_expect_q(engine, as.vector(ld(matrix(theta, ncol = 1L))) - log_p))
    },
    grad = function(theta) {
      ratio <- exp(as.vector(ld(matrix(theta, ncol = 1L))) - log_p)
      as.vector(crossprod(score(family, theta, engine@nodes), w * ratio))
    },
    value_batch = function(theta_mat) {
      exp(col_logsumexp(ld(theta_mat) - log_p + engine@log_w))
    }
  )
}


#' Estimate the Frank--Wolfe duality gap
#'
#' Maximises (locally) the linear oracle over the null.
#'
#' @keywords internal
#' @noRd
linear_gap <- function(state, log_p, ld, seeds) {
  ctl <- state@control
  obj <- linear_oracle(state, log_p, ld)
  values <- vapply(
    state@null@subnulls,
    function(s) {
      maximise_over(
        s,
        obj,
        seeds = seeds,
        n_seeds = ctl$n_seeds,
        n_restarts = ctl$n_restarts
      )$value
    },
    numeric(1)
  )
  max(values) - 1
}

#' Directions the Li--Barron inner optimisation may use
#'
#' `"away"` is dropped. `path_away` never touches the candidate column, so an
#' away path yields the same mixture for every candidate, and the gradient is
#' exactly zero. Leaving it in lays a plateau over the objective and strands
#' at its seed every restart landing on it. Dropping it costs nothing: the
#' outer step still weighs the away direction against the candidate the oracle
#' returns, and when away wins the candidate is discarded regardless of what it
#' was.
#' @keywords internal
#' @noRd
inner_directions <- function(directions) {
  d <- setdiff(directions, "away")
  if (length(d)) d else "forward"
}


#' What a step towards `theta` would do, without doing it
#'
#' Returns a function of `theta` reporting the weights, mixture and KL that a
#' step towards it would produce. Nothing is written to the state, so an oracle
#' can score candidates with the same machinery that later takes the step.
#'
#' The candidate enters at `at` with weight zero. `at = NULL` puts it last,
#' which is what an oracle wants, since it never writes back. A caller that does
#' intend to commit passes `insert_index()`, so the weights come back in the
#' state's own flat ordering and need no permutation.
#'
#' The returned closure memoises its last call. `optim()` asks for the value and
#' the gradient at the same point through separate slots of [objective()], and
#' each evaluation here costs a line search -- or, under `correct`, an entire
#' weight solve.
#'
#' @param at Index the candidate should occupy, or `NULL` for last.
#' @param correct Re-solve every weight after the step, as Li and Barron allow.
#' @return A function of `theta` giving `weights`, `log_p`, `kl`, `gamma`,
#'   `direction`, `uses_candidate`, plus `ld_new` and `new_idx` for the caller.
#' @keywords internal
#' @noRd
plan_step <- function(
  state,
  log_p,
  ld,
  directions = "forward",
  size = "line-search",
  gamma_fixed = NULL,
  correct = FALSE,
  at = NULL
) {
  engine <- state@engine
  ctl <- state@control
  ld_atoms <- ld(flat_atoms(state))
  w_now <- flat_weights(state)

  new_idx <- if (is.null(at)) length(w_now) + 1L else at
  w <- append(w_now, 0, after = new_idx - 1L)

  last_theta <- NULL
  last <- NULL

  function(theta) {
    # Caches previous theta so that nonlinear_oracle doesn't recompute twice
    # for gradient + objective.
    if (!is.null(last_theta) && identical(theta, last_theta)) {
      return(last)
    }
    ld_new <- as.vector(ld(matrix(theta, ncol = 1L)))
    ld_all <- insert_col(ld_atoms, ld_new, new_idx)

    res <- apply_step(
      ld_all,
      w,
      new_idx,
      log_p,
      engine,
      directions = directions,
      size = size,
      gamma_fixed = gamma_fixed
    )
    if (correct) {
      res$weights <- solve_weights(
        ld_all,
        res$weights,
        engine,
        tol = ctl$lb_fc_tol,
        max_iter = ctl$lb_fc_max_iter
      )
      res$log_p <- mixture_log_p(ld_all, res$weights)
      res$kl <- expect_q(engine, engine@log_q - res$log_p)
    }
    res$ld_new <- ld_new
    res$new_idx <- new_idx

    last_theta <<- theta
    last <<- res
    res
  }
}


#' Gradient of the Li--Barron objective
#'
#' By the envelope theorem the derivative with respect to the inner variables
#' vanishes at their optimum, so the weight can be held fixed and only the
#' candidate's own dependence on \eqn{\theta}{theta} differentiated:
#' \deqn{\partial_\theta E_Q[\log P] = E_Q\!\left[\frac{w(\theta) P_\theta}{P}\, s_\theta\right].}{d/dtheta E_Q[log P] = E_Q[(w(theta) P_theta / P) s_theta].}
#' No differentiation through the line search or the weight solve is needed.
#'
#' Zero when the step left the candidate unweighted: it is then absent from the
#' mixture, so moving it changes nothing.
#' @param planned One result from a `plan_step` closure.
#' @keywords internal
#' @noRd
lb_gradient <- function(state, theta, planned) {
  w_new <- planned$weights[planned$new_idx]
  if (w_new <= 0) {
    return(numeric(length(theta)))
  }
  engine <- state@engine
  share <- w_new * exp(planned$ld_new - planned$log_p)
  as.vector(crossprod(
    score(engine@family, theta, engine@nodes),
    exp(engine@log_w) * share
  ))
}

#' The Li--Barron nonlinear oracle
#'
#' See [oracles] for the user-facing account. Implementation notes only here.
#'
#' Li and Barron score a candidate by the KL it leaves behind *after* its weight
#' has been chosen, so the step rule is this oracle's inner optimisation and
#' `directions`, `size` and `correct` are consumed here rather than by the
#' caller. `value` returns `-kl`, since maximising that is minimising KL, which
#' is the definition in [oracles].
#'
#' `value_batch` scores every seed with a full step rather than ranking them by
#' a cheap proxy. That is expensive -- `n_seeds` line searches per oracle call,
#' or `n_seeds` weight solves under `correct` -- and deliberate: this oracle
#' exists to be compared against [fw_step()], and a baseline weakened by an
#' approximate seeding heuristic would not be worth the comparison.
#' @keywords internal
#' @noRd
nonlinear_oracle <- function(
  state,
  log_p,
  ld,
  directions = "forward",
  size = "line-search",
  gamma_fixed = NULL,
  correct = FALSE
) {
  after <- plan_step(
    state,
    log_p,
    ld,
    directions = inner_directions(directions),
    size = size,
    gamma_fixed = gamma_fixed,
    correct = correct
  )

  objective(
    value = function(theta) -after(theta)$kl,
    grad = function(theta) lb_gradient(state, theta, after(theta)),
    value_batch = function(theta_mat) {
      -vapply(
        seq_len(ncol(theta_mat)),
        \(i) after(theta_mat[, i])$kl,
        numeric(1)
      )
    }
  )
}


# --- Step rules ---------------------------------------------------------------

#' Log density of a mixture with weights `w` over the atoms of `ld_all`
#' @keywords internal
#' @noRd
mixture_log_p <- function(ld_all, w) {
  row_logsumexp(add_by_col(ld_all, log(pmax(w, 0))))
}


#' The active atom with the smallest `G`, or `NULL` if fewer than two are active
#'
#' Fewer than two leaves nothing to move away from: emptying the only active
#' atom would empty the mixture.
#' @keywords internal
#' @noRd
worst_atom <- function(ld_all, w, log_p, engine) {
  active <- which(w > 0)
  if (length(active) < 2L) {
    return(NULL)
  }
  g <- exp(col_logsumexp(ld_all[, active, drop = FALSE] - log_p + engine@log_w))
  active[which.min(g)]
}


#' A one-parameter path through weight space
#'
#' `gamma -> w`, plus the largest admissible `gamma` and the log density along
#' the way. One shape for every direction, so one line search serves all of them.
#' @keywords internal
#' @noRd
weight_path <- function(direction, gamma_max, w_of, log_p_at) {
  list(
    direction = direction,
    gamma_max = gamma_max,
    w_of = w_of,
    log_p_at = log_p_at
  )
}


#' Forward: transfer mass to the candidate from everything else
#'
#' The ordinary standard direction,
#' \eqn{w \leftarrow (1-\gamma) w + \gamma e_{new}}{w <- (1 - gamma) w + gamma e_new}.
#' @keywords internal
#' @noRd
path_forward <- function(ld_all, w, new_idx, log_p) {
  ld_new <- ld_all[, new_idx]
  weight_path(
    direction = "forward",
    gamma_max = 1,
    w_of = function(gamma) {
      out <- (1 - gamma) * w
      out[new_idx] <- gamma
      out
    },
    # Two columns rather than C + 1: this direction is a convex combination of
    # the current mixture with one new atom, so `log_p` is reused. Worth the
    # special case, since the line search calls this tens of times per step.
    log_p_at = function(gamma) {
      if (gamma <= 0) {
        return(log_p)
      }
      if (gamma >= 1) {
        return(ld_new)
      }
      row_logsumexp(cbind(log_p + log1p(-gamma), ld_new + log(gamma)))
    }
  )
}


#' Pairwise: transfer mass to the candidate from the worst active atom
#'
#' Lacoste-Julien & Jaggi's pairwise step. Capped at the worst atom's weight,
#' so at the cap that atom empties.
#'
#' With one active atom `worst` is that atom and this reduces to
#' `path_forward` exactly.
#' @references
#'   \insertRef{LacosteJulienJaggi2015}{ripr}
#' @keywords internal
#' @noRd
path_pairwise <- function(ld_all, w, new_idx, log_p, worst) {
  # No two-column shortcut here: mass leaves one atom and arrives at another, so
  # the mixture is rebuilt.
  move <- function(gamma) {
    out <- w
    out[worst] <- out[worst] - gamma
    out[new_idx] <- gamma
    out
  }
  weight_path(
    direction = "pairwise",
    gamma_max = w[worst],
    w_of = move,
    log_p_at = function(gamma) mixture_log_p(ld_all, move(gamma))
  )
}


#' Away: remove mass from the worst active atom, spread evenly across all others
#'
#' \eqn{w \leftarrow (1+\gamma) w - \gamma e_v}{w <- (1 + gamma) w - gamma e_v}.
#' The candidate goes unused, which is why this is only ever offered alongside
#' another direction.
#' @references
#'   \insertRef{LacosteJulienJaggi2015}{ripr}
#' @keywords internal
#' @noRd
path_away <- function(ld_all, w, new_idx, log_p, worst) {
  # Cap is w_v/(1 - w_v): beyond it the worst atom's weight would go negative.
  w_of <- function(gamma) {
    out <- (1 + gamma) * w
    out[worst] <- out[worst] - gamma
    out
  }
  weight_path(
    direction = "away",
    gamma_max = w[worst] / (1 - w[worst]),
    w_of = w_of,
    log_p_at = function(gamma) mixture_log_p(ld_all, w_of(gamma))
  )
}


#' The directions on offer this step
#'
#' One path per requested direction, unavailable ones dropped. `apply_step`
#' line searches each and keeps the lowest KL, so the set names a choice rather
#' than a single move: `"forward"` alone is vanilla Frank--Wolfe,
#' `c("forward", "away")` is away-step Frank--Wolfe, `"pairwise"` is pairwise
#' Frank--Wolfe.
#'
#' Only `"away"` can be unavailable, and only below two active atoms. Asking for
#' it alone in that state is an error rather than a silent fallback -- a
#' fallback is what previously made `"away"` behave as pairwise while its
#' documentation claimed forward.
#' @return A non-empty list of `weight_path`s.
#' @keywords internal
#' @noRd
step_paths <- function(directions, ld_all, w, new_idx, log_p, engine) {
  worst <- worst_atom(ld_all, w, log_p, engine)
  paths <- lapply(directions, function(d) {
    switch(
      d,
      forward = path_forward(ld_all, w, new_idx, log_p),
      # `worst` is NULL below two active atoms; pairwise then uses the single
      # active atom and coincides with forward. See `path_pairwise`.
      pairwise = path_pairwise(
        ld_all,
        w,
        new_idx,
        log_p,
        if (is.null(worst)) which(w > 0)[1L] else worst
      ),
      away = if (!is.null(worst)) path_away(ld_all, w, new_idx, log_p, worst)
    )
  })
  paths <- Filter(Negate(is.null), paths)
  if (!length(paths)) {
    stop(
      "`away` needs a second active atom to move mass to. ",
      "Add `\"forward\"` to `directions`, or take another step first.",
      call. = FALSE
    )
  }
  paths
}


#' The open-loop step size, for `size = "fixed"`
#'
#' The alternative to `line_search`. The only rule here that can increase KL,
#' kept to reproduce open-loop results in the literature rather than as a good
#' default.
#'
#' `k` must count only steps that took a Frank--Wolfe direction. EM sweeps
#' between two steps must not advance the schedule, or it is no longer the
#' sequence the published rates are proved for.
#' @references
#'   \insertRef{Jaggi2013}{ripr}
#'
#'   \insertRef{LiBarron1999}{ripr}
#' @keywords internal
#' @noRd
schedule_gamma <- function(k) {
  # Jaggi's 2/(k+2) from k = 0 and Li--Barron's 2/(k+1) from k = 1 are the same
  # sequence; the lineages differ in the oracle, not the schedule.
  2 / (k + 1)
}


#' Minimise KL along a step path
#'
#' `gamma = 0` is always in range, so no step can increase KL.
#' @keywords internal
#' @noRd
line_search <- function(log_p_at, gamma_max, engine) {
  if (!is.finite(gamma_max) || gamma_max <= 0) {
    return(0)
  }
  stats::optimize(
    \(gamma) expect_q(engine, engine@log_q - log_p_at(gamma)),
    interval = c(0, gamma_max),
    tol = 1e-12
  )$minimum
}


#' Take one step, without a state and without a trace
#'
#' Compares every path on offer from `step_paths` and keeps the one with
#' lowest KL. `size = "fixed"` does not search, taking `gamma_fixed` capped
#' at the path's own maximum -- the cap matters, since pairwise and away cap
#' below 1 and an uncapped schedule value would leave the simplex.
#' @return `list(weights, log_p, kl, gamma, direction, uses_candidate)`, with
#'   `weights` of length `C + 1`.
#' @keywords internal
#' @noRd
apply_step <- function(
  ld_all,
  w,
  new_idx,
  log_p,
  engine,
  directions = "forward",
  size = "line-search",
  gamma_fixed = NULL
) {
  take <- function(path) {
    gamma <- if (size == "fixed") {
      min(gamma_fixed, path$gamma_max)
    } else {
      line_search(path$log_p_at, path$gamma_max, engine)
    }
    stepped <- path$log_p_at(gamma)
    weights <- pmax(path$w_of(gamma), 0)
    list(
      weights = weights,
      log_p = stepped,
      kl = expect_q(engine, engine@log_q - stepped),
      gamma = gamma,
      direction = path$direction,
      # Derived, not declared: the candidate is used exactly when it ends with
      # weight. False for away, which never touches it, and for any search that
      # put nothing there.
      uses_candidate = weights[new_idx] > 0
    )
  }

  taken <- lapply(
    step_paths(directions, ld_all, w, new_idx, log_p, engine),
    take
  )
  taken[[which.min(vapply(taken, \(r) r$kl, numeric(1)))]]
}


# --- Fully corrective weights -------------------------------------------------

#' One multiplicative sweep on the mixture weights
#'
#' \eqn{w_c \leftarrow w_c G(\theta_c)}{w_c <- w_c G(theta_c)}: the exact M-step
#' for the weights with the atoms fixed, and the MM algorithm for the convex
#' problem of minimising KL over the simplex. Monotone by construction, so
#' nothing needs checking afterwards.
#'
#' `residual` is \eqn{\max_c G(\theta_c) - 1}{max_c G(theta_c) - 1}, measured
#' *before* the sweep. It is the Frank--Wolfe gap of the restricted problem over
#' the current support, so it is a bona fide upper-bound for how much KL is
#' still available by optimising weights.
#' @keywords internal
#' @noRd
weight_sweep <- function(ld_all, w, log_p, engine) {
  g <- exp(col_logsumexp(ld_all - log_p + engine@log_w))
  # `sum_c w_c G_c = 1` identically, so the result needs no renormalisation.
  list(weights = w * g, residual = max(g) - 1)
}


#' Re-optimise every weight over the current atoms
#'
#' Stops on the KKT residual, not on the change in KL. A `dKL` rule cannot tell
#' convergence from crawling.
#'
#' It is likely that `max_iter` will be reached. The convergence rate depends on
#' how close the atoms are to one another, and we often place new atoms near old
#' ones, so it degrades as the fit proceeds. This is a budget, not a correctness
#' condition: fully-corrective Frank--Wolfe re-corrects on the next step either
#' way.
#' @param tol Stop when the residual falls below this.
#' @param max_iter Cap on sweeps.
#' @keywords internal
#' @noRd
solve_weights <- function(ld_all, w, engine, tol, max_iter) {
  log_p <- mixture_log_p(ld_all, w)
  for (k in seq_len(max_iter)) {
    sweep <- weight_sweep(ld_all, w, log_p, engine)
    if (sweep$residual < tol) {
      break
    }
    w <- sweep$weights
    # The atoms do not move, so only the mixture is rebuilt, not `ld_all`.
    log_p <- mixture_log_p(ld_all, w)
  }
  w
}


# --- Updating state -----------------------------------------------------------

#' Mix the new atom into the current iterate
#'
#' The stateful half: assemble the candidate into the flat arithmetic, hand off
#' to `apply_step`, write the result back and record it. Every rule here
#' changes only *how fast* the iterate moves, never what it converges to.
#'
#' The candidate is inserted at the index `add_atom` would place it, so the
#' weights come back in the state's own flat ordering and need no permutation.
#' @keywords internal
#' @noRd
step_update <- function(
  state,
  theta_new,
  subnull_index,
  gap,
  oracle_value,
  log_p,
  ld,
  directions = "forward",
  size = "line-search",
  gamma_fixed = NULL
) {
  engine <- state@engine

  new_idx <- insert_index(state, subnull_index)
  ld_new <- as.vector(ld(matrix(theta_new, ncol = 1L)))
  ld_all <- insert_col(ld(flat_atoms(state)), ld_new, new_idx)
  w <- append(flat_weights(state), 0, after = new_idx - 1L)

  res <- apply_step(
    ld_all,
    w,
    new_idx,
    log_p,
    engine,
    directions = directions,
    size = size,
    gamma_fixed = gamma_fixed
  )

  state <- if (res$uses_candidate) {
    add_atom(state, theta_new, subnull_index, res$weights)
  } else {
    set_weights(state, res$weights[-new_idx])
  }

  state <- record(
    state,
    phase = "step",
    inner = 0L,
    kl = res$kl,
    gap = gap,
    oracle_value = oracle_value,
    subnull = if (res$uses_candidate) subnull_index else NA_integer_,
    step_size = res$gamma,
    direction = res$direction,
    ld = ld
  )

  state
}


# --- EM -----------------------------------------------------------------------

#' One EM sweep: weights, then atoms: weights, then atoms
#'
#' Weights first, because the atom M-step conditions on the responsibilities and
#' those are sharper once the weights have been updated.
#'
#' Always both halves. Moving only the weights is the corrective step, reached
#' through `solve_weights` with its own convergence test rather than by
#' running this a fixed number of times.
#' @keywords internal
#' @noRd
em_sweep <- function(state, ld) {
  wt <- exp(state@engine@log_w) * em_responsibilities(state, ld)
  em_atom_step(em_weight_step(state, wt), ld, wt)
}


#' Responsibility of each atom for each node
#'
#' `r_ic = w_c p_c(x_i) / P(x_i)`, so rows sum to 1.
#' @keywords internal
#' @noRd
em_responsibilities <- function(state, ld) {
  log_comp <- add_by_col(ld(flat_atoms(state)), log(flat_weights(state)))
  exp(log_comp - row_logsumexp(log_comp))
}


#' The M-step for the weights
#'
#' \eqn{w_c \leftarrow E_Q[r_c]}{w_c <- E_Q[r_c]}, which is the column sums of
#' the weighted responsibilities. Equal to the multiplicative update
#' `weight_sweep` performs, reached from the other direction:
#' \eqn{E_Q[r_c] = w_c G(\theta_c)}{E_Q[r_c] = w_c G(theta_c)}.
#'
#' @param wt `(M, C)` responsibilities scaled by the quadrature weights.
#' @keywords internal
#' @noRd
em_weight_step <- function(state, wt) {
  new_w <- colSums(wt)
  # Sums to 1 already, since the rows of `wt` sum to the quadrature weights.
  set_weights(state, new_w / sum(new_w))
}


#' The M-step for the atoms
#'
#' Each atom maximises its own responsibility-weighted log-likelihood over its
#' own subnull, so an atom cannot migrate between subnulls however the
#' responsibilities fall.
#'
#' Local only: seeded at the current atom with no random starts. Exploration is
#' the oracle's job, and a global search here would let atoms teleport between
#' sweeps.
#'
#' @inheritParams em_weight_step
#' @keywords internal
#' @noRd
em_atom_step <- function(state, ld, wt) {
  engine <- state@engine
  family <- engine@family
  atoms_flat <- flat_atoms(state)
  if (ncol(atoms_flat) == 0L) {
    return(state)
  }
  idx <- flat_subnull(state)

  moved <- vapply(
    seq_len(ncol(atoms_flat)),
    function(c_i) {
      w_c <- wt[, c_i]
      obj <- objective(
        value = function(theta) {
          sum(w_c * as.vector(ld(matrix(theta, ncol = 1L))))
        },
        grad = function(theta) {
          as.vector(crossprod(score(family, theta, engine@nodes), w_c))
        },
        value_batch = function(theta_mat) {
          as.vector(crossprod(ld(theta_mat), w_c))
        }
      )
      maximise_over(
        state@null@subnulls[[idx[c_i]]],
        obj,
        seeds = atoms_flat[, c_i, drop = FALSE],
        n_seeds = 0L,
        n_restarts = 1L
      )$theta
    },
    numeric(nrow(atoms_flat))
  )
  state@atoms <- unflatten_atoms(state, matrix(moved, nrow = nrow(atoms_flat)))
  state
}


# --- Step verbs ---------------------------------------------------------------

#' Search every subnull and keep the best candidate
#'
#' Returns the subnull index, its maximiser, and the attained value. `seeds` is
#' the current atoms: the identity `sum_c w_c G(theta_c) = 1` forces
#' `max_c G(theta_c) >= 1`, so including them stops the search reporting a
#' maximum below one.
#' @keywords internal
#' @noRd
search_null <- function(state, obj) {
  ctl <- state@control
  found <- lapply(
    state@null@subnulls,
    \(s) {
      maximise_over(
        s,
        obj,
        seeds = flat_atoms(state),
        n_seeds = ctl$n_seeds,
        n_restarts = ctl$n_restarts
      )
    }
  )
  best <- which.max(vapply(found, \(f) f$value, numeric(1)))
  c(found[[best]], list(subnull = best))
}


#' Write a planned step back to the state
#'
#' The candidate becomes an atom only if the step gave it weight; `"away"` never
#' does, and neither does a search that settled on `gamma = 0`. Incumbent atoms
#' driven towards zero stay, since only an oracle can grow the support back.
#' @keywords internal
#' @noRd
commit_step <- function(state, theta, subnull, planned) {
  if (planned$uses_candidate) {
    add_atom(state, theta, subnull, planned$weights)
  } else {
    set_weights(state, planned$weights[-planned$new_idx])
  }
}


# --- Support identification ---------------------------------------------------

#' Log density with atom `c_i` removed and the rest renormalised
#'
#' A rank-one downdate, `O(M)` rather than the `O(MC)` of rebuilding.
#' @keywords internal
#' @noRd
log_p_without <- function(log_p, ld_c, w_c) {
  # `w_c p_c <= P` holds exactly but can fail by an ulp when one atom carries
  # nearly all the mass, sending `log1p` to NaN rather than -Inf.
  share <- pmin(w_c * exp(ld_c - log_p), 1)
  log_p + log1p(-share) - log1p(-min(w_c, 1 - .Machine$double.eps))
}


#' Zero every atom whose removal does not increase KL
#'
#' One pass in ascending weight order, each removal applied before the next is
#' tested. Lightest first because those are the likeliest removals, and applying
#' as we go means the later tests see the mass already redistributed.
#'
#' Exact where a weight threshold is arbitrary: an atom carrying weight far
#' above any sensible floor can still be strictly better removed, and the KKT
#' screen \eqn{G < 1}{G < 1} over-identifies when the iterate has not converged.
#' The screen only shortlists; the KL comparison decides.
#'
#' **Only safe once the atoms have stopped moving.** Nothing can restore a
#' zeroed atom -- the multiplicative update maps zero to zero, and an oracle
#' will not re-propose a point whose `G` is below one. Run mid-fit it ratchets
#' the support down and raises KL.
#' @keywords internal
#' @noRd
identify_support <- function(ld_all, w, engine) {
  log_p <- mixture_log_p(ld_all, w)
  kl <- expect_q(engine, engine@log_q - log_p)

  for (c_i in order(w)) {
    active <- which(w > 0)
    if (w[c_i] <= 0 || length(active) < 2L) {
      next
    }
    trial <- log_p_without(log_p, ld_all[, c_i], w[c_i])
    kl_trial <- expect_q(engine, engine@log_q - trial)
    if (kl_trial <= kl) {
      w[c_i] <- 0
      w <- w / sum(w)
      log_p <- trial
      kl <- kl_trial
    }
  }
  w
}
