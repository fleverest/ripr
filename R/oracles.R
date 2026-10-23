# Oracles: where the next atom should go (see ?oracles).

#' KL Minimisation Oracles
#'
#' The RIPr fit minimises \eqn{KL(Q \| P)}{KL(Q || P)} over mixtures `P` whose
#' components lie in the null. An *oracle* proposes the next component to
#' bring in; [fw_step()] and [lb_step()] differ only in how far ahead they look.
#'
#' Write \eqn{P_i}{P_i} for the current mixture and
#' \eqn{G(\theta) = E_Q[P_\theta / P_i] = E_\theta[Q / P_i]}{G(theta) = E_Q[P_theta / P_i] = E_theta[Q / P_i]}.
#'
#' The Frank--Wolfe \insertCite{Jaggi2013}{ripr} linear oracle [fw_step()] asks
#' how fast KL falls as an infinitesimal amount of mass moves towards
#' \eqn{P_\theta}{P_theta}:
#'
#' \deqn{\left.\frac{d}{d\epsilon} KL\!\left(Q \,\|\, (1-\epsilon) P_i + \epsilon P_\theta\right)\right|_{\epsilon=0} = 1 - G(\theta).}{d/d(eps) KL(Q || (1 - eps) P_i + eps P_theta) =  1 - G(theta) at eps=0.}
#'
#' Equivalently, it finds where in the null the likelihood ratio
#' \eqn{R_i=Q/P_i}{R_i=Q/P_i} most fails to be an e-variable. The attained
#' \eqn{\sup_{\theta \in \Theta_0} G(\theta) = 1 + \mathrm{gap}}{sup G = 1 + gap}
#' bounds the KL suboptimality by Frank--Wolfe duality. Since the search is
#' non-convex, the reported gap is only a lower bound on the true gap. It is
#' the gap of the mixture the step started *from*, so it fills the previous
#' row's `gap_after`.
#'
#' The Li--Barron \insertCite{LiBarron1999}{ripr} greedy oracle [lb_step()]
#' does not linearise; it picks the atom minimising KL after choosing its
#' weight:
#'
#' \deqn{\theta^* = \arg\min_{\theta \in \Theta_0} KL\!\left(Q \,\|\, (1 - w(\theta)) P_i + w(\theta) P_\theta\right),}{theta* = argmin_theta KL(Q || (1 - w(theta)) P_i + w(theta) P_theta),}
#'
#' where \eqn{w(\theta)}{w(theta)} comes from the weight update. Each
#' candidate may therefore cost a line search, or a full weight solve under
#' correction. It has no duality bound; `record_gap = TRUE` estimates the
#' Frank--Wolfe gap with an extra linear-oracle sweep, without adding its atom.
#'
#' # What a trace row records
#'
#' Every verb writes one row per step to `state@trace`. A row spans a step, so
#' its columns split by which mixture they measure.
#'
#' \describe{
#'   \item{`step`, `phase`}{The row's index ([ripr_init()]'s row is 0) and the
#'   verb that wrote it (`"init"`, `"fw"`, `"lb"`, `"em"` or `"weight"`).
#'   Count steps of a kind with e.g. `cumsum(trace$phase == "fw")`.}
#'   \item{`oracle_value`, `oracle_theta`}{What the step's own search found at
#'   the mixture it stepped *from*. For [fw_step()], `oracle_theta` is the
#'   linear oracle's maximiser and `oracle_value` is `NA` (its value is the
#'   previous row's `gap_after + 1`). For [lb_step()] `oracle_value` is the
#'   Li--Barron objective, not a gap; for [weight_step()] it is the pre-sweep
#'   support gap plus one, with no `oracle_theta`. When `part` is not `NA`, the
#'   atom went to `oracle_theta`.}
#'   \item{`kl`, `gap_after`, `gap_after_theta`, `gap_after_part`}{The mixture
#'   the row *produced*, and the linear oracle's gap, maximiser and part for
#'   it. Filled by `record_gap = TRUE` or the next [fw_step()]; otherwise `NA`.}
#'   \item{`gap_after_elapsed`}{Seconds the search that filled `gap_after`
#'   took.}
#'   \item{`part`, `step_size`, `direction`, `support_size`, `max_weight`}{
#'   Where an oracle step put its atom (`NA` if none), the step length and
#'   direction, and the size and heaviest weight of the produced mixture.}
#'   \item{`elapsed`}{Seconds of the step rule's work, excluding diagnostics
#'   (`record_gap`, `snapshot`, `until`). An [fw_step()] row includes the
#'   linear-oracle search it consumed, even if cached. Per step; use
#'   `cumsum(trace$elapsed)` for the total.}
#' }
#'
#' Both `theta` columns are list columns holding the family's parameter, with
#' `NA` for rows that recorded none.
#'
#' `state@oracle` caches the linear oracle's result for the current mixture
#' only: a step that changes the mixture clears it and `record_gap = TRUE`
#' sets it. [fw_step()] uses the cache when present, so an fw run pays one
#' search per step. The final row's `gap_after` stays `NA` unless
#' `record_gap = TRUE` or `until` stopped the run; `ripr_finish(record_gap =
#' TRUE)` measures the mixture it returns.
#'
#' # Which to use
#'
#' Prefer [fw_step()]: it is cheaper per iteration with similar guarantees
#' (none of which strictly hold, since the oracle is heuristic).
#'
#' Expect the two to disagree sharply on the gap even if the KL is similar.
#' Frank--Wolfe places its atom exactly at the worst-case \eqn{\theta}{theta},
#' absorbing the point that defines the gap; Li--Barron places it wherever
#' post-step KL is smallest, leaving that point untouched.
#'
#' @name oracles
#' @references
#'   \insertAllCited{}
#' @seealso [fw_step()], [lb_step()]
NULL


#' The linear oracle: maximise `G(theta)` over the null (see [oracles])
#'
#' Gradient is \eqn{E_Q[\nabla_\theta p_\theta / P_i]}{E_Q[grad p_theta / P_i]},
#' differentiating under the integral; see `compile_density_grad()`.
#' @keywords internal
#' @noRd
linear_oracle <- function(state, log_p, ld) {
  engine <- state@engine
  dg <- compile_engine_grad(engine)
  log_m <- engine@log_w - log_p

  objective(
    value = function(theta) {
      exp(log_expect_q(engine, as.vector(ld(matrix(theta, nrow = 1L))) - log_p))
    },
    grad = function(theta) as.vector(dg(theta, log_m)),
    value_batch = function(theta_mat) atom_g(ld(theta_mat), log_p, engine)
  )
}


#' What a step towards `theta` would do, without doing it
#'
#' Returns a memoised function of `theta` giving `weights`, `log_p`, `kl`,
#' `gamma`, `direction` and `uses_candidate`, so an oracle can score
#' candidates with the machinery that later takes the step.
#'
#' `correct` re-solves every weight after the step, which may zero the
#' candidate or an incumbent.
#' @keywords internal
#' @noRd
plan_step <- function(
  state,
  log_p,
  ld,
  directions = "forward",
  size = "line-search",
  gamma_fixed = NULL,
  correct = FALSE
) {
  engine <- state@engine
  ctl <- state@control
  ld_atoms <- ld(state@mixing@atoms)
  w <- c(state@mixing@weights, 0)

  memoise_last(function(theta) {
    ld_new <- as.vector(ld(matrix(theta, nrow = 1L)))
    ld_all <- cbind(ld_atoms, ld_new, deparse.level = 0)

    res <- apply_step(
      ld_all,
      w,
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
        tol = ctl$fc_tol,
        max_iter = ctl$fc_max_iter
      )
      res$log_p <- mixture_log_p(ld_all, res$weights)
      res$kl <- kl_at(engine, res$log_p)
      res$uses_candidate <- res$weights[length(res$weights)] > 0
    }
    res
  })
}


#' Gradient of the Li--Barron objective
#'
#' By the envelope theorem the weight can be held fixed:
#' \deqn{\partial_\theta E_Q[\log P] = w(\theta)\, E_Q\!\left[\frac{\nabla_\theta P_\theta}{P}\right],}{d/dtheta E_Q[log P] = w(theta) E_Q[grad P_theta / P],}
#' so no differentiation through the line search or weight solve is needed.
#' Zero when the step left the candidate unweighted. `dg` is the engine's
#' `compile_engine_grad()`.
#' @keywords internal
#' @noRd
lb_gradient <- function(state, theta, planned, dg) {
  w_new <- utils::tail(planned$weights, 1L)
  if (w_new <= 0) {
    return(numeric(length(theta)))
  }
  w_new * as.vector(dg(theta, state@engine@log_w - planned$log_p))
}

#' The Li--Barron nonlinear oracle (see [oracles]); `value` is `-kl`
#'
#' `value_batch` scores every seed with a full step (`n_seeds` line searches or
#' weight solves).
#' @keywords internal
#' @noRd
nonlinear_oracle <- function(
  state,
  log_p,
  ld,
  size = "line-search",
  gamma_fixed = NULL,
  correct = FALSE
) {
  after <- plan_step(
    state,
    log_p,
    ld,
    size = size,
    gamma_fixed = gamma_fixed,
    correct = correct
  )

  dg <- compile_engine_grad(state@engine)

  objective(
    value = function(theta) -after(theta)$kl,
    grad = function(theta) lb_gradient(state, theta, after(theta), dg),
    value_batch = function(theta_mat) {
      -vapply(
        seq_len(nrow(theta_mat)),
        \(i) after(theta_mat[i, ])$kl,
        numeric(1)
      )
    }
  )
}
