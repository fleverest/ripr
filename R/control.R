#' @include ripr-package.R
NULL

# Computational knobs that control how the answer is computed, never what is
# computed. The algorithm is defined at the call site via the step verbs,
# not here.
#
# Nothing in this file depends on anything else in the package.

#' Tuning for the RIPr optimiser
#'
#' `ripr_control` changes how well the answer is computed, never what is being
#' computed. `alternative`, `null`, and the engine define the problem; anything
#' in this object may be varied without changing the target. Operations
#' that alter the returned mixture (e.g. pruning) belong to [ripr_finish()].
#'
#' @param n_seeds Random chart seeds per oracle call.
#' @param n_restarts Number of top seeds to be refined by SLSQP per oracle
#'   call.
#' @param snapshot How often to record the fitted mixture alongside the trace.
#'   `"none"` (default) never records the state; `"step"` records the atoms and
#'   weights once per step call (e.g. [fw_step()]); `"all"` once per iteration
#'   within a call, so `times` snapshots for a verb called with `times`.
#'
#'   Note `"step"` counts calls, not iterations, so `fw_step(times = 10)` yields
#'   one snapshot while ten separate `fw_step()` calls yield ten. Using `"all"`
#'   yields ten snapshots in both cases.
#'
#'   The trace is always recorded as it costs only a row of scalars per event.
#'   Snapshots copy the whole mixture, so the memory cost grows with support
#'   size, so snapshotting is `"none"` by default.
#' @param lb_fc_tol Convergence tolerance for a corrective weight solve in
#'   the Li--Barron greedy oracle inner loop.
#' @param lb_fc_max_iter Cap on corrective weight sweeps within the Li--Barron
#'   inner loop.
#' @return A list of control settings for [ripr_init()].
#' @examples
#' ripr_control(n_seeds = 50L, snapshot = "step")
#' @export
ripr_control <- function(
  n_seeds = 200L,
  n_restarts = 25L,
  lb_fc_tol = 1e-10,
  lb_fc_max_iter = 500L,
  snapshot = c("none", "step", "all")
) {
  snapshot <- rlang::arg_match(snapshot)
  rlang::check_number_whole(n_seeds, min = 0, max = 2147483647)
  rlang::check_number_whole(n_restarts, min = 1, max = 2147483647)
  rlang::check_number_decimal(lb_fc_tol, min = 0)
  rlang::check_number_whole(lb_fc_max_iter, min = 1, max = 2147483647)
  list(
    n_seeds = as.integer(n_seeds),
    n_restarts = as.integer(n_restarts),
    lb_fc_tol = lb_fc_tol,
    lb_fc_max_iter = as.integer(lb_fc_max_iter),
    snapshot = snapshot
  )
}
