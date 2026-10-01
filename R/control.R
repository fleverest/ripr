#' @include ripr-package.R
NULL

#' Tuning for the RIPr optimiser
#'
#' @param n_seeds Random chart seeds per oracle call.
#' @param n_restarts Number of top seeds to be refined by SLSQP per oracle
#'   call.
#' @param snapshot How often to record the fitted mixture alongside the trace.
#'   `"none"` (default) never records the state; `"step"` records the atoms and
#'   weights once per step call (e.g. [fw_step()]), so `fw_step(times = 10)`
#'   yields one snapshot; `"all"` records once per iteration, yielding ten.
#'   The trace is always recorded; snapshots copy the whole mixture, so their
#'   memory cost grows with support size.
#' @param fc_tol Relative tolerance on the weights for a corrective weight
#'   solve (`correct = TRUE`, or `ripr_finish(reoptimise = TRUE)`).
#' @param fc_max_iter Cap on the objective evaluations such a solve may take.
#' @return A list of control settings for [ripr_init()].
#' @examples
#' ripr_control(n_seeds = 50L, snapshot = "step")
#' @export
ripr_control <- function(
  n_seeds = 200L,
  n_restarts = 25L,
  fc_tol = 1e-10,
  fc_max_iter = 500L,
  snapshot = c("none", "step", "all")
) {
  snapshot <- rlang::arg_match(snapshot)
  rlang::check_number_whole(n_seeds, min = 0, max = 2147483647)
  rlang::check_number_whole(n_restarts, min = 1, max = 2147483647)
  rlang::check_number_decimal(fc_tol, min = 0)
  rlang::check_number_whole(fc_max_iter, min = 1, max = 2147483647)
  list(
    n_seeds = as.integer(n_seeds),
    n_restarts = as.integer(n_restarts),
    fc_tol = fc_tol,
    fc_max_iter = as.integer(fc_max_iter),
    snapshot = snapshot
  )
}
