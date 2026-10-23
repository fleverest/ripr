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
#' @param weight_solver How a corrective weight solve (`correct = TRUE`,
#'   [weight_step()], or `ripr_finish(reoptimise = TRUE)`) finds the optimal
#'   weights for fixed atoms. `"sqp"` (default) is an exact-Hessian sequential
#'   quadratic programming method, cf. mix-SQP \insertCite{Kim2020}{ripr}, which
#'   stops once the Frank--Wolfe gap over the current support is below `kkt_tol`
#'   and falls back on `"slsqp"` if it doesn't converge. `"slsqp"` is nloptr's
#'   quasi-Newton SLSQP, which stops when the KL stops changing (much) in
#'   floating point, typically the gap is larger than the change in KL at that
#'   point.
#' @param kkt_tol For `weight_solver = "sqp"`: the target KKT residual, the
#'   largest of \eqn{|G_c - 1|}{|G_c - 1|} over atoms with weight and
#'   \eqn{G_c - 1}{G_c - 1} over atoms without.
#' @param fc_tol For `"slsqp"`, and the fallback: its relative tolerance on the
#'   weights.
#' @param fc_max_iter Cap on the iterations (`"sqp"`) or objective evaluations
#'   (`"slsqp"`) a solve may take.
#' @param merge_atol,merge_rtol When two atoms of the same part agree in every
#'   coordinate to within `merge_atol + merge_rtol * max(|a|, |b|)`, they are
#'   merged into one at their weight-averaged position: after every step that
#'   adds or moves atoms, and before [weight_step()]'s solve. Two copies of one
#'   atom make the weight problem degenerate (its Hessian singular), so this
#'   keeps the solve well-posed. Meant to catch atoms that coincide up to
#'   rounding; large values change the fit.
#' @return A list of control settings for [ripr_init()].
#' @references
#'   \insertAllCited{}
#' @examples
#' ripr_control(n_seeds = 50L, snapshot = "step")
#' @export
ripr_control <- function(
  n_seeds = 200L,
  n_restarts = 25L,
  weight_solver = c("sqp", "slsqp"),
  kkt_tol = 1e-12,
  fc_tol = 1e-10,
  fc_max_iter = 500L,
  merge_atol = 1e-9,
  merge_rtol = 1e-9,
  snapshot = c("none", "step", "all")
) {
  snapshot <- rlang::arg_match(snapshot)
  weight_solver <- rlang::arg_match(weight_solver)
  rlang::check_number_decimal(kkt_tol, min = 0)
  rlang::check_number_whole(n_seeds, min = 0, max = 2147483647)
  rlang::check_number_whole(n_restarts, min = 1, max = 2147483647)
  rlang::check_number_decimal(fc_tol, min = 0)
  rlang::check_number_whole(fc_max_iter, min = 1, max = 2147483647)
  rlang::check_number_decimal(merge_atol, min = 0)
  rlang::check_number_decimal(merge_rtol, min = 0)
  list(
    n_seeds = as.integer(n_seeds),
    n_restarts = as.integer(n_restarts),
    weight_solver = weight_solver,
    kkt_tol = kkt_tol,
    fc_tol = fc_tol,
    fc_max_iter = as.integer(fc_max_iter),
    merge_atol = merge_atol,
    merge_rtol = merge_rtol,
    snapshot = snapshot
  )
}
