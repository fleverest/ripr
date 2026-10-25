#' @include ripr-package.R
NULL

#' Tuning for the RIPr optimiser
#'
#' These settings often will not require adjustment on your part, they are
#' mainly provided here in one place for troubleshooting.
#'
#' The defaults that can be overridden by `ripr_control()` are related to
#' one of four tasks:
#'
#' @section Search:
#' The linear oracle, the Li--Barron oracle and the EM atom step each maximise
#' over the null by multi-start SLSQP in each part's chart.
#'
#' @section Slide:
#' [sfw_step()]'s joint SLSQP descent over every weight and atom position.
#'
#' @section Weight solves:
#' The corrective solve over the weights with fixed support (atom locations),
#' `correct_weights = TRUE`, [weight_step()], [sfw_step()] and
#' `ripr_finish(correct_weights = TRUE)` all run solves for the optimal weights.
#'
#' @section Merging:
#' Two copies of one atom make the weight problem degenerate (its Hessian
#' singular), so atoms of a part that coincide are merged into one at their
#' weight-averaged position.
#'
#' @param n_seeds Search: random chart seeds per search.
#' @param n_restarts Search: number of the best seeds refined by SLSQP per
#'   search.
#' @param search_xtol Search: SLSQP's relative tolerance on chart coordinates
#'   when refining a seed. Near a maximum the objective is flat, so its value
#'   is accurate to roughly the square of this.
#' @param search_max_eval Search: cap on objective evaluations per refined
#'   seed.
#' @param slide_xtol Slide: SLSQP's relative tolerance on the joint weights
#'   and chart coordinates.
#' @param slide_max_eval Slide: cap on objective evaluations per slide.
#' @param weight_solver Weight solve: `"sqp"` (default) is an exact-Hessian
#'   sequential quadratic programming method
#'   \insertCite{NocedalWright2006}{ripr}, cf. mix-SQP
#'   \insertCite{Kim2020}{ripr}, which stops once the Frank--Wolfe gap over the
#'   current support is below `weight_sqp_tol` and falls back on `"slsqp"` if
#'   it doesn't converge. `"slsqp"` is nloptr's quasi-Newton SLSQP, which stops
#'   once the KL stops changing in floating point, which can leave the gap over
#'   the support much larger than `weight_sqp_tol`.
#' @param weight_sqp_tol Weight solve, `"sqp"`: the target KKT residual, the
#'   largest of \eqn{|G_c - 1|}{|G_c - 1|} over atoms with weight and the
#'   positive part of \eqn{G_c - 1}{G_c - 1} over atoms without.
#' @param weight_slsqp_tol Weight solve, `"slsqp"` and the fallback: its
#'   relative tolerance on the weights.
#' @param weight_max_iter Weight solve: cap on the iterations (`"sqp"`) or
#'   objective evaluations (`"slsqp"`) a solve may take.
#' @param merge_atol,merge_rtol Merging: two atoms of the same part are the same
#'   point when every coordinate agrees to within
#'   `merge_atol + merge_rtol * max(|a|, |b|)`. Applied after every step that
#'   adds or moves atoms and before every weight solve. Meant to catch atoms
#'   that coincide up to rounding; large values change the fit.
#' @param candidate_merge_tol Merging: `fw_step(insert = "multi")` inserts
#'   every distinct local maximum the search found; restarts that reached the
#'   same maximum end within the search's precision of each other, so
#'   candidates within this (absolute and relative) tolerance are one maximum,
#'   and a candidate this close to a current atom is dropped.
#' @param snapshot How often to record the fitted mixture alongside the trace.
#'   `"none"` (default) never records the state; `"step"` records the atoms and
#'   weights once per step call (e.g. [fw_step()]), so `fw_step(times = 10)`
#'   yields one snapshot; `"all"` records once per iteration, yielding ten.
#'   The trace is always recorded; snapshots copy the whole mixture, so their
#'   memory cost grows with support size.
#' @return A list of control settings for [ripr_init()].
#' @references
#'   \insertAllCited{}
#' @examples
#' ripr_control(n_seeds = 50L, snapshot = "step")
#' @export
ripr_control <- function(
  n_seeds = 200L,
  n_restarts = 25L,
  search_xtol = 1e-8,
  search_max_eval = 1000L,
  slide_xtol = 1e-10,
  slide_max_eval = 2000L,
  weight_solver = c("sqp", "slsqp"),
  weight_sqp_tol = 1e-12,
  weight_slsqp_tol = 1e-10,
  weight_max_iter = 500L,
  merge_atol = 1e-9,
  merge_rtol = 1e-9,
  candidate_merge_tol = 1e-6,
  snapshot = c("none", "step", "all")
) {
  # Checks that report against ripr_control() and the argument's own name.
  whole <- function(x, min, arg = rlang::caller_arg(x)) {
    rlang::check_number_whole(
      x,
      min = min,
      max = 2147483647,
      arg = arg,
      call = rlang::caller_env()
    )
  }
  tol <- function(x, arg = rlang::caller_arg(x)) {
    rlang::check_number_decimal(
      x,
      min = 0,
      arg = arg,
      call = rlang::caller_env()
    )
  }
  snapshot <- rlang::arg_match(snapshot)
  weight_solver <- rlang::arg_match(weight_solver)
  whole(n_seeds, min = 0)
  whole(n_restarts, min = 1)
  tol(search_xtol)
  whole(search_max_eval, min = 1)
  tol(slide_xtol)
  whole(slide_max_eval, min = 1)
  tol(weight_sqp_tol)
  tol(weight_slsqp_tol)
  whole(weight_max_iter, min = 1)
  tol(merge_atol)
  tol(merge_rtol)
  tol(candidate_merge_tol)
  list(
    n_seeds = as.integer(n_seeds),
    n_restarts = as.integer(n_restarts),
    search_xtol = search_xtol,
    search_max_eval = as.integer(search_max_eval),
    slide_xtol = slide_xtol,
    slide_max_eval = as.integer(slide_max_eval),
    weight_solver = weight_solver,
    weight_sqp_tol = weight_sqp_tol,
    weight_slsqp_tol = weight_slsqp_tol,
    weight_max_iter = as.integer(weight_max_iter),
    merge_atol = merge_atol,
    merge_rtol = merge_rtol,
    candidate_merge_tol = candidate_merge_tol,
    snapshot = snapshot
  )
}
