#' @include fit.R
NULL

#' Fit the RIPr in one call
#'
#' Runs one of the recommended fitting methods until the Frank--Wolfe gap is
#' below `tol` (or `max_steps` steps have been taken), and returns the result.
#' It is a convenient shorthand for the [ripr_init()], step verbs (e.g.
#' [sfw_step()]) with `until = gap_below(tol)`, followed by
#' `ripr_finish(record_gap = TRUE)`; use those directly to combine methods or
#' to record more specific details.
#'
#' Each of the methods available here re-solve the weights after every step,
#' only differing in how the support itself evolves:
#' - `"sfw"`, sliding Frank--Wolfe ([sfw_step()]): adds the oracle's atom, then
#'   "slides" every atom and weight jointly. This typically converges very fast,
#'   both in fewer total iterations and the shortest CPU time.
#' - `"pdap"` (`fw_step(correct_weights = TRUE, insert = "multi")`): adds every
#'   distinct local maximiser the oracle found, then re-solves the weights,
#'   but the support never moves. Each iteration is cheaper than `"sfw"` steps,
#'   and can sometimes outperform it as a result.
#' - `"fcfw"`, fully-corrective Frank--Wolfe (`fw_step(correct_weights =
#'   TRUE)`): adds just the oracle's atom (the single largest maximiser). Atoms
#'   never move, so it converges relatively slowly.
#'
#' Passing a [ripr_state] or [ripr_fit] as `x` resumes the fit: iterations
#' pick up from the last state. `max_steps` counts the new iterations only. Name
#' the other arguments when resuming, since the second position is `null`.
#'
#' @param x The alternative \eqn{Q}{Q}, a [distribution]; or a [ripr_state] or
#'   [ripr_fit] to resume.
#' @param null A [null_model]. Required for a new fit, and must be `NULL` when
#'   resuming, since a state keeps its own.
#' @param method `"sfw"` (default), `"pdap"` or `"fcfw"`; see Details.
#' @param tol Stop once the Frank--Wolfe gap, \eqn{\sup_\theta G(\theta) -
#'   1}{sup G(theta) - 1}, is below this.
#' @param max_steps The most steps to take.
#' @param engine An engine spec, e.g. [exact_engine()]. New fits only.
#' @param control From [ripr_control()]. For a resumed fit, replaces the
#'   state's settings if given.
#' @return A [ripr_fit]. Its `gap_final` is the gap over the returned mixture,
#'   searched afresh; when it is not below `tol`, a warning of class
#'   `ripr_budget_warning` says so.
#' @seealso [ripr_init()] and the step verbs, e.g. [sfw_step()], for fitting
#'   step by step; [certify()] to bound the gap.
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
#' fit <- ripr(Q, plurality, tol = 1e-10)
#' fit
#'
#' # A fit that ran out of steps can be resumed where it stopped.
#' short <- suppressWarnings(ripr(Q, plurality, method = "fcfw", max_steps = 2L))
#' medium <- ripr(short, method = "fcfw", tol = 1e-6)
#' long <- ripr(medium, method = "sfw", tol = 1e-10)
#' @export
ripr <- function(
  x,
  null = NULL,
  method = c("sfw", "pdap", "fcfw"),
  tol = 1e-8,
  max_steps = 100L,
  engine = exact_engine(),
  control = ripr_control()
) {
  method <- rlang::arg_match(method)
  rlang::check_number_decimal(tol, min = 0)
  rlang::check_number_whole(max_steps, min = 1, max = 2147483647)

  if (S7_inherits(x, ripr_fit)) {
    x <- x@state
  }
  if (S7_inherits(x, ripr_state)) {
    if (is.character(null)) {
      stop(
        "a fit being resumed takes no `null`; name the method, as in ",
        "`ripr(fit, method = \"",
        null[[1L]],
        "\")`.",
        call. = FALSE
      )
    }
    if (!is.null(null) || !missing(engine)) {
      stop(
        "a fit being resumed keeps its own `null` and `engine`; leave them ",
        "unset.",
        call. = FALSE
      )
    }
    state <- x
    if (!missing(control)) {
      state@control <- control
    }
  } else if (S7_inherits(x, distribution)) {
    if (is.null(null)) {
      stop("`null` is required to begin a fit.", call. = FALSE)
    }
    state <- ripr_init(x, null, engine, control = control)
  } else {
    stop(
      "`x` must be the alternative (a distribution), or a `ripr_state` or ",
      "`ripr_fit` to resume.",
      call. = FALSE
    )
  }

  rows_before <- nrow(state@trace)
  until <- gap_below(tol)
  state <- switch(
    method,
    sfw = sfw_step(state, max_steps, until = until),
    pdap = fw_step(
      state,
      max_steps,
      correct_weights = TRUE,
      insert = "multi",
      until = until
    ),
    fcfw = fw_step(state, max_steps, correct_weights = TRUE, until = until)
  )
  fit <- ripr_finish(state, record_gap = TRUE)

  if (!(fit@gap_final < tol)) {
    steps <- nrow(state@trace) - rows_before
    rlang::warn(
      c(
        sprintf(
          "The gap is %s after %s, not below `tol` = %s.",
          format(signif(fit@gap_final, 3L)),
          count_label(steps, "step"),
          tol
        ),
        i = "Resume with `ripr(fit, max_steps = ...)` to take more."
      ),
      class = "ripr_budget_warning"
    )
  }
  fit
}
