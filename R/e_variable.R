#' RIPr e-variable
#'
#' Bundles a numerator distribution `Q` and its reverse information projection
#' `P*` into a single object that evaluates the e-variable on data. The rescaled
#' e-variable `(Q / P*) / (1 + gap)` has expectation at most 1 under every null
#' distribution; `gap` is the certified duality gap the projection achieved
#' (`0` leaves the raw likelihood ratio uncorrected).
#'
#' Usually obtained by fitting with [run_ripr()] (as the `$e_variable` element of
#' its result), but can also be assembled directly from the pieces.
#'
#' @param numerator The numerator `Q` as an [outcome_distribution].
#' @param projection The reverse information projection `P*` as an
#'   [outcome_distribution] (typically a [marginal] whose mixing measure is a
#'   [finite_mixing]).
#' @param gap The certified duality gap used to rescale. Default `0`. Floored at
#'   `0`, since a negative gap is numerical noise and would otherwise inflate the
#'   e-value.
#' @param gap_certified Is `gap` a *proven* upper bound on `sup_theta G - 1`
#'   over the null? Default `FALSE`, which is the honest default: it means the
#'   rescaling rests on a heuristic oracle ([oracle()], multi-start BFGS on a
#'   generally non-convex objective) whose error is one-sided in the **unsafe**
#'   direction -- a missed global optimum makes `gap` too small and the
#'   e-variable under-corrected, so its null expectation may exceed 1. `TRUE`
#'   means the gap came from [oracle_bound()] and the expectation bound holds as
#'   an inequality, not as evidence. Take it from `certify()`'s `gap_certified`;
#'   do not set it by hand.
#' @return An `e_variable`.
#' @export
e_variable <- new_class(
  "e_variable",
  properties = list(
    numerator = outcome_distribution,
    projection = outcome_distribution,
    gap = class_numeric,
    gap_certified = class_logical
  ),
  constructor = function(numerator, projection, gap = 0, gap_certified = FALSE) {
    gap <- as.numeric(gap)
    gap_certified <- isTRUE(gap_certified)
    # Under certification the floor below is provably inert: the mixture weights
    # satisfy sum_j w_j G(theta_j) = 1, so sup G >= 1 and the gap is >= 0 as an
    # inequality rather than as an approximation. A negative one therefore
    # indicates a bug in the bound, not noise, and must not pass quietly.
    if (gap_certified && gap < 0) {
      warning(sprintf(
        paste0(
          "certified gap is negative (%g). `sum_j w_j G(theta_j) = 1` forces ",
          "`sup_theta G >= 1`, so a certified gap cannot be below 0: this is ",
          "a bug in the bound, not numerical noise. Flooring it at 0 anyway."
        ),
        gap
      ), call. = FALSE)
    }
    # A duality gap is non-negative by definition; a slightly negative value is
    # numerical noise near convergence. Floor it at 0 so the correction never
    # divides by `1 + gap < 1` (which would inflate, not deflate, the e-value).
    new_object(
      S7_object(),
      numerator = numerator,
      projection = projection,
      gap = max(gap, 0),
      gap_certified = gap_certified
    )
  }
)

#' Evaluate an e-variable on outcomes
#'
#' Computes the e-value `e(x) = (Q(x) / P*(x)) / (1 + gap)` at each outcome `x`.
#'
#' @param ev An `e_variable`.
#' @param x Outcomes: an `(N, d)` matrix, or a length-d vector for a single
#'   outcome.
#' @param log Return the log e-value instead? Default `FALSE`.
#' @param corrected Apply the `1 / (1 + gap)` rescaling that makes the e-variable
#'   valid under the null? Default `TRUE`. With `FALSE` the raw likelihood ratio
#'   `Q(x) / P*(x)` (the e-variable before the duality-gap correction) is
#'   returned.
#' @return Numeric vector of (log) e-values, one per outcome.
#' @export
e_value <- new_generic(
  "e_value",
  "ev",
  function(ev, x, log = FALSE, corrected = TRUE) {
    S7::S7_dispatch()
  }
)

method(e_value, e_variable) <- function(ev, x, log = FALSE, corrected = TRUE) {
  log_e <- dist_log_density(ev@numerator, x) -
    dist_log_density(ev@projection, x)
  if (corrected) {
    log_e <- log_e - log1p(ev@gap)
  }
  if (log) log_e else exp(log_e)
}

method(print, e_variable) <- function(x, ...) {
  # An outcome_distribution is described by its mixing measure when it is a
  # marginal; otherwise (a direct outcome law) by its own class.
  describe <- function(d) {
    if (S7_inherits(d, marginal)) S7_class(d@mixing)@name else S7_class(d)@name
  }
  proj <- x@projection
  n_atoms <- if (
    S7_inherits(proj, marginal) && S7_inherits(proj@mixing, finite_mixing)
  ) {
    ncol(proj@mixing@components)
  } else {
    NA_integer_
  }
  cat("<e_variable>  e(x) = (Q / P*) / (1 + gap)\n")
  cat(sprintf("  numerator  Q  : %s\n", describe(x@numerator)))
  cat(sprintf(
    "  projection P* : %s%s\n",
    describe(proj),
    if (is.na(n_atoms)) "" else sprintf(" (%d atom%s)", n_atoms, if (n_atoms == 1L) "" else "s")
  ))
  cat(sprintf("  gap           : %.4g   (correction 1 + gap = %.6g)\n", x@gap, 1 + x@gap))
  cat(sprintf(
    "  gap certified : %s\n",
    if (x@gap_certified) "TRUE" else "FALSE  (heuristic oracle; see ?certify)"
  ))
  invisible(x)
}
